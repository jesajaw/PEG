/*
Copyright (C) 2026 Jesaja Weintritt (jesaja.weintritt@stud.eah-jena.de) and 2012 Mark Boots (mark.boots@usask.ca).

This program was originally implemented as a part of the Parallel Efficiency of Gratings project PEG and got reworked in 2026. PEG is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License, version 3, as published by the Free Software Foundation.
See <http://www.gnu.org/licenses/> for details.

This reworked version contains substantial modifications by Jesaja Weintritt (2026) and has not been independently verified against the original. It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; use at your own risk and verify results independently.
*/

#include "TMSolver.h"

#include <cmath>
#include <algorithm>
#include <limits>
#include <iostream>

#include <omp.h>
#include <boost/numeric/odeint.hpp>

namespace {
	/// Thrown from TMSolver::odeFunction() when the grating expansion cannot be computed at
	/// some y. Caught in integrateTrialSolutionAlongY() and mapped to Result::InvalidGratingFailure.
	/// (Own copy, mirroring TESolver.cpp's anonymous-namespace struct -- duplication is intentional
	/// for now, see the pending SolverSupport extraction.)
	struct GratingExpansionError {};
}

TMSolver::TMSolver(const Grating& grating, const MathOptions& mo, int numThreads, bool measureTiming)
	: g_(grating)
{
	numThreads_ = numThreads;
	measureTiming_ = measureTiming;
	time_ = omp_get_wtime();

	N_ = mo.N;
	integrationTolerance_ = mo.integrationTolerance;
	twoNp1_ = 2*N_ + 1;
	fourNp2_ = 4*N_ + 2;
	eightNp4_ = 8*N_ + 4;

	wVectors_.assign(fourNp2_, std::vector<double>(eightNp4_));

	T11_.resize(twoNp1_, twoNp1_);
	T12_.resize(twoNp1_, twoNp1_);
	T21_.resize(twoNp1_, twoNp1_);
	T22_.resize(twoNp1_, twoNp1_);
	S12_.resize(twoNp1_, twoNp1_);
	S22_.resize(twoNp1_, twoNp1_);
	Zinv_.resize(twoNp1_, twoNp1_);
	Z_.resize(twoNp1_, twoNp1_);
	workMatrix_.resize(twoNp1_, twoNp1_);

	alpha_.resize(twoNp1_);
	betaM_.resize(twoNp1_);
	beta1_.resize(twoNp1_);
	BM_.resize(twoNp1_);

	// one k2_ array per thread. NOTE: sized twoNp1_ (2N+1), matching TESolver's actual (not its
	// docstring's claimed) buffer size -- see accompanying note on the banded |n-m|<=N coupling.
	k2_.assign(numThreads_, std::vector<std::complex<double>>(twoNp1_));

	// one (2N+1 x 2N+1) Toeplitz(k^2) workspace per thread, reused across odeFunction() calls.
	toeplitz_.assign(numThreads_, Eigen::MatrixXcd(twoNp1_, twoNp1_));

	timing_[0] = omp_get_wtime() - time_;
}

Result TMSolver::getEffTM(double incidenceDeg, double wl, double rmsRoughnessNm, bool printDebugOutput) {

	time_ = omp_get_wtime();

	// 1. Setup incidence variables and constants
	/////////////////////////////////////////

	wl_ = wl;
	v_1_ = g_.substrateRefractiveIndex(wl_);
	if(v_1_.real() == 0.0 && v_1_.imag() == 0.0) {
		return Result::MissingRefractiveDataFailure;
	}
	if(g_.coatingThickness() > 0) {
		v_c_ = g_.coatingRefractiveIndex(wl_);
		if(v_c_.real() == 0.0 && v_c_.imag() == 0.0) {
			return Result::MissingRefractiveDataFailure;
		}
	}

	timing_[1] = time_;
	time_ = omp_get_wtime();
	timing_[1] = time_ - timing_[1];

	// 2. compute all alpha_n and beta1_n, betaM_n.
	///////////////////////////////////

	computeAlphaAndBeta(incidenceDeg);
	computeLayers();

	timing_[2] = time_;
	time_ = omp_get_wtime();
	timing_[2] = time_ - timing_[2];

	if(printDebugOutput) {
		std::cout << "\nWavelength wl (um): " << wl_ << std::endl;
		std::cout << "Refractive index: " << v_1_.real() << ", " << v_1_.imag() << std::endl;
		std::cout << "Grating height a (um): " << g_.totalHeight() << std::endl;
		std::cout << "Grating period (um): " << g_.period() << std::endl;
		std::cout << "Number of layers: " << numLayers_ << std::endl;
		for(int m=1; m<M_; ++m) {
			std::cout << "   y_" << m << " = " << y_[m] << std::endl;
		}

		std::cout << "\nbeta1_n:" << std::endl;
		for(int i=0; i<twoNp1_; ++i) {
			std::cout << i - N_ << ":\t" << beta1_[i].real() << "\t\t" << beta1_[i].imag() << std::endl;
		}

		std::cout << "\nbetaM_n:" << std::endl;
		for(int i=0; i<twoNp1_; ++i) {
			std::cout << i - N_ << ":\t" << betaM_[i].real() << "\t\t" << betaM_[i].imag() << std::endl;
		}
	}

	// 3. Recursive computation of S-matrix below each layer.
	/////////////////////////////////////////////////////////////

	Result::Code status = computeTMatrixBelowLayer(2, printDebugOutput);
	if(status != Result::Success)
		return status;

	timing_[3] = time_;
	time_ = omp_get_wtime();
	timing_[3] = time_ - timing_[3];

	Zinv_ = T11_;
	{
		Eigen::PartialPivLU<Eigen::MatrixXcd> lu(Zinv_);
		S22_ = lu.inverse();
	}
	S12_ = T21_ * S22_;

	timing_[4] = time_;
	time_ = omp_get_wtime();
	timing_[4] = time_ - timing_[4];

	for(int m=3; m<M_; ++m) {

		time_ = omp_get_wtime();

		status = computeTMatrixBelowLayer(m, printDebugOutput);
		if(status != Result::Success) return status;

		timing_[3] += omp_get_wtime() - time_;

		time_ = omp_get_wtime();

		Zinv_ = T11_ + T12_ * S12_;

		{
			Eigen::PartialPivLU<Eigen::MatrixXcd> lu(Zinv_);
			Z_ = lu.inverse();
		}

		workMatrix_ = T21_ + T22_ * S12_;
		S12_ = workMatrix_ * Z_;
		S22_ = S22_ * Z_;

		timing_[4] += omp_get_wtime() - time_;
	}

	// 4. Calculate B_n^M from center column of S matrix * exp(...).
	//////////////////////////////////////////////
	time_ = omp_get_wtime();

	computeBMFromSMatrix();

	timing_[5] = time_;
	time_ = omp_get_wtime();
	timing_[5] = time_ - timing_[5];

	if(printDebugOutput) {
		std::cout << "\nBM_:" << std::endl;
		for(int i=0; i<twoNp1_; ++i) {
			std::cout << i - N_ << ":\t" << BM_[i].real() << "\t\t" << BM_[i].imag() << std::endl;
		}
	}

	// 5. Compute efficiency and put into result structure.
	////////////////////////////////////////
	// NOTE: unlike a naive port of the TM efficiency formula (eta_n = |B_n|^2 * (Re(beta_n)/eps_r,n)
	// / (Re(beta_0)/eps_r,0), TMSolver_alg_ref.md \S5.2), the eps_r ratio here is exactly 1: every
	// reflected order (including the reference order N_) lives in the same vacuum superstrate
	// (eps_r=1), so this reduces to literally the same formula as TESolver::getEffTE().

	Result result(N_);
	result.wavelength = wl_;
	result.incidenceDeg = incidenceDeg;

	double effSum = 0;
	for(int i=0; i<twoNp1_; ++i) {
		result.eff[i] = std::norm(BM_[i]) * betaM_[i].real() / betaM_[N_].real();
		effSum += result.eff[i];
	}

	timing_[6] = time_;
	time_ = omp_get_wtime();
	timing_[6] = time_ - timing_[6];

	if(measureTiming_) {
		std::cout << "Timing Profile:" << std::endl;
		std::cout << "   Allocate Memory: " << timing_[0] << std::endl;
		std::cout << "   Look up refractive index: " << timing_[1] << std::endl;
		std::cout << "   Compute alpha, beta values and layers: " << timing_[2] << std::endl;
		std::cout << "   Numerically integrating trial solutions: " << timing_[3] << std::endl;
		std::cout << "   Matrix operations: " << timing_[4] << std::endl;
		std::cout << "   Computing Rayleigh coeffients B_n: " << timing_[5] << std::endl;
		std::cout << "   Compute and package efficiencies: " << timing_[6] << std::endl;
		time_ = timing_[0] + timing_[1] + timing_[2] + timing_[3] + timing_[4] + timing_[5] + timing_[6];
		std::cout << "   Total (solver) time: " << time_ << std::endl << std::endl;
	}

	if(printDebugOutput) {
		std::cout << "Sum of reflected efficiencies: " << effSum << std::endl;

		// Debug-only transmission check. Unlike reflection, transmitted order i (substrate,
		// eps_r = v_1_^2) and the reference order N_ (vacuum, eps_r = 1) do NOT share the same
		// medium, so the eps_r ratio from TMSolver_alg_ref.md \S5.2/\S7.7 does not cancel here.
		double a = g_.totalHeight();
		std::vector<std::complex<double>> A1(twoNp1_);
		std::vector<double> e_t(twoNp1_);
		double sumTransmitted = 0;
		std::complex<double> epsSub = v_1_ * v_1_;
		for(int i=0; i<twoNp1_; i++) {
			A1[i] = S22_(i, N_) * std::exp(std::complex<double>(0, -a) * betaM_[N_]);
			std::complex<double> ratio_i = beta1_[i] / epsSub;
			e_t[i] = std::norm(A1[i]) * ratio_i.real() / betaM_[N_].real();
			sumTransmitted += e_t[i];
		}
		std::cout << "Sum of transmitted efficiencies: " << sumTransmitted << std::endl;
		std::cout << "Total efficiency (should be <= 1): " << sumTransmitted + effSum << std::endl;
	}

	if(rmsRoughnessNm > 0) {
		double roughness = g_.roughnessFactor(rmsRoughnessNm/1000., wl, g_.coatingThickness() > 0 ? v_c_ : v_1_, incidenceDeg);
		for(int i=0; i<twoNp1_; ++i)
			result.eff[i] = roughness*result.eff.at(i);
	}

	return result;
}

std::complex<double> TMSolver::complex_sqrt_upperComplexPlane(std::complex<double> z) {
	std::complex<double> w = std::sqrt(z);
	if(w.imag() < 0)
		w = -w;
	return w;
}

double TMSolver::conditionNumber(const Eigen::MatrixXcd& A) {
	Eigen::JacobiSVD<Eigen::MatrixXcd> svd(A);
	const auto& sv = svd.singularValues();
	double maxSv = sv(0);
	double minSv = sv(sv.size()-1);
	if(minSv == 0)
		return std::numeric_limits<double>::infinity();
	return maxSv / minSv;
}

void TMSolver::computeAlphaAndBeta(double incidenceDeg)
{
	double theta_2 = incidenceDeg * M_PI / 180;

	double k_2 = 2 * M_PI / wl_;
	std::complex<double> k_1 = v_1_ * k_2;

	double d = g_.period();

#pragma omp parallel for num_threads(numThreads_)
	for(int i=0; i<twoNp1_; i++) {
		int n = i - N_;

		double alpha = k_2 * sin(theta_2) + 2 * M_PI * n / d;
		alpha_[i] = alpha;

		double k22minusAn2 = k_2*k_2 - alpha*alpha;
		if(k22minusAn2 >= 0)
			betaM_[i] = std::complex<double>(sqrt(k22minusAn2), 0);
		else
			betaM_[i] = std::complex<double>(0, sqrt(-k22minusAn2));

		std::complex<double> k12minusAn2 = k_1*k_1 - alpha*alpha;
		beta1_[i] = complex_sqrt_upperComplexPlane(k12minusAn2);
	}
}

void TMSolver::computeLayers()
{
	double a = g_.totalHeight();
	double magicNumber = 3;

	numLayers_ = std::max( std::abs(betaM_[0])*a/magicNumber, std::abs(betaM_[2*N_])*a/magicNumber );
	if(numLayers_ < 1)
		numLayers_ = 1;

	M_ = numLayers_+2;

	y_.resize(M_);

	for(int m=1; m<M_; ++m) {
		y_[m] = double(m-1)/numLayers_*a;
	}
}

Result::Code TMSolver::computeGratingExpansion(double y, std::complex<double>* k2) const {

	std::complex<double> k_M(2 * M_PI / wl_, 0);
	std::complex<double> k_1 = v_1_ * k_M;
	std::complex<double> k_c = v_c_ * k_M;

	std::complex<double> k2_M = k_M * k_M;
	std::complex<double> k2_1 = k_1 * k_1;
	std::complex<double> k2_c = k_c * k_c;

	double stepsX[PEG_MAX_PROFILE_CROSSINGS];
	std::complex<double> stepsK2[PEG_MAX_PROFILE_CROSSINGS];

	int numSteps = g_.computeK2StepsAtY(y, k2_M, k2_1, k2_c, stepsX, stepsK2);
	if(numSteps < 1)
		return Result::InvalidGratingFailure;

	computeGratingExpansion(stepsX, stepsK2, numSteps, k2);

	return Result::Success;
}

void TMSolver::computeGratingExpansion(const double *stepsX, const std::complex<double> *stepsK2, int numSteps, std::complex<double> *k2) const
{
	double d = g_.period();
	double K = 2*M_PI/d;

	if(numSteps == 1) {
		for(int i=0; i<twoNp1_; ++i)
			k2[i] = std::complex<double>(0,0);
		k2[N_] = stepsK2[0];
		return;
	}

	std::complex<double> sigma[PEG_MAX_PROFILE_CROSSINGS];
	for(int p=0;p<numSteps-1; ++p)
		sigma[p] = stepsK2[p+1] - stepsK2[p];
	sigma[numSteps-1] = stepsK2[0] - stepsK2[numSteps-1];

	for(int i=0; i<twoNp1_; ++i) {
		int n = i - N_;

		if(n == 0) {
			std::complex<double> f0 = stepsK2[0] * d;
			for(int p=0; p<numSteps; ++p)
				f0 -= sigma[p] * stepsX[p];
			k2[i] = f0 / d;
		}
		else {
			std::complex<double> fn(0,0);
			for(int p=0; p<numSteps; ++p) {
				double nKx = n*K*stepsX[p];
				fn += sigma[p] * std::complex<double>(sin(nKx), cos(nKx));
			}
			k2[i] = fn / (-2*M_PI*n);
		}
	}
}

void TMSolver::buildToeplitz(const std::complex<double>* coeffs, int N, Eigen::MatrixXcd& toeplitz)
{
	// coeffs[p + N] for p in [-N, N] (twoNp1_-sized buffer). Entries with |n-m| > N are set to
	// zero -- replicating the same banded-coupling restriction already present, undocumented, in
	// TESolver::odeFunction (see accompanying note).
	int twoNp1 = 2*N + 1;
	for(int i=0; i<twoNp1; ++i) {
		int n = i - N;
		for(int j=0; j<twoNp1; ++j) {
			int m = j - N;
			int p = n - m;
			toeplitz(i, j) = (p >= -N && p <= N) ? coeffs[p + N] : std::complex<double>(0,0);
		}
	}
}

std::complex<double>* TMSolver::k2ForCurrentThread() {
	return k2_[omp_get_thread_num()].data();
}

Eigen::MatrixXcd& TMSolver::toeplitzForCurrentThread() {
	return toeplitz_[omp_get_thread_num()];
}

void TMSolver::odeFunction(double y, const std::vector<double>& w_arr, std::vector<double>& f)
{
	std::complex<double>* localK2 = k2ForCurrentThread();
	if(computeGratingExpansion(y, localK2) != Result::Success) {
		std::cout << "ODE: Function Error: Cannot compute grating expansion at y = " << y << std::endl;
		throw GratingExpansionError{};
	}

	Eigen::MatrixXcd& toeplitzK2 = toeplitzForCurrentThread();
	buildToeplitz(localK2, N_, toeplitzK2);

	Eigen::VectorXcd u_vec(twoNp1_), w_vec(twoNp1_), alphaU(twoNp1_);
	for(int i=0; i<twoNp1_; ++i) {
		u_vec(i) = std::complex<double>(w_arr[2*i], w_arr[2*i+1]);
		w_vec(i) = std::complex<double>(w_arr[fourNp2_ + 2*i], w_arr[fourNp2_ + 2*i + 1]);
		alphaU(i) = alpha_[i] * u_vec(i);
	}

	double k0 = 2 * M_PI / wl_;
	double k0sq = k0 * k0;

	// du/dy = (1/k0^2) * Toeplitz(k^2) * w  -- plain multiply.
	Eigen::VectorXcd duDy = (toeplitzK2 * w_vec) / k0sq;

	// dw/dy = k0^2 * (D_alpha * v - u), where Toeplitz(k^2) * v = D_alpha * u  -- one LU solve.
	Eigen::PartialPivLU<Eigen::MatrixXcd> lu(toeplitzK2);
	Eigen::VectorXcd v_vec = lu.solve(alphaU);

	for(int i=0; i<twoNp1_; ++i) {
		std::complex<double> dwDy_i = k0sq * (alpha_[i] * v_vec(i) - u_vec(i));

		f[2*i] = duDy(i).real();
		f[2*i+1] = duDy(i).imag();
		f[fourNp2_ + 2*i] = dwDy_i.real();
		f[fourNp2_ + 2*i + 1] = dwDy_i.imag();
	}
}

Result::Code TMSolver::integrateTrialSolutionAlongY(std::vector<double>& w_arr, double yStart, double yEnd) {

	namespace odeint = boost::numeric::odeint;

	odeint::bulirsch_stoer<std::vector<double>> stepper(integrationTolerance_, integrationTolerance_);

	auto system = [this](const std::vector<double>& state, std::vector<double>& dwdy, double y) {
		odeFunction(y, state, dwdy);
	};

	double hStart = (yEnd - yStart)/200;

	try {
		odeint::integrate_adaptive(stepper, system, w_arr, yStart, yEnd, hStart);
	}
	catch(const GratingExpansionError&) {
		return Result::InvalidGratingFailure;
	}
	catch(...) {
		std::cout << "ODE: Integration failure between y = " << yStart << " and y = " << yEnd << std::endl;
		return Result::ConvergenceFailure;
	}

	return Result::Success;
}

void TMSolver::setIntegrationStartingValues(std::vector<double>& w_arr, int j, int m)
{
	std::fill(w_arr.begin(), w_arr.end(), 0.0);

	bool secondRound = false;
	if(j >= twoNp1_) {
		secondRound = true;
		j -= twoNp1_;
	}

	w_arr[2*j] = 1.0;

	// w = eps_r^{-1} * u'. At m==1 (substrate boundary) eps_r = v_1_^2; at m>1 (vacuum-referenced
	// intermediate boundary, matching TESolver's convention of using betaM_ there) eps_r = 1.
	std::complex<double> wStart = (m == 1 ? beta1_[j] : betaM_[j]) * std::complex<double>(0, secondRound ? 1 : -1);
	if(m == 1)
		wStart /= (v_1_ * v_1_);

	w_arr[fourNp2_ + 2*j] = wStart.real();
	w_arr[fourNp2_ + 2*j+1] = wStart.imag();
}

void TMSolver::computeBMFromSMatrix()
{
	double a = g_.totalHeight();

#pragma omp parallel for num_threads(numThreads_)
	for(int i=0; i<twoNp1_; i++) {
		BM_[i] = S12_(i, N_) * std::exp(std::complex<double>(0,-a) * (betaM_[i] + betaM_[N_]));
	}
}

Result::Code TMSolver::computeTMatrixBelowLayer(int m, bool printDebugOutput)
{
	bool integrationFailureOccurred = false;

#pragma omp parallel for num_threads(numThreads_) schedule(dynamic) reduction(||:integrationFailureOccurred)
	for(int j=0; j<fourNp2_; ++j) {

		std::vector<double>& w_arr = wVectorForP(j);

		setIntegrationStartingValues(w_arr, j, m-1);

		if(printDebugOutput && omp_get_thread_num() == 0) {
			std::cout << "Initial value u_{p=" << j-N_ << "}(yStart):" <<std::endl;
			std::cout << "     ";
			for(int n=0; n<twoNp1_; ++n)
				std::cout << w_arr[2*n] << "," << w_arr[2*n+1] << "    ";
			std::cout << std::endl;
			std::cout << "Initial value w_{p=" << j-N_ << "}(yStart):" <<std::endl;
			std::cout << "     ";
			for(int n=0; n<twoNp1_; ++n)
				std::cout << w_arr[fourNp2_ + 2*n] << "," << w_arr[fourNp2_ + 2*n+1] << "    ";
			std::cout << std::endl;
			std::cout << std::endl;
		}

		Result::Code status = integrateTrialSolutionAlongY(w_arr, y_[m-1], y_[m]);

		if(printDebugOutput && omp_get_thread_num() == 0) {
			std::cout << "Final value u_{p=" << j-N_ << "}(yEnd):" <<std::endl;
			std::cout << "     ";
			for(int n=0; n<twoNp1_; ++n)
				std::cout << w_arr[2*n] << "," << w_arr[2*n+1] << "    ";
			std::cout << std::endl;
			std::cout << "Final value w_{p=" << j-N_ << "}(yEnd):" <<std::endl;
			std::cout << "     ";
			for(int n=0; n<twoNp1_; ++n)
				std::cout << w_arr[fourNp2_ + 2*n] << "," << w_arr[fourNp2_ + 2*n+1] << "    ";
			std::cout << std::endl;
			std::cout << std::endl;
		}

		if(status != Result::Success)
			integrationFailureOccurred = true;
		else {
			// Transform formula is unchanged from TESolver's: this construction is always
			// referenced against betaM_ (vacuum, eps_r=1), so the eps_r factor from
			// TMSolver_alg_ref.md \S3.2 evaluates to 1 here -- see class-level comment in TMSolver.h.
			if(j >= twoNp1_) {
				int jj = j - twoNp1_;

				for(int i=0; i<twoNp1_; ++i) {
					std::complex<double> u_ij(w_arr[2*i], w_arr[2*i+1]);
					std::complex<double> w_ij(w_arr[fourNp2_ + 2*i], w_arr[fourNp2_ + 2*i + 1]);
					std::complex<double> temp = w_ij / (betaM_[i] * std::complex<double>(0,1));

					T12_(i, jj) = 0.5*(u_ij - temp);
					T22_(i, jj) = 0.5*(u_ij + temp);
				}
			}
			else {
				for(int i=0; i<twoNp1_; ++i) {
					std::complex<double> u_ij(w_arr[2*i], w_arr[2*i+1]);
					std::complex<double> w_ij(w_arr[fourNp2_ + 2*i], w_arr[fourNp2_ + 2*i + 1]);
					std::complex<double> temp = w_ij / (betaM_[i] * std::complex<double>(0,1));

					T11_(i, j) = 0.5*(u_ij - temp);
					T21_(i, j) = 0.5*(u_ij + temp);
				}
			}
		}
	}

	if(integrationFailureOccurred)
		return Result::ConvergenceFailure;
	else
		return Result::Success;
}