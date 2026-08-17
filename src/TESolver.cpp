/*
Copyright (C) 2026 Jesaja Weintritt (jesaja.weintritt@stud.eah-jena.de) and 2012 Mark Boots (mark.boots@usask.ca).

This program was originally implemented as a part of the Parallel Efficiency of Gratings project PEG and got reworked in 2026. PEG is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License, version 3, as published by the Free Software Foundation.
See <http://www.gnu.org/licenses/> for details.

This reworked version contains substantial modifications by Jesaja Weintritt (2026) and has not been independently verified against the original. It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; use at your own risk and verify results independently.
*/

#include "TESolver.h"

#include <cmath>
#include <algorithm>
#include <limits>
#include <iostream>

#include <omp.h>
#include <boost/numeric/odeint.hpp>

namespace {
	/// Thrown from TESolver::odeFunction() when the grating expansion cannot be computed at
	/// some y (invalid geometry, or y above the profile height). Caught in
	/// integrateTrialSolutionAlongY() and mapped to Result::InvalidGratingFailure.
	struct GratingExpansionError {};
}

TESolver::TESolver(const Grating& grating, const MathOptions& mo, int numThreads, bool measureTiming)
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

	// All storage below is RAII (std::vector / Eigen); nothing needs manual freeing anymore.
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

	// one k2_ array per thread, since they will be used simultaneously.
	k2_.assign(numThreads_, std::vector<std::complex<double>>(twoNp1_));

	timing_[0] = omp_get_wtime() - time_;		// time to allocate memory.
}

Result TESolver::getEffTE(double incidenceDeg, double wl, double rmsRoughnessNm, bool printDebugOutput) {

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
	timing_[1] = time_ - timing_[1];	// time to look up refractive index.

	// 2. compute all alpha_n and beta1_n, betaM_n.
	///////////////////////////////////

	computeAlphaAndBeta(incidenceDeg);	// (This uses a parallel loop).
	computeLayers();

	timing_[2] = time_;
	time_ = omp_get_wtime();
	timing_[2] = time_ - timing_[2];	// time to calculate alpha, beta, and layers.

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

	// Handle first layer separately, as a special case.
	Result::Code status = computeTMatrixBelowLayer(2, printDebugOutput);
	if(status != Result::Success)
		return status;

	timing_[3] = time_;
	time_ = omp_get_wtime();
	timing_[3] = time_ - timing_[3];	// timing_[3]: numerical integration.

	// For the first layer, we have Zinv_ = T11_.
	// S12_ = T21_ Zinv_^{-1}
	// S22_ = Zinv_^{-1}
	///////////////
	Zinv_ = T11_;
	{
		Eigen::PartialPivLU<Eigen::MatrixXcd> lu(Zinv_);
		S22_ = lu.inverse();	// S22_ = Zinv_^{-1}, ie, Z.
	}
	S12_ = T21_ * S22_;

	timing_[4] = time_;
	time_ = omp_get_wtime();
	timing_[4] = time_ - timing_[4];	// timing_[4]: matrix calcs.


	// First layer done. Handle subsequent layers
	for(int m=3; m<M_; ++m) {

		time_ = omp_get_wtime();

		status = computeTMatrixBelowLayer(m, printDebugOutput);
		if(status != Result::Success) return status;

		timing_[3] += omp_get_wtime() - time_;

		time_ = omp_get_wtime();

		// Zinv_ = T11_ + T12_ S12_.
		Zinv_ = T11_ + T12_ * S12_;

		// Invert Zinv_...
		{
			Eigen::PartialPivLU<Eigen::MatrixXcd> lu(Zinv_);
			Z_ = lu.inverse();
		}

		// S12 = (T21 + T22 S12) Z
		workMatrix_ = T21_ + T22_ * S12_;
		S12_ = workMatrix_ * Z_;
		// S22 = S22 Z
		S22_ = S22_ * Z_;

		timing_[4] += omp_get_wtime() - time_;
	}

	// 4.  Calculate B_n^M from center column of S matrix * exp(...).
	//////////////////////////////////////////////
	time_ = omp_get_wtime();

	computeBMFromSMatrix();

	timing_[5] = time_;
	time_ = omp_get_wtime();
	timing_[5] = time_ - timing_[5];	// timing_[5]: compute BM_ from S-matrix.

	if(printDebugOutput) {
		std::cout << "\nBM_:" << std::endl;
		for(int i=0; i<twoNp1_; ++i) {
			std::cout << i - N_ << ":\t" << BM_[i].real() << "\t\t" << BM_[i].imag() << std::endl;
		}
	}

	// 8. Now we have BM_. Compute efficiency and put into result structure.
	////////////////////////////////////////

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
	timing_[6] = time_ - timing_[6];	// timing_[6]: calculate efficiencies from BM_

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

		// in debugging, let's also compute and sum the transmitted efficiencies:
		double a = g_.totalHeight();
		std::vector<std::complex<double>> A1(twoNp1_);
		std::vector<double> e_t(twoNp1_);
		double sumTransmitted = 0;
		for(int i=0; i<twoNp1_; i++) {
			A1[i] = S22_(i, N_) * std::exp(std::complex<double>(0, -a) * betaM_[N_]);
			e_t[i] = std::norm(A1[i]) * beta1_[i].real() / betaM_[N_].real();
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

std::complex<double> TESolver::complex_sqrt_upperComplexPlane(std::complex<double> z) {

	std::complex<double> w = std::sqrt(z); // returns w in the right half of the complex plane.

	if(w.imag() < 0) {
		// flip to the upper half-plane: e^{i*pi} = -1, same square, opposite branch.
		w = -w;
	}

	return w;
}

Result::Code TESolver::computeGratingExpansion(double y, std::complex<double>* k2) const {

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

Result::Code TESolver::integrateTrialSolutionAlongY(std::vector<double>& w, double yStart, double yEnd) {

	namespace odeint = boost::numeric::odeint;

	// The original GSL implementation used gsl_odeiv2_step_msadams -- an Adams-Moulton
	// multistep method in P(EC)^m *functional iteration* mode, which never actually used
	// the Jacobian. bulirsch_stoer is the closest Boost.Odeint equivalent: adaptive,
	// high-order, and likewise doesn't require a Jacobian.
	odeint::bulirsch_stoer<std::vector<double>> stepper(integrationTolerance_, integrationTolerance_);

	auto system = [this](const std::vector<double>& state, std::vector<double>& dwdy, double y) {
		odeFunction(y, state, dwdy);
	};

	double hStart = (yEnd - yStart)/200;

	try {
		odeint::integrate_adaptive(stepper, system, w, yStart, yEnd, hStart);
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

void TESolver::odeFunction(double y, const std::vector<double>& w, std::vector<double>& f) {

	// w contains the last values of u_n{re, im} and u'_n{re, im}, in that order.
	// need to compute f = dw/dy = u'_n{re, im} followed by u''_n{re, im}

	std::complex<double>* localK2 = k2ForCurrentThread();
	if(computeGratingExpansion(y, localK2) != Result::Success) {
		std::cout << "ODE: Function Error: Cannot compute grating expansion at y = " << y << std::endl;
		throw GratingExpansionError{};	// invalid profile, or y above the profile height.
	}

	for(int i=0; i<fourNp2_; ++i) {
		f[i] = w[i + fourNp2_];
	}

	for(int i=fourNp2_; i<eightNp4_; i+=2) {
		int n = (i - fourNp2_)/2 - N_;
		std::complex<double> upp_n(0,0);

		double alpha2 = alpha_[n + N_];
		alpha2 *= alpha2;

		for(int j=0; j<fourNp2_; j+=2) {
			int m = j/2 - N_;

			std::complex<double> u_m(w[j], w[j+1]);

			std::complex<double> M_nm(0,0);
			if(n-m >= -N_ && n-m <= N_)
				M_nm = -localK2[n-m + N_];
			if(n == m)
				M_nm += alpha2;

			upp_n += M_nm * u_m;
		}

		f[i] = upp_n.real();
		f[i+1] = upp_n.imag();
	}
}

///\todo DISABLED (see TESolver.h). Preserved for reference in case an implicit stepper (e.g. rosenbrock4) is
/// wanted later. dfdw is the Jacobian in row-major order (size eightNp4_ x eightNp4_); dfdy is the explicit
/// y-dependence of odeFunction(), left at 0 here (only exactly true for a rectangular grating -- see original \todo).
/*
void TESolver::odeJacobian(double y, const std::vector<double>& w, std::vector<double>& dfdw, std::vector<double>& dfdy)
{
	(void)w;

	std::complex<double>* localK2 = k2ForCurrentThread();
	if(computeGratingExpansion(y, localK2) != Result::Success) {
		std::cout << "ODE: Jacobian Error: Cannot compute grating expansion at y = " << y << std::endl;
		throw GratingExpansionError{};
	}

	std::fill(dfdw.begin(), dfdw.end(), 0.0);
	std::fill(dfdy.begin(), dfdy.end(), 0.0);

	// upper right-hand block: identity matrix.
	for(int i=0; i<fourNp2_; ++i) {
		dfdw[i*eightNp4_+fourNp2_+i] = 1.0;
	}

	// lower left-hand block: 2x2 submatrices [M_re, -M_im; M_im, M_re] at (i,j)..(i+1,j+1).
	for(int i=fourNp2_; i<eightNp4_; i+=2) {
		int n = (i - fourNp2_)/2 - N_;

		double alpha2 = alpha_[n + N_];
		alpha2 *= alpha2;

		for(int j=0; j<fourNp2_; j+=2) {
			int m = j/2 - N_;

			std::complex<double> M_nm(0,0);
			if(n-m >= -N_ && n-m <= N_)
				M_nm = -localK2[n-m + N_];
			if(n == m)
				M_nm += alpha2;

			dfdw[i*eightNp4_ + j] = M_nm.real();
			dfdw[i*eightNp4_ + j + 1] = -M_nm.imag();
			dfdw[(i+1)*eightNp4_ + j] = M_nm.imag();
			dfdw[(i+1)*eightNp4_ + j + 1] = M_nm.real();
		}
	}

	/// \todo IMPORTANT! Leaving dfdy = 0 for now. This is only true in case of rectangular grating...
}
*/

std::complex<double>* TESolver::k2ForCurrentThread() {
	return k2_[omp_get_thread_num()].data();
}

double TESolver::conditionNumber(const Eigen::MatrixXcd& A) {
	// cond(A) = max(singular values) / min(singular values), 2-norm. Previously incomplete/broken
	// in the GSL version (there was no straightforward way to get a complex matrix norm); Eigen's
	// JacobiSVD gives this directly.
	Eigen::JacobiSVD<Eigen::MatrixXcd> svd(A);
	const auto& sv = svd.singularValues();	// already sorted descending.
	double maxSv = sv(0);
	double minSv = sv(sv.size()-1);
	if(minSv == 0)
		return std::numeric_limits<double>::infinity();
	return maxSv / minSv;
}

void TESolver::setIntegrationStartingValues(std::vector<double>& w, int j, int m)
{
	std::fill(w.begin(), w.end(), 0.0);

	bool secondRound = false;
	if(j >= twoNp1_) {
		secondRound = true;
		j -= twoNp1_;
	}

	w[2*j] = 1.0;

	std::complex<double> uprime = (m == 1 ? beta1_[j] : betaM_[j]) * std::complex<double>(0, secondRound ? 1 : -1);
	w[fourNp2_ + 2*j] = uprime.real();
	w[fourNp2_ + 2*j+1] = uprime.imag();
}

void TESolver::computeAlphaAndBeta(double incidenceDeg)
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

		// betaM_: rayleigh expansion above grating.
		double k22minusAn2 = k_2*k_2 - alpha*alpha;
		if(k22minusAn2 >= 0)
			betaM_[i] = std::complex<double>(sqrt(k22minusAn2), 0);
		else
			betaM_[i] = std::complex<double>(0, sqrt(-k22minusAn2));

		// beta1_: rayleigh expansion inside grating
		std::complex<double> k12minusAn2 = k_1*k_1 - alpha*alpha;
		beta1_[i] = complex_sqrt_upperComplexPlane(k12minusAn2);
	}
}

void TESolver::computeLayers()
{
	double a = g_.totalHeight();

	double magicNumber = 3;	// should be ln(1e15). However, emperically this is not enough to maintain stability (ex: REIXS LEG).  7 = ln(1e3) seems stable for all tests so far.

	numLayers_ = std::max( std::abs(betaM_[0])*a/magicNumber, std::abs(betaM_[2*N_])*a/magicNumber );
	if(numLayers_ < 1)
		numLayers_ = 1;

	M_ = numLayers_+2;

	y_.resize(M_);	// y_[0] unused, so that y_m = y_[m], lowest m=1, highest m=M-1.

	for(int m=1; m<M_; ++m) {
		y_[m] = double(m-1)/numLayers_*a;
	}
}

void TESolver::computeBMFromSMatrix()
{
	double a = g_.totalHeight();

#pragma omp parallel for num_threads(numThreads_)
	for(int i=0; i<twoNp1_; i++) {
		BM_[i] = S12_(i, N_) * std::exp(std::complex<double>(0,-a) * (betaM_[i] + betaM_[N_]));
	}
}

Result::Code TESolver::computeTMatrixBelowLayer(int m, bool printDebugOutput)
{
	bool integrationFailureOccurred = false;

#pragma omp parallel for num_threads(numThreads_) schedule(dynamic) reduction(||:integrationFailureOccurred)
	for(int j=0; j<fourNp2_; ++j) {

		std::vector<double>& w = wVectorForP(j);

		setIntegrationStartingValues(w, j, m-1);

		if(printDebugOutput && omp_get_thread_num() == 0) {
			std::cout << "Initial value u_{p=" << j-N_ << "}(yStart):" <<std::endl;
			std::cout << "     ";
			for(int n=0; n<twoNp1_; ++n)
				std::cout << w[2*n] << "," << w[2*n+1] << "    ";
			std::cout << std::endl;
			std::cout << "Initial value u'_{p=" << j-N_ << "}(yStart):" <<std::endl;
			std::cout << "     ";
			for(int n=0; n<twoNp1_; ++n)
				std::cout << w[fourNp2_ + 2*n] << "," << w[fourNp2_ + 2*n+1] << "    ";
			std::cout << std::endl;
			std::cout << std::endl;
		}

		Result::Code status = integrateTrialSolutionAlongY(w, y_[m-1], y_[m]);

		if(printDebugOutput && omp_get_thread_num() == 0) {
			std::cout << "Final value u_{p=" << j-N_ << "}(yEnd):" <<std::endl;
			std::cout << "     ";
			for(int n=0; n<twoNp1_; ++n)
				std::cout << w[2*n] << "," << w[2*n+1] << "    ";
			std::cout << std::endl;
			std::cout << "Final value u'_{p=" << j-N_ << "}(yEnd):" <<std::endl;
			std::cout << "     ";
			for(int n=0; n<twoNp1_; ++n)
				std::cout << w[fourNp2_ + 2*n] << "," << w[fourNp2_ + 2*n+1] << "    ";
			std::cout << std::endl;
			std::cout << std::endl;
		}

		if(status != Result::Success)
			integrationFailureOccurred = true;
		else {
			// Fill T-matrix at this column. If j >= 2*N+1, we're dealing with the right-side blocks T12_, T22_.
			if(j >= twoNp1_) {
				int jj = j - twoNp1_;

				for(int i=0; i<twoNp1_; ++i) {
					std::complex<double> u_ij(w[2*i], w[2*i+1]);
					std::complex<double> uprime_ij(w[fourNp2_ + 2*i], w[fourNp2_ + 2*i + 1]);
					std::complex<double> temp = uprime_ij / (betaM_[i] * std::complex<double>(0,1)); // = u'_ij / (i*betaM_n)

					T12_(i, jj) = 0.5*(u_ij - temp);
					T22_(i, jj) = 0.5*(u_ij + temp);
				}
			}
			// Otherwise we're filling T11, T21. (left-side blocks).
			else {
				for(int i=0; i<twoNp1_; ++i) {
					std::complex<double> u_ij(w[2*i], w[2*i+1]);
					std::complex<double> uprime_ij(w[fourNp2_ + 2*i], w[fourNp2_ + 2*i + 1]);
					std::complex<double> temp = uprime_ij / (betaM_[i] * std::complex<double>(0,1));

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

// computes the fourier expansion of the multistep function given by values stepsK2 at x-axis locations stepsX, and stores in k2.
void TESolver::computeGratingExpansion(const double *stepsX, const std::complex<double> *stepsK2, int numSteps, std::complex<double> *k2) const
{
	double d = g_.period();
	double K = 2*M_PI/d;

	/// \warning assumes numSteps is in [1, PEG_MAX_PROFILE_CROSSINGS]

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