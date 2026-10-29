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
	struct GratingExpansionError {};
}

TMSolver::TMSolver(const Grating& grating, const MathOptions& mo, int numThreads, bool measureTiming)
	: SolverBase(grating, mo, numThreads, measureTiming)
{
	// *** 2026-09-29 performance fix: single buffers, not one-per-thread -- odeFunctionMatrix() is
	// called by one shared stepper now (see class-level / integrateAllTrialSolutionsAlongY() notes).
	k2_.resize(twoNp1_);
	toeplitz_.resize(twoNp1_, twoNp1_);
	k2Inv_.resize(twoNp1_);
	toeplitzInv_.resize(twoNp1_, twoNp1_);

	// concatenated [u,w] blocks for all fourNp2_ trial columns -- see combinedState_'s header comment.
	combinedState_.resize(fourNp2_ * eightNp4_);

	timing_[0] = omp_get_wtime() - time_;
}

Result TMSolver::getEffTM(double incidenceDeg, double wl, double rmsRoughnessNm, bool printDebugOutput) {
	return getEff(incidenceDeg, wl, rmsRoughnessNm, printDebugOutput);
}

Result::Code TMSolver::computeInvEpsExpansion(double y, std::complex<double>* invEps) const {

	// Same step geometry as computeGratingExpansion(double,...), but Fourier-expand 1/eps_r(y) instead of
	// k^2(y) -- see 2026-09-29 class-level note in TMSolver.h for why the u' equation needs this second
	// expansion, computed via Li's inverse rule.

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

	// eps_r = k2/k2_M (k2_M = k_0^2 is material-independent), so 1/eps_r = k2_M/k2 -- dimensionless, as
	// needed: computeGratingExpansion()'s Fourier-coefficient math doesn't care what physical quantity the
	// step values represent.
	std::complex<double> stepsInvEps[PEG_MAX_PROFILE_CROSSINGS];
	for(int p=0; p<numSteps; ++p)
		stepsInvEps[p] = k2_M / stepsK2[p];

	computeGratingExpansion(stepsX, stepsInvEps, numSteps, invEps);

	return Result::Success;
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

void TMSolver::odeFunctionMatrix(double y, const std::vector<double>& state, std::vector<double>& f)
{
	// Build both Toeplitz matrices and factorize them ONCE for this y -- shared across all fourNp2_
	// trial columns below. This one-per-y (not one-per-column-per-y) cost is the entire point of
	// batching; see the performance note on this function in TMSolver.h.
	if(computeGratingExpansion(y, k2_.data()) != Result::Success) {
		std::cout << "ODE: Function Error: Cannot compute grating expansion at y = " << y << std::endl;
		throw GratingExpansionError{};
	}
	buildToeplitz(k2_.data(), N_, toeplitz_);

	if(computeInvEpsExpansion(y, k2Inv_.data()) != Result::Success) {
		std::cout << "ODE: Function Error: Cannot compute 1/eps grating expansion at y = " << y << std::endl;
		throw GratingExpansionError{};
	}
	buildToeplitz(k2Inv_.data(), N_, toeplitzInv_);

	Eigen::PartialPivLU<Eigen::MatrixXcd> luInvEps(toeplitzInv_);
	Eigen::PartialPivLU<Eigen::MatrixXcd> lu(toeplitz_);

	// Unpack all fourNp2_ columns' [u,w] out of the flat, concatenated state (see combinedState_'s
	// header comment for the layout -- identical to a single column's w_arr, repeated fourNp2_ times)
	// into (twoNp1_ x fourNp2_) matrices, one column of each matrix per trial solution.
	Eigen::MatrixXcd U(twoNp1_, fourNp2_), W(twoNp1_, fourNp2_), alphaU(twoNp1_, fourNp2_);
	for(int j=0; j<fourNp2_; ++j) {
		const double* col = state.data() + j*eightNp4_;
		for(int i=0; i<twoNp1_; ++i) {
			U(i,j) = std::complex<double>(col[2*i], col[2*i+1]);
			W(i,j) = std::complex<double>(col[fourNp2_ + 2*i], col[fourNp2_ + 2*i + 1]);
		}
	}
	for(int i=0; i<twoNp1_; ++i)
		alphaU.row(i) = alpha_[i] * U.row(i);

	double k0 = 2 * M_PI / wl_;
	double k0sq = k0 * k0;

	// Same two equations as the old per-column odeFunction() (see class-level note for the physics),
	// now with a MatrixXcd right-hand side: each solve() below reuses ONE LU factorization across
	// every column at once, rather than refactorizing per column as the old code did.
	Eigen::MatrixXcd dUdy = luInvEps.solve(W);
	Eigen::MatrixXcd V = lu.solve(alphaU);

	for(int j=0; j<fourNp2_; ++j) {
		double* col = f.data() + j*eightNp4_;
		for(int i=0; i<twoNp1_; ++i) {
			std::complex<double> dwDy_ij = k0sq * (alpha_[i] * V(i,j) - U(i,j));
			col[2*i]              = dUdy(i,j).real();
			col[2*i+1]            = dUdy(i,j).imag();
			col[fourNp2_ + 2*i]   = dwDy_ij.real();
			col[fourNp2_ + 2*i+1] = dwDy_ij.imag();
		}
	}
}

Result::Code TMSolver::integrateAllTrialSolutionsAlongY(std::vector<double>& state, double yStart, double yEnd) {

	namespace odeint = boost::numeric::odeint;

	odeint::bulirsch_stoer<std::vector<double>> stepper(integrationTolerance_, integrationTolerance_);

	auto system = [this](const std::vector<double>& s, std::vector<double>& dsdy, double y) {
		odeFunctionMatrix(y, s, dsdy);
	};

	double hStart = (yEnd - yStart)/200;

	try {
		odeint::integrate_adaptive(stepper, system, state, yStart, yEnd, hStart);
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

double TMSolver::transmittedOrderWeight(int i) const
{
	// Transmitted order i (substrate, eps_r = v_1_^2) and the reference order N_ (vacuum, eps_r = 1) do NOT share the
	// same medium, so the eps_r ratio (TMSolver_alg_ref.md \S5.2/\S7.7) does not cancel here.
	std::complex<double> epsSub = v_1_ * v_1_;
	std::complex<double> ratio_i = beta1_[i] / epsSub;
	return ratio_i.real();
}

Result::Code TMSolver::computeTMatrixBelowLayer(int m, bool printDebugOutput)
{
	// *** 2026-09-29 performance fix: fill every column's starting values (setIntegrationStartingValues()
	// itself is unchanged), then integrate all fourNp2_ trial solutions as ONE combined ODE instead of
	// fourNp2_ independent ones under a #pragma omp parallel for -- see the header comment on
	// integrateAllTrialSolutionsAlongY() / odeFunctionMatrix() for why this is faster even on one core.
	for(int j=0; j<fourNp2_; ++j) {
		std::vector<double>& w_arr = wVectorForP(j);
		setIntegrationStartingValues(w_arr, j, m-1);
		std::copy(w_arr.begin(), w_arr.end(), combinedState_.begin() + j*eightNp4_);

		if(printDebugOutput && j == 0)
			printTrialSolution(true, "w", w_arr, j);
	}

	Result::Code status = integrateAllTrialSolutionsAlongY(combinedState_, y_[m-1], y_[m]);
	if(status != Result::Success)
		return status;

	for(int j=0; j<fourNp2_; ++j) {
		std::vector<double>& w_arr = wVectorForP(j);
		std::copy(combinedState_.begin() + j*eightNp4_, combinedState_.begin() + (j+1)*eightNp4_, w_arr.begin());

		if(printDebugOutput && j == 0)
			printTrialSolution(false, "w", w_arr, j);

		// Transform formula is unchanged from TESolver's: this construction is always referenced
		// against betaM_ (vacuum, eps_r=1), so the eps_r factor from TMSolver_alg_ref.md \S3.2
		// evaluates to 1 here -- see class-level comment in TMSolver.h.
		storeTMatrixColumn(w_arr, j);
	}

	return Result::Success;
}