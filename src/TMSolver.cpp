/*
Copyright (C) 2026 Jesaja Weintritt (jesaja.weintritt@stud.eah-jena.de) and 2012 Mark Boots (mark.boots@usask.ca).

This program was originally implemented as a part of the Parallel Efficiency of Gratings project PEG and got reworked in 2026. PEG is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License, version 3, as published by the Free Software Foundation.
See <http://www.gnu.org/licenses/> for details.

This reworked version contains substantial modifications by Jesaja Weintritt (2026) and has not been independently verified against the original. It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; use at your own risk and verify results independently.
*/


#include "TMSolver.h"

#include <cmath>
#include <iostream>
#include <omp.h>

TMSolver::TMSolver(const Grating& grating, const MathOptions& mo, int numThreads, bool measureTiming)
	: SolverBase(grating, mo, numThreads, measureTiming)
{
	k2Inv_.resize(twoNp1_);
	toeplitzInv_.resize(twoNp1_, twoNp1_);

	timing_[0] = omp_get_wtime() - time_;
}

Result::Code TMSolver::computeInvEpsExpansion(double y, std::complex<double>* invEps) const {

	// Same step geometry as computeGratingExpansion(double,...), but Fourier-expand 1/eps_r(y) instead of k^2(y) -- see 2026-09-29 class-level note in TMSolver.h for why the u' equation needs this second expansion, computed via Li's inverse rule.

	double stepsX[PEG_MAX_PROFILE_CROSSINGS];
	std::complex<double> stepsK2[PEG_MAX_PROFILE_CROSSINGS];

	int numSteps = computeStepsAtY(y, stepsX, stepsK2);
	if(numSteps < 1)
		return Result::InvalidGratingFailure;

	// eps_r = k2/k2_M (k2_M = k_0^2 is material-independent), so 1/eps_r = k2_M/k2 -- dimensionless, as needed: computeGratingExpansion()'s Fourier-coefficient math doesn't care what physical quantity the step values represent.
	std::complex<double> k_M(2 * M_PI / wl_, 0);
	std::complex<double> k2_M = k_M * k_M;

	std::complex<double> stepsInvEps[PEG_MAX_PROFILE_CROSSINGS];
	for(int p=0; p<numSteps; ++p)
		stepsInvEps[p] = k2_M / stepsK2[p];

	computeGratingExpansion(stepsX, stepsInvEps, numSteps, invEps);

	return Result::Success;
}

void TMSolver::prepareAtY(double y)
{
	if(computeInvEpsExpansion(y, k2Inv_.data()) != Result::Success) {
		std::cout << "ODE: Function Error: Cannot compute 1/eps grating expansion at y = " << y << std::endl;
		throw GratingExpansionError{};
	}
	buildToeplitz(k2Inv_.data(), N_, toeplitzInv_);
}

void TMSolver::factorizeAtY()
{
	// the two independent LU factorizations run as two sections; the implicit barrier at the end of the sections block makes sure both are complete before any thread starts solving.
#pragma omp sections
	{
#pragma omp section
		luInvEps_.compute(toeplitzInv_);
#pragma omp section
		lu_.compute(toeplitz_);
	}
}

void TMSolver::columnDerivatives(const Eigen::MatrixXcd& U, const Eigen::MatrixXcd& W, Eigen::MatrixXcd& dU, Eigen::MatrixXcd& dW) const
{
	const double k0 = 2 * M_PI / wl_;
	const double k0sq = k0 * k0;
	const int nc = static_cast<int>(U.cols());

	Eigen::MatrixXcd alphaU(twoNp1_, nc);
	for(int i=0; i<twoNp1_; ++i)
		alphaU.row(i) = alpha_[i] * U.row(i);

	// the two equations, reusing the ONE LU factorization from factorizeAtY() for all columns
	dU = luInvEps_.solve(W);
	Eigen::MatrixXcd V = lu_.solve(alphaU);

	dW.resize(twoNp1_, nc);
	for(int jj=0; jj<nc; ++jj)
		for(int i=0; i<twoNp1_; ++i)
			dW(i,jj) = k0sq * (alpha_[i] * V(i,jj) - U(i,jj));
}

double TMSolver::transmittedOrderWeight(int i) const
{
	// Transmitted order i (substrate, eps_r = v_1_^2) and the reference order N_ (vacuum, eps_r = 1) do NOT share the same medium, so the eps_r ratio (TMSolver_alg_ref.md \S5.2/\S7.7) does not cancel here.
	std::complex<double> epsSub = v_1_ * v_1_;
	std::complex<double> ratio_i = beta1_[i] / epsSub;
	return ratio_i.real();
}