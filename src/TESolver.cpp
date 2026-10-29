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
	/// Thrown from TESolver::odeFunction() when the grating expansion cannot be computed at some y (invalid geometry, or y above the profile height). Caught in integrateTrialSolutionAlongY() and mapped to Result::InvalidGratingFailure.
	struct GratingExpansionError {};
}

TESolver::TESolver(const Grating& grating, const MathOptions& mo, int numThreads, bool measureTiming)
	: SolverBase(grating, mo, numThreads, measureTiming)
{
	// one k2_ array per thread, since they will be used simultaneously.
	k2_.assign(numThreads_, std::vector<std::complex<double>>(twoNp1_));

	timing_[0] = omp_get_wtime() - time_;		// time to allocate memory.
}

Result TESolver::getEffTE(double incidenceDeg, double wl, double rmsRoughnessNm, bool printDebugOutput) {
	return getEff(incidenceDeg, wl, rmsRoughnessNm, printDebugOutput);
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

double TESolver::transmittedOrderWeight(int i) const
{
	return beta1_[i].real();
}

Result::Code TESolver::computeTMatrixBelowLayer(int m, bool printDebugOutput)
{
	bool integrationFailureOccurred = false;

#pragma omp parallel for num_threads(numThreads_) schedule(dynamic) reduction(||:integrationFailureOccurred)
	for(int j=0; j<fourNp2_; ++j) {

		std::vector<double>& w = wVectorForP(j);

		setIntegrationStartingValues(w, j, m-1);

		if(printDebugOutput && omp_get_thread_num() == 0)
			printTrialSolution(true, "u'", w, j);

		Result::Code status = integrateTrialSolutionAlongY(w, y_[m-1], y_[m]);

		if(printDebugOutput && omp_get_thread_num() == 0)
			printTrialSolution(false, "u'", w, j);

		if(status != Result::Success)
			integrationFailureOccurred = true;
		else
			storeTMatrixColumn(w, j);
	}

	if(integrationFailureOccurred)
		return Result::ConvergenceFailure;
	else
		return Result::Success;
}
