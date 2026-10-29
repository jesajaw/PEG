/*
Copyright (C) 2026 Jesaja Weintritt (jesaja.weintritt@stud.eah-jena.de) and 2012 Mark Boots (mark.boots@usask.ca).

This program was originally implemented as a part of the Parallel Efficiency of Gratings project PEG and got reworked in 2026. PEG is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License, version 3, as published by the Free Software Foundation.
See <http://www.gnu.org/licenses/> for details.

This reworked version contains substantial modifications by Jesaja Weintritt (2026) and has not been independently verified against the original. It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; use at your own risk and verify results independently.
*/


#ifndef TESolver_H
#define TESolver_H

#include "SolverBase.h"

#include <complex>
#include <vector>

#include <Eigen/Dense>

/// Contains the context (memory structures, etc.) and algorithm for solving the grating efficiency.

class TESolver : public SolverBase{
public:
	/// Construct a solver context for the given \c grating and math options \c mo.  \c numThreads specifies how many threads to use for fine parallelization
	TESolver(const Grating& grating, const MathOptions& mo = MathOptions(), int numThreads = 1, bool measureTiming = false);
	/// Destroy a solver context. All storage is owned by std::vector / Eigen types with automatic memory management, so there is nothing left to free manually.
	~TESolver() = default;

	/// Kept non-copyable to avoid accidentally copying large per-thread buffers/matrices; nothing about the storage itself requires this anymore (unlike before, where copying would have caused a double-free of the GSL-owned memory). Revisit if a genuine use case for copying/moving arises.
	TESolver(const TESolver&) = delete;
	TESolver& operator=(const TESolver&) = delete;
	
	/// Calculates the efficiency at incidence angle \c incidenceDeg and wavelength \c wl.  Side effects: sets the refractive index member variable v_1_; modifies the contents of u_, uprime_, alpha_, beta_, etc.
	Result getEffTE(double incidenceDeg, double wl, double rmsRoughnessNm = 0, bool printDebugOutput = false);


	/// Computes the blocks of the T matrix (T11_, T12_, T21_, T22_) for the layer below \c y_[m].  Since \c m = 1 is the top of the substrate, \c m can range from [2, M-1].
	Result::Code computeTMatrixBelowLayer(int m, bool printDebugOutput = false);

	/// Initializes an 8N+4 array of double [\c u, \c uprime] to contain the starting integration values of the electric field Fourier components. The first half of the array \c w contains \c u, the second half contains \c uprime, with each entry in {re,im} order.  The u value is set to $\delta_{n,p}$ and the u' value is set to $-i \beta_n^{(M)} \delta_{n,p}$ or $i \beta_n^{(M)} \delta_{n,p}$, depending on whether p > 2N_.  (Note n,p here are using numbering from 0, and that for the delta functions, p is aliased back onto [0, 2N] once it reaches 2N_+1.
	/// ! If layer \c m = 1, then uses beta1_n instead of betaM_n for the derivative.
	void setIntegrationStartingValues(std::vector<double>& w, int p, int m);

	/// Integrates the electric field Fourier component vectors contained in \c w from y = \c yStart to y = \c yEnd, using Boost.Odeint. \c w should contain vector \c u followed by \c uprime, with each entry in {re,im} order. Calls computeGratingExpansion() at each y value, so reads member variables N_, v_1_, and g_. Results are returned in-place.
	/// \todo Exact state/stepper type still open -- see accompanying question about which Boost.Odeint stepper to use.
	Result::Code integrateTrialSolutionAlongY(std::vector<double>& w, double yStart, double yEnd);

	/// Computes the derivative dw/dy at a given \c y for the ODE stepper, in the same [u, uprime] layout as integrateTrialSolutionAlongY().
	/// \todo Final signature depends on the chosen Boost.Odeint stepper/state type (see accompanying question).
	void odeFunction(double y, const std::vector<double>& w, std::vector<double>& f);

	///\todo DISABLED. Only needed if an implicit stepper (e.g. rosenbrock4) is used instead of bulirsch_stoer;
	/// the current bulirsch_stoer-based integrateTrialSolutionAlongY() never calls this. Kept as a reference
	/// declaration -- the math is preserved as a commented-out implementation in TESolver.cpp.
	///void odeJacobian(double y, const std::vector<double>& w, std::vector<double>& dfdw, std::vector<double>& dfdy);


	/// Computes the \c BM_ outgoing reflected Rayleigh coefficients, based on a finished S matrix (S12_ block).
	/// void computeBMFromSMatrix();


protected:

	/// Weight of substrate order \c i in the debug transmission sum (see SolverBase::transmittedOrderWeight()).
	double transmittedOrderWeight(int i) const override;

	/// pre-allocated storage for the grating k^2 fourier coefficients.  There is one inner array (per outer vector entry) for each thread to use.  Each inner array has size 2N+1.
	std::vector<std::vector<std::complex<double>>> k2_;
	/// Helper function: returns the k^2 array that should be used by a given thread.
	std::complex<double>* k2ForCurrentThread();

	/// Returns the electric field Fourier component derivative \c u' for order \c n (index \c i) and trial solution \c p (index \c j), out of wVectors_.  i and j are numbered from 0.
	std::complex<double>* uprime(int i, int j) {
		return reinterpret_cast<std::complex<double>*>(wVectorForP(j).data() + fourNp2_ + 2*i);
	}

};

#endif