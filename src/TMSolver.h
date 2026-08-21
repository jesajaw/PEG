/*
Copyright (C) 2026 Jesaja Weintritt (jesaja.weintritt@stud.eah-jena.de) and 2012 Mark Boots (mark.boots@usask.ca).

This program was originally implemented as a part of the Parallel Efficiency of Gratings project PEG and got reworked in 2026. PEG is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License, version 3, as published by the Free Software Foundation.
See <http://www.gnu.org/licenses/> for details.

This reworked version contains substantial modifications by Jesaja Weintritt (2026) and has not been independently verified against the original. It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; use at your own risk and verify results independently.
*/


#ifndef TMSolver_H
#define TMSolver_H

#include "PEG.h"

#include <complex>
#include <vector>

#include <Eigen/Dense>

/// Contains the context (memory structures, etc.) and algorithm for solving the grating efficiency
/// for TM polarization. Structurally mirrors TESolver; see TMSolver_alg_ref.md for the physics/math
/// that differs, and TESolver_alg_ref.md for everything that is shared (Floquet expansion, transfer-/
/// S-matrix formalism, layering, mode truncation, energy conservation).
///
/// Key structural difference from TESolver: the state carried through the ODE integration and the
/// T-matrix construction is [u, w] rather than [u, u'], where w = eps_r^{-1} * u' is the quantity
/// physically proportional to the tangential E-field (E_x) and hence the one that must be continuous
/// across layer boundaries (see TMSolver_alg_ref.md \S3.1).
///
/// H(y) = Toeplitz(eps_r(y))^{-1} is used throughout (Li's "inverse rule" factorization, chosen for
/// correct convergence at material discontinuities -- see TMSolver_alg_ref.md \S1.3/\S7.3). In terms of
/// k^2 = k_0^2 * eps_r (the same quantity TESolver already computes via computeGratingExpansion()):
///   u'(y) = (1/k_0^2) * Toeplitz(k^2(y)) * w(y)                         -- plain multiply, no solve
///   w'(y) = k_0^2 * ( D_alpha * v(y) - u(y) ),  Toeplitz(k^2(y)) * v(y) = D_alpha * u(y)   -- one LU solve
/// per odeFunction() evaluation.

class TMSolver {
public:
	/// Construct a solver context for the given \c grating and math options \c mo.  \c numThreads specifies how many threads to use for fine parallelization; ideally it should be <= the number of processor cores on your computer / on a single cluster node.
	TMSolver(const Grating& grating, const MathOptions& mo = MathOptions(), int numThreads = 1, bool measureTiming = false);
	/// Destroy a solver context. All storage is owned by std::vector / Eigen types with automatic memory management, so there is nothing left to free manually.
	~TMSolver() = default;

	/// Kept non-copyable to avoid accidentally copying large per-thread buffers/matrices (same rationale as TESolver).
	TMSolver(const TMSolver&) = delete;
	TMSolver& operator=(const TMSolver&) = delete;

	/// Calculates the efficiency at incidence angle \c incidenceDeg and wavelength \c wl.  Side effects: sets the refractive index member variable v_1_; modifies the contents of u_, w_, alpha_, beta_, etc.
	Result getEffTM(double incidenceDeg, double wl, double rmsRoughnessNm = 0, bool printDebugOutput = false);


	// Solving implementation functions
	////////////////////////////////////////

	/// Computes alpha_, beta1_, and betaM_ for all n. Identical in derivation and result to TESolver::computeAlphaAndBeta() -- see TMSolver_alg_ref.md \S2.2/\S7.4 (beta_n^2 = k^2 - alpha_n^2 is polarization-independent).
	void computeAlphaAndBeta(double incidenceDeg);

	/// Calculates how many vertical layers (numLayers_ and M_) are sufficient to keep exponentials from contamination.  Fills y_ with the vertical coordinate at each layer. Identical to TESolver::computeLayers().
	void computeLayers();

	/// Computes the blocks of the T matrix (T11_, T12_, T21_, T22_) for the layer below \c y_[m]. Amplitude-transform formula is the same as TESolver's (the eps_r factor from TMSolver_alg_ref.md \S3.2 evaluates to 1 here, since this construction is always referenced to the vacuum superstrate via betaM_ -- see class-level comment). \c m can range from [2, M-1].
	Result::Code computeTMatrixBelowLayer(int m, bool printDebugOutput = false);

	/// Calculates the grating fourier expansion for k^2_m at a given \c y value and wavelength \c wl, and stores in \c k2. Identical to TESolver::computeGratingExpansion(double, std::complex<double>*) -- reused as-is, since H(y) is built from the same k^2 = k_0^2*eps_r coefficients (see class-level comment). \c k2 must have space for 4*N_ + 2 coefficients.
	Result::Code computeGratingExpansion(double y, std::complex<double>* k2) const;

	/// Computes the Fourier components of the grating expansion k^2_m into \c k2. Identical to TESolver's overload of the same name.
	void computeGratingExpansion(const double* stepsX, const std::complex<double>* stepsK2, int numSteps, std::complex<double>* k2) const;

	/// Builds the full (2N_+1 x 2N_+1) Toeplitz matrix from a Fourier coefficient buffer of size 4*N_+2 (as produced by computeGratingExpansion()), for use in odeFunction()'s linear solve against Toeplitz(k^2(y)).
	static void buildToeplitz(const std::complex<double>* coeffs, int N, Eigen::MatrixXcd& toeplitz);

	/// Initializes an 8N+4 array of double [\c u, \c w] to contain the starting integration values.  The first half of the array \c w_arr contains \c u, the second half contains \c w = eps_r^{-1}*u'. The u value is set to $\delta_{n,p}$; the w value is set to $(\mp i \beta_n^{(M)} \delta_{n,p})$ if layer \c m > 1 (vacuum-referenced, eps_r=1), or to $(\mp i \beta_n^{(1)} \delta_{n,p}) / v_1^2$ if \c m == 1 (substrate-referenced, eps_r = v_1^2) -- see TMSolver_alg_ref.md \S7.6.
	void setIntegrationStartingValues(std::vector<double>& w_arr, int p, int m);

	/// Integrates the [u, w] state contained in \c w_arr from y = \c yStart to y = \c yEnd, using Boost.Odeint. Structurally identical to TESolver::integrateTrialSolutionAlongY() (same stepper, same error handling); only odeFunction() differs.
	Result::Code integrateTrialSolutionAlongY(std::vector<double>& w_arr, double yStart, double yEnd);

	/// Computes the derivative d[u,w]/dy at a given \c y for the ODE stepper. Unlike TESolver::odeFunction() (a direct O(N^2) sum), this requires building Toeplitz(k^2(y)) and performing one (2N_+1)-dimensional complex linear solve per call -- see class-level comment for the exact equations.
	void odeFunction(double y, const std::vector<double>& w_arr, std::vector<double>& f);

	/// Computes the \c BM_ outgoing reflected Rayleigh coefficients, based on a finished S matrix (S12_ block). Identical to TESolver::computeBMFromSMatrix() -- the phase-referenced S-matrix-column extraction (TMSolver_alg_ref.md \S0) does not depend on polarization.
	void computeBMFromSMatrix();


	// General Mathematical Helper functions:
	///////////////////////////////

	/// Returns the square root \c w of a complex number \c z, choosing the branch cut so that Im(w) >= 0. Identical to TESolver's version.
	static std::complex<double> complex_sqrt_upperComplexPlane(std::complex<double> z);

	/// Returns the condition number of a complex square matrix \c A. Identical to TESolver's version.
	static double conditionNumber(const Eigen::MatrixXcd& A);


protected:

	int numThreads_;
	int N_;
	int twoNp1_, fourNp2_, eightNp4_;
	double integrationTolerance_;

	std::vector<double> alpha_;
	std::vector<std::complex<double>> betaM_, beta1_;
	std::vector<std::complex<double>> BM_;

	int numLayers_;
	int M_;

	/// pre-allocated storage for the grating k^2 fourier coefficients (size 4N_+2, one buffer per thread). Same role/layout as TESolver::k2_.
	std::vector<std::vector<std::complex<double>>> k2_;
	std::complex<double>* k2ForCurrentThread();

	/// pre-allocated (2N_+1 x 2N_+1) Toeplitz(k^2(y)) workspace, one per thread, rebuilt at each odeFunction() call from k2ForCurrentThread()'s buffer via buildToeplitz(). Avoids re-allocating the matrix on every call.
	std::vector<Eigen::MatrixXcd> toeplitz_;
	Eigen::MatrixXcd& toeplitzForCurrentThread();

	/// This block of storage contains the [u, w] Fourier component vectors, laid out identically to TESolver::wVectors_ (u followed by w, each entry {re,im}, repeated per trial solution).
	std::vector<std::vector<double>> wVectors_;

	std::vector<double>& wVectorForP(int p) { return wVectors_[p]; }
	std::complex<double>* u(int i, int j) {
		return reinterpret_cast<std::complex<double>*>(wVectorForP(j).data() + 2*i);
	}
	/// Returns w_n = eps_r^{-1} u'_n for order \c n (index \c i) and trial solution \c p (index \c j). Named wComp() (not uprime()) since it is not the plain field derivative -- see class-level comment.
	std::complex<double>* wComp(int i, int j) {
		return reinterpret_cast<std::complex<double>*>(wVectorForP(j).data() + fourNp2_ + 2*i);
	}

	Eigen::MatrixXcd T11_, T12_, T21_, T22_;
	Eigen::MatrixXcd S12_, S22_;
	Eigen::MatrixXcd Zinv_;
	Eigen::MatrixXcd Z_, workMatrix_;

	double wl_;
	std::complex<double> v_1_;
	std::complex<double> v_c_ = std::complex<double>(0, 0);

	std::vector<double> y_;

	const Grating& g_;

	bool measureTiming_;
	double timing_[12];
	double time_;
};



#endif