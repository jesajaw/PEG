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
/// *** 2026-09-29 correctness fix: the two odeFunction() equations below each couple a discontinuous
/// eps_r(x,y) to a field with a complementary discontinuity, and Li's factorization rules (Li, JOSA A
/// 13, 1996, "Use of Fourier series in the analysis of discontinuous periodic structures") require
/// each equation to use whichever Toeplitz
/// (eps_r, or its reciprocal) is inverted, NOT the same one for both. Only the w' equation was using the
/// correct rule; the u' equation was previously using Toeplitz(k^2(y)) directly (equivalent to the
/// *direct* rule for eps_r), which converges very poorly for absorbing/high-contrast (e.g. metallic)
/// materials. Both are now correct:
///   u'(y) = H(y) * w(y),           H(y)   = Toeplitz(1/eps_r(y))^{-1}     -- *inverse* rule.
///           eps_r * w is proportional to D_x = eps_r * E_x, the continuous electric displacement, so
///           its coefficients must come from inverting Toeplitz(1/eps_r), not from Toeplitz(eps_r).
///   w'(y) = k_0^2 * ( D_alpha * v(y) - u(y) ),  Toeplitz(k^2(y)) * v(y) = D_alpha * u(y)   -- *direct*
///           rule (unchanged): this term is proportional to the continuous tangential field E_y, and
///           here it is eps_r itself (not its reciprocal) whose Toeplitz matrix gets inverted.
/// (TE has no analogous issue: TESolver::odeFunction() never inverts a Toeplitz matrix at all.)

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

	/// Same role as computeGratingExpansion(double, complex<double>*), but expands 1/eps_r(y) instead of
	/// k^2(y) = k_0^2*eps_r(y) -- needed for Li's inverse-rule H(y) = Toeplitz(1/eps_r(y))^{-1} in odeFunction()'s
	/// u' equation (see 2026-09-29 class-level note). \c invEps must have space for 4*N_ + 2 coefficients,
	/// matching computeGratingExpansion()'s convention.
	Result::Code computeInvEpsExpansion(double y, std::complex<double>* invEps) const;

	/// Builds the full (2N_+1 x 2N_+1) Toeplitz matrix from a Fourier coefficient buffer of size 4*N_+2 (as produced by computeGratingExpansion()), for use in odeFunction()'s linear solve against Toeplitz(k^2(y)).
	static void buildToeplitz(const std::complex<double>* coeffs, int N, Eigen::MatrixXcd& toeplitz);

	/// Initializes an 8N+4 array of double [\c u, \c w] to contain the starting integration values.  The first half of the array \c w_arr contains \c u, the second half contains \c w = eps_r^{-1}*u'. The u value is set to $\delta_{n,p}$; the w value is set to $(\mp i \beta_n^{(M)} \delta_{n,p})$ if layer \c m > 1 (vacuum-referenced, eps_r=1), or to $(\mp i \beta_n^{(1)} \delta_{n,p}) / v_1^2$ if \c m == 1 (substrate-referenced, eps_r = v_1^2) -- see TMSolver_alg_ref.md \S7.6.
	void setIntegrationStartingValues(std::vector<double>& w_arr, int p, int m);

	/// *** 2026-09-29 performance fix: integrates ALL fourNp2_ trial solutions for the layer between
	/// y = \c yStart and y = \c yEnd as ONE combined ODE, rather than fourNp2_ independent ones (see
	/// odeFunctionMatrix() for why). Replaces the old per-column integrateTrialSolutionAlongY() /
	/// odeFunction() pair, which is why computeTMatrixBelowLayer() no longer needs a
	/// \c \#pragma \c omp \c parallel \c for over columns: the batching already does the work sharing
	/// that parallelization was compensating for, and does more of it than \c numThreads_ threads
	/// plausibly could (see note below). \c state is the fourNp2_ per-column [u,w] blocks (see
	/// wVectorForP() layout) concatenated end-to-end, length fourNp2_*eightNp4_.
	Result::Code integrateAllTrialSolutionsAlongY(std::vector<double>& state, double yStart, double yEnd);

	/// Computes d(state)/dy for ALL fourNp2_ trial solutions at once. The key saving over the old
	/// per-column odeFunction(): building the two (2N_+1 x 2N_+1) Toeplitz matrices and LU-factorizing
	/// them is O(N^3) and was previously repeated once per column (fourNp2_ = 4N_+2 times) at every y
	/// the stepper visited; here it happens ONCE per y, and each factorization's solve() is then
	/// reused across all fourNp2_ columns as one matrix right-hand side (each such solve is only
	/// O(N^2) per extra column). For the N this solver normally runs at, the O(N^3) factorization
	/// dominates the O(N^2) per-column solves, so batching cuts the total work roughly fourNp2_-fold
	/// -- typically far more than \c numThreads_ way parallelism recovers, which is why this replaces
	/// (rather than runs alongside) the old \c \#pragma \c omp \c parallel \c for over columns.
	void odeFunctionMatrix(double y, const std::vector<double>& state, std::vector<double>& f);

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

	/// pre-allocated storage for the grating k^2 fourier coefficients (actual size 2N_+1, matching the
	/// banded |n-m|<=N_ coupling in buildToeplitz() -- see accompanying note there). *** 2026-09-29:
	/// no longer one-buffer-per-thread -- odeFunctionMatrix() is called by a single shared stepper
	/// (integrateAllTrialSolutionsAlongY() runs one combined integration, not fourNp2_ parallel ones),
	/// so only one buffer is ever in use at a time now.
	std::vector<std::complex<double>> k2_;
	/// pre-allocated (2N_+1 x 2N_+1) Toeplitz(k^2(y)) workspace, rebuilt at each odeFunctionMatrix()
	/// call from k2_ via buildToeplitz(). Avoids re-allocating the matrix on every call.
	Eigen::MatrixXcd toeplitz_;

	/// pre-allocated storage for the grating 1/eps_r fourier coefficients (size 2N_+1, parallel to
	/// k2_) -- see computeInvEpsExpansion() and the 2026-09-29 correctness-fix note above.
	std::vector<std::complex<double>> k2Inv_;
	/// pre-allocated (2N_+1 x 2N_+1) Toeplitz(1/eps_r(y)) workspace. Parallel to toeplitz_.
	Eigen::MatrixXcd toeplitzInv_;

	/// *** 2026-09-29 performance fix: scratch space for the batched integration -- the fourNp2_
	/// per-column [u,w] blocks written by setIntegrationStartingValues() into wVectors_, concatenated
	/// end-to-end (length fourNp2_*eightNp4_) so integrateAllTrialSolutionsAlongY() can integrate all
	/// of them as one ODE. See odeFunctionMatrix().
	std::vector<double> combinedState_;

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