/*
Copyright (C) 2026 Jesaja Weintritt (jesaja.weintritt@stud.eah-jena.de) and 2012 Mark Boots (mark.boots@usask.ca).

This program was originally implemented as a part of the Parallel Efficiency of Gratings project PEG and got reworked in 2026. PEG is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License, version 3, as published by the Free Software Foundation.
See <http://www.gnu.org/licenses/> for details.

This reworked version contains substantial modifications by Jesaja Weintritt (2026) and has not been independently verified against the original. It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; use at your own risk and verify results independently.
*/

/*
Shared base class for the TE and TM solvers, due their partial similar calculation steps:
- storage
- setup of alpha/beta
- and layer stack
- Fourier expansion of the grating (and the Toeplitz matrix built from it)
- the combined integration of all trial solutions through one layer, and the T-matrix columns
- S-matrix recursion
- conversion of the result into efficiencies

What is polarization-specific stays in TESolver / TMSolver and is plugged in through the virtual functions in the protected section: the equations themselves (columnDerivatives), their per-y preparation (prepareAtY, factorizeAtY), and three small details (secondFieldName, secondFieldIsEpsScaled, transmittedOrderWeight).
*/

#ifndef SolverBase_H
#define SolverBase_H

#include "PEG.h"

#include <complex>
#include <memory>
#include <vector>

#include <Eigen/Dense>

#define PEG_RCOND_FAIL_THRESHOLD 1e-12	// below this, inverse is untrustworthy -> AlgebraFailure
#define PEG_RCOND_WARN_THRESHOLD 1e-8	// below this, just a debug-output warning
#define PEG_ENERGY_TOLERANCE 1e-2		// reflected efficiencies sum to at most =1; + this before the result is rejected


/// Common steps of the solvers
class SolverBase {
public:
	/// Construct the shared solver context for the given \c grating and math options \c mo. \c numThreads specifies how many threads to use for fine parallelization; ideally it should be <= the number of processor cores on your computer / on a single cluster node.
	SolverBase(const Grating& grating, const MathOptions& mo = MathOptions(), int numThreads = 1, bool measureTiming = false);
	virtual ~SolverBase() = default;

	/// Kept non-copyable to avoid accidentally copying large per-thread buffers/matrices.
	SolverBase(const SolverBase&) = delete;
	SolverBase& operator=(const SolverBase&) = delete;

	/// Computes alpha_, beta1_, and betaM_ for all n, based on v_1_ (material refractive index) and wl_. beta_n^2 = k^2 - alpha_n^2 is polarization-independent.
	void computeAlphaAndBeta(double incidenceDeg);

	/// Calculates how many vertical layers (numLayers_ and M_) are sufficient to keep exponentials from contamination. Fills y_ with the vertical coordinate at each layer.
	void computeLayers();

	/// Calculates the grating fourier expansion for k^2_m at a given \c y value and wavelength \c wl, and stores in \c k2 (must have space for 4*N_ + 1 coefficients, since we will be computing from n = -2N_ to 2N.). Reads member variables N_, wavelength wl_, grating refractive index \c v_1_, and grating geometry from \c g_. Returns Result::Success, or Result::InvalidGratingFailure if the profile is not supported or \c y is larger than the groove height.
	Result::Code computeGratingExpansion(double y, std::complex<double>* k2) const;

	/// Computes the Fourier components of the grating expansion k^2_m into \c k2, based on an array of x crossing (step) values \c stepsX and corresponding k^2 values \c stepsK2 immediately to the left of those x values. \c numSteps is the number of steps [usually two or four, if there are interpenetrating coatings)]. Valid only up to PEG_MAX_PROFILE_CROSSINGS (60) to avoid allocating memory, since this function is called repeatedly.
	void computeGratingExpansion(const double* stepsX, const std::complex<double>* stepsK2, int numSteps, std::complex<double>* k2) const;

	/// Computes the \c BM_ outgoing reflected Rayleigh coefficients, based on a finished S matrix (S12_ block).
	void computeBMFromSMatrix();

	/// Returns the square root \c w of a complex number \c z, choosing the branch cut so that Im(w) >= 0.
	/// Note: This is different than the usual branch cut defined for the "principal square root", which uses the negative real axis so that Re(w) >= 0 [ie: w is in the right complex plane].
	static std::complex<double> complex_sqrt_upperComplexPlane(std::complex<double> z);

	/// Returns the condition number of a complex square matrix \c A (ratio of largest to smallest singular value), computed via Eigen's JacobiSVD.
	static double conditionNumber(const Eigen::MatrixXcd& A);

	/// Creates the solver for the given \c polarization (a TESolver or a TMSolver). Callers only need this class: create(), then getEff().
	static std::unique_ptr<SolverBase> create(Polarization polarization, const Grating& grating, const MathOptions& mo = MathOptions(), int numThreads = 1, bool measureTiming = false);

	/// The whole efficiency calculation shared by TE and TM: the efficiencies of all orders at incidence angle \c incidenceDeg and wavelength \c wl. One solver object is meant for one calculation (it overwrites its working storage).
	Result getEff(double incidenceDeg, double wl, double rmsRoughnessNm, bool printDebugOutput);

protected:
	/// Thrown while evaluating the ODE when the grating expansion cannot be computed at some y (invalid geometry, or y above the profile height). Caught in integrateAllTrialSolutionsAlongY() and mapped to Result::InvalidGratingFailure.
	struct GratingExpansionError {};

	// ---- What is polarization-specific (TESolver / TMSolver implement these) ----

	/// Weight used for the transmitted efficiency of substrate order \c i in the debug output (the quantity multiplying |A_i|^2 before dividing by Re(betaM_[N_])). For TE this is Re(beta1_[i]); for TM the eps_r of the substrate has to be divided out.
	virtual double transmittedOrderWeight(int i) const = 0;

	/// Name of the second field of the trial solutions, for the debug output: "u'" (TE) or "w" (TM).
	virtual const char* secondFieldName() const = 0;

	/// true if the second field is w = eps_r^-1 u' (TM), whose start value at the substrate boundary (m == 1) is therefore divided by v_1_^2; false if it is u' itself (TE).
	virtual bool secondFieldIsEpsScaled() const = 0;

	/// Called once per y, before the ODE right-hand side is evaluated, after k2_ and toeplitz_ (the Toeplitz matrix of k^2(y)) have been filled. For extra per-y preparation that is not shared (TM: the Toeplitz matrix of 1/eps_r). May throw GratingExpansionError. Default: nothing.
	virtual void prepareAtY(double y) { (void)y; }

	/// Called once per y by every thread of the team that evaluates the right-hand side, before any column is processed. May contain an OpenMP \c sections block (TM factorizes two matrices here). Default: nothing.
	virtual void factorizeAtY() {}

	/// The polarization-specific equations. \c U and \c W hold a block of trial columns (2N+1 rows each): U = u, W = the second field. Computes d(u)/dy into \c dU and d(second field)/dy into \c dW (same size). Must only read members that prepareAtY() / factorizeAtY() filled.
	virtual void columnDerivatives(const Eigen::MatrixXcd& U, const Eigen::MatrixXcd& W, Eigen::MatrixXcd& dU, Eigen::MatrixXcd& dW) const = 0;

	// ---- What both polarizations share ----

	/// Computes the blocks of the T matrix (T11_, T12_, T21_, T22_) for the layer below \c y_[m]. Since \c m = 1 is the top of the substrate, \c m can range from [2, M-1]. Fills the starting values of all 4N+2 trial solutions, integrates them as ONE combined ODE, and calls storeTMatrixColumn() for each.
	Result::Code computeTMatrixBelowLayer(int m, bool printDebugOutput);

	/// Initializes an 8N+4 array of double [\c u, second field] with the starting values of trial solution \c j for layer \c m: u = delta_{n,p}, second field = -i*beta_n*delta_{n,p} (+i*beta_n for the second round, j >= 2N+1), with beta1_ if m == 1 and betaM_ otherwise. At m == 1 an eps-scaled second field (TM) is divided by v_1_^2.
	void setIntegrationStartingValues(std::vector<double>& w, int j, int m);

	/// Integrates ALL 4N+2 trial solutions from y = \c yStart to y = \c yEnd as one combined ODE (Boost.Odeint, bulirsch_stoer). \c state holds the per-column [u, second field] blocks end-to-end (length fourNp2_*eightNp4_). Returns Success, InvalidGratingFailure (GratingExpansionError), or ConvergenceFailure (any other exception).
	Result::Code integrateAllTrialSolutionsAlongY(std::vector<double>& state, double yStart, double yEnd);

	/// Computes d(state)/dy for all trial solutions at once: the grating expansion and the Toeplitz matrices are built ONCE per y and shared by all columns, then the columns are processed in blocks (one per thread) with the polarization's columnDerivatives().
	void odeFunctionMatrix(double y, const std::vector<double>& state, std::vector<double>& f);

	/// Builds the full (2N+1 x 2N+1) Toeplitz matrix from the Fourier coefficients \c coeffs[p + N], p in [-N, N]. Entries with |n-m| > N are zero.
	static void buildToeplitz(const std::complex<double>* coeffs, int N, Eigen::MatrixXcd& toeplitz);

	/// The x crossings and k^2 values (vacuum, substrate, coating at wl_) of the grating at height \c y; returns the number of steps (< 1: invalid). \c stepsX and \c stepsK2 need room for PEG_MAX_PROFILE_CROSSINGS entries.
	int computeStepsAtY(double y, double* stepsX, std::complex<double>* stepsK2) const;

	/// Stores trial solution \c w (layout [u, second field], see wVectorForP()) of column \c j into the T-matrix blocks: columns j < 2N+1 go to T11_/T21_, the others to T12_/T22_. The second field (u' for TE, w = eps_r^-1 u' for TM) is divided by i*betaM_n.
	void storeTMatrixColumn(const std::vector<double>& w, int j);

	/// Debug output of the initial (\c initial = true) or final values of trial solution \c w of column \c j. \c secondName is the name of the second field ("u'" or "w").
	void printTrialSolution(bool initial, const char* secondName, const std::vector<double>& w, int j) const;

	/// Number of threads to use for this calculation
	int numThreads_;

	/// The number of Fourier coefficients
	int N_;
	/// 2*N_ + 1, 4*N_+2, and 8*N_+4, since these are used a lot
	int twoNp1_, fourNp2_, eightNp4_;

	/// The accuracy to use for numerical integration. Default 1e-8, gets it from mathoption
	double integrationTolerance_;

	/// alpha array (size 2N+1)
	std::vector<double> alpha_;
	/// beta array (size 2N+1). betaM_ is for the superstrate, beta1_ is for the substrate.
	std::vector<std::complex<double>> betaM_, beta1_;
	/// B_n^{M} array: outgoing reflected Rayleigh coefficients. (size 2N+1)
	std::vector<std::complex<double>> BM_;

	/// Number of layers to use in the vertical stack to keep the growing exponentials from numerical contamination.
	int numLayers_;
	/// M-2 = numLayers_.
	int M_;

	/// This block of storage contains the [u, second field] Fourier component vectors. They are arranged with each value {re,im}, from order [-N to N], the u vector followed by the second vector... repeated for each trial solution. The size of one w vector is (8*N_+4), and there are (4*N_+2) trial solutions.
	std::vector<std::vector<double>> wVectors_;

	/// Returns the wVector for a trial solution \c p, where \c p is numbered from [0, 4*N_+1].
	std::vector<double>& wVectorForP(int p) { return wVectors_[p]; }
	/// Fourier coefficients k^2(y) of the grating (size 2N+1) and the Toeplitz matrix built from them, refilled at every y of the ODE.
	std::vector<std::complex<double>> k2_;
	Eigen::MatrixXcd toeplitz_;
	/// Scratch space of the combined integration: the 4N+2 per-column [u, second field] blocks end-to-end.
	std::vector<double> combinedState_;

	/// Blocks of T matrix, used in computation of a single layer.
	Eigen::MatrixXcd T11_, T12_, T21_, T22_;
	/// Blocks of S matrix, used in recursive computation of everything up to current layer.
	Eigen::MatrixXcd S12_, S22_;
	/// Inverse of Z-matrix, used in computation of S. (Note: we don't actually compute any inverses; Zinv_ is directly calculated from Zinv^{q+1} = T11^{q+1} + T12^{q+1} S12^{1}, and then we solve Zinv_ * Z_ = I via Eigen's LU decomposition to avoid loss of precision from an explicit matrix inversion.)
	Eigen::MatrixXcd Zinv_;
	/// This is a 2*N_+1 x 2*N_+1 matrix used as a workspace matrix.
	Eigen::MatrixXcd Z_, workMatrix_;

	/// wavelength for the current calculation
	double wl_;
	/// refractive index of the grating substrate material, at wl_
	std::complex<double> v_1_;
	/// refractive index of the coating material, at wl_ (only meaningful if the grating has a coating)
	std::complex<double> v_c_ = std::complex<double>(0, 0);

	/// The y-coordinate of the infinitely-thin Rayleigh layer at y_m, with m = [1, M_ - 1]. y_[0] is unused, so that we can take y_m = y_[m].
	std::vector<double> y_;

	/// a reference to the grating we're solving
	const Grating& g_;

	/// A flag that indicates that we should measure the time required for all related blocks of operations
	bool measureTiming_;

	/// Stores the timing results:
	double timing_[12];
	/// Used for calculating computation times.
	double time_;
};

#endif