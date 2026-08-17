/*
Copyright (C) 2026 Jesaja Weintritt (jesaja.weintritt@stud.eah-jena.de) and 2012 Mark Boots (mark.boots@usask.ca).

This program was originally implemented as a part of the Parallel Efficiency of Gratings project PEG and got reworked in 2026. PEG is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License, version 3, as published by the Free Software Foundation.
See <http://www.gnu.org/licenses/> for details.

This reworked version contains substantial modifications by Jesaja Weintritt (2026) and has not been independently verified against the original. It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; use at your own risk and verify results independently.
*/


#ifndef TESolver_H
#define TESolver_H

#include "PEG.h"

#include <complex>
#include <vector>

#include <Eigen/Dense>

/// Contains the context (memory structures, etc.) and algorithm for solving the grating efficiency.

class TESolver {
public:
	/// Construct a solver context for the given \c grating and math options \c mo.  \c numThreads specifies how many threads to use for fine parallelization; ideally it should be <= the number of processor cores on your computer / on a single cluster node.
	TESolver(const Grating& grating, const MathOptions& mo = MathOptions(), int numThreads = 1, bool measureTiming = false);
	/// Destroy a solver context. All storage is owned by std::vector / Eigen types with automatic memory management, so there is nothing left to free manually.
	~TESolver() = default;

	/// Kept non-copyable to avoid accidentally copying large per-thread buffers/matrices; nothing about the storage itself requires this anymore (unlike before, where copying would have caused a double-free of the GSL-owned memory). Revisit if a genuine use case for copying/moving arises.
	TESolver(const TESolver&) = delete;
	TESolver& operator=(const TESolver&) = delete;
	
	/// Calculates the efficiency at incidence angle \c incidenceDeg and wavelength \c wl.  Side effects: sets the refractive index member variable v_1_; modifies the contents of u_, uprime_, alpha_, beta_, etc.
	Result getEffTE(double incidenceDeg, double wl, double rmsRoughnessNm = 0, bool printDebugOutput = false);


	// Solving implementation functions
	////////////////////////////////////////

	/// Computes alpha_, beta1_, and betaM_ for all n, based on v_1_ (material refractive index) and wl_.  Also calculates
	void computeAlphaAndBeta(double incidenceDeg);

	/// Calculates how many vertical layers (numLayers_ and M_) are sufficient to keep exponentials from contamination.  Fills y_ with the vertical coordinate at each layer.
	void computeLayers();

	/// Computes the blocks of the T matrix (T11_, T12_, T21_, T22_) for the layer below \c y_[m].  Since \c m = 1 is the top of the substrate, \c m can range from [2, M-1].
	Result::Code computeTMatrixBelowLayer(int m, bool printDebugOutput = false);

	/// Calculates the grating fourier expansion for k^2_m at a given \c y value and wavelength \c wl, and stores in \c k2.  \c k2 must have space for 4*N_ + 1 coefficients, since we will be computing from n = -2N_ to 2N.   Reads member variables N_, wavelength wl_, grating refractive index \c v_1_, and grating geometry from \c g_.  Returns Result::Success, or Result::InvalidGratingFailure if the profile is not supported or \c y is larger than the groove height.
	Result::Code computeGratingExpansion(double y, std::complex<double>* k2) const;

	/// Computes the Fourier components of the grating expansion k^2_m into \c k2, based on an array of x crossing (step) values \c stepsX and corresponding k^2 values \c stepsK2 immediately to the left of those x values. \c numSteps is the number of steps [usually two or four, if there are interpenetrating coatings)].  Valid only up to PEG_MAX_PROFILE_CROSSINGS (60) to avoid allocating memory, since this function is called repeatedly.
	void computeGratingExpansion(const double* stepsX, const std::complex<double>* stepsK2, int numSteps, std::complex<double>* k2) const;

	/// Initializes an 8N+4 array of double [\c u, \c uprime] to contain the starting integration values of the electric field Fourier components. The first half of the array \c w contains \c u, the second half contains \c uprime, with each entry in {re,im} order.  The u value is set to $\delta_{n,p}$ and the u' value is set to $-i \beta_n^{(M)} \delta_{n,p}$ or $i \beta_n^{(M)} \delta_{n,p}$, depending on whether p > 2N_.  (Note n,p here are using numbering from 0, and that for the delta functions, p is aliased back onto [0, 2N] once it reaches 2N_+1.
	/*! If layer \c m = 1, then uses beta1_n instead of betaM_n for the derivative. */
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
	void computeBMFromSMatrix();




	// General Mathematical Helper functions:
	///////////////////////////////
	
	/// Returns the square root \c w of a complex number \c z, choosing the branch cut so that Im(w) >= 0.  
	/*! Note: This is different than the usual branch cut defined for the "principal square root", which uses the negative real axis so that Re(w) >= 0 [ie: w is in the right complex plane].*/
	static std::complex<double> complex_sqrt_upperComplexPlane(std::complex<double> z);

	/// Returns the condition number of a complex square matrix \c A (ratio of largest to smallest singular value), computed via Eigen's JacobiSVD.
	static double conditionNumber(const Eigen::MatrixXcd& A);


protected:

	/// Number of threads to use for this calculation
	int numThreads_;
	
	/// The number of Fourier coefficients
	int N_;
	/// 2*N_ + 1, 4*N_+2, and 8*N_+4, since these are used a lot
	int twoNp1_, fourNp2_, eightNp4_;
	
	/// The accuracy to use for numerical integration. Default 1e-8. \todo Get from math options.
	double integrationTolerance_;
	
	
	/// alpha array (size 2N+1)
	std::vector<double> alpha_;
	/// beta array (size 2N+1).  betaM_ is for the superstrate, beta1_ is for the substrate.
	std::vector<std::complex<double>> betaM_, beta1_;
	/// B_n^{M} array: outgoing reflected Rayleigh coefficients. (size 2N+1)
	std::vector<std::complex<double>> BM_;

	/// Number of layers to use in the vertical stack to keep the growing exponentials from numerical contamination.
	int numLayers_;
	/// M-2 = numLayers_.
	int M_;
	
	/// pre-allocated storage for the grating k^2 fourier coefficients.  There is one inner array (per outer vector entry) for each thread to use.  Each inner array has size 2N+1.
	std::vector<std::vector<std::complex<double>>> k2_;
	/// Helper function: returns the k^2 array that should be used by a given thread.
	std::complex<double>* k2ForCurrentThread();
	
	/// This block of storage contains the [u, uprime] electric field Fourier component vectors. They are arranged with each value {re,im}, from order [-N to N], the u vector followed by uprime vector... repeated for each trial solution.  Access the [u, uprime] vector for a given trial solution with wVectorForP().   The size of one w vector is (8*N_+4), and there are (4*N_+2) trial solutions.
	std::vector<std::vector<double>> wVectors_;

	/// Returns the wVector for a trial solution \c p at index \c j, where \c j is numbered from [0, 4*N_+1].
	std::vector<double>& wVectorForP(int p) { return wVectors_[p]; }
	/// Returns the electric field Fourier component \c u for order \c n (index \c i) and trial solution \c p (index \c j), out of wVectors_.  i and j are numbered from 0.
	std::complex<double>* u(int i, int j) {
		return reinterpret_cast<std::complex<double>*>(wVectorForP(j).data() + 2*i);
	}
	/// Returns the electric field Fourier component derivative \c u' for order \c n (index \c i) and trial solution \c p (index \c j), out of wVectors_.  i and j are numbered from 0.
	std::complex<double>* uprime(int i, int j) {
		return reinterpret_cast<std::complex<double>*>(wVectorForP(j).data() + fourNp2_ + 2*i);
	}

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

	/// The y-coordinate of the infinitely-thin Rayleigh layer at y_m, with m = [1, M_ - 1].  y_[0] is unused, so that we can take y_m = y_[m].
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