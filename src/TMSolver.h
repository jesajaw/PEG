/*
Copyright (C) 2026 Jesaja Weintritt (jesaja.weintritt@stud.eah-jena.de) and 2012 Mark Boots (mark.boots@usask.ca).

This program was originally implemented as a part of the Parallel Efficiency of Gratings project PEG and got reworked in 2026. PEG is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License, version 3, as published by the Free Software Foundation.
See <http://www.gnu.org/licenses/> for details.

This reworked version contains substantial modifications by Jesaja Weintritt (2026) and has not been independently verified against the original. It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; use at your own risk and verify results independently.
*/


#ifndef TMSolver_H
#define TMSolver_H

#include "SolverBase.h"

class TMSolver : public SolverBase {
public:
	/// Construct a solver context for the given \c grating and math options \c mo. \c numThreads specifies how many threads to use for fine parallelization
	TMSolver(const Grating& grating, const MathOptions& mo = MathOptions(), int numThreads = 1, bool measureTiming = false);

	/// Calculates the grating fourier expansion of 1/eps_r(y) at a given \c y value (same step geometry as computeGratingExpansion(), but expanding 1/eps_r instead of k^2). \c invEps must have space for 2*N_ + 1 coefficients.
	Result::Code computeInvEpsExpansion(double y, std::complex<double>* invEps) const;

protected:
	/// Weight of substrate order \c i in the debug transmission sum (see SolverBase::transmittedOrderWeight()).
	double transmittedOrderWeight(int i) const override;

	const char* secondFieldName() const override { return "w"; }
	bool secondFieldIsEpsScaled() const override { return true; }

	/// Fills k2Inv_ and toeplitzInv_ (the Toeplitz matrix of 1/eps_r(y)) for this y.
	void prepareAtY(double y) override;

	/// Factorizes the two Toeplitz matrices (LU), in parallel as two OpenMP sections. Called by every thread of the team in SolverBase::odeFunctionMatrix().
	void factorizeAtY() override;

	/// The TM equations for a block of trial columns, reusing the ONE LU factorization per y for all columns: dU = Toeplitz(1/eps_r)^{-1} W, dW = k_0^2 ( D_alpha * Toeplitz(k^2)^{-1} D_alpha U - U ) (see the class-level comment above).
	void columnDerivatives(const Eigen::MatrixXcd& U, const Eigen::MatrixXcd& W, Eigen::MatrixXcd& dU, Eigen::MatrixXcd& dW) const override;

	/// Fourier coefficients of 1/eps_r(y) (size 2N_+1), the Toeplitz matrix built from them, and the LU factorizations of both Toeplitz matrices (refilled at every y).
	std::vector<std::complex<double>> k2Inv_;
	Eigen::MatrixXcd toeplitzInv_;
	Eigen::PartialPivLU<Eigen::MatrixXcd> luInvEps_, lu_;
};
#endif