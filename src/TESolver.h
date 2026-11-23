/*
Copyright (C) 2026 Jesaja Weintritt (jesaja.weintritt@stud.eah-jena.de) and 2012 Mark Boots (mark.boots@usask.ca).

This program was originally implemented as a part of the Parallel Efficiency of Gratings project PEG and got reworked in 2026. PEG is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License, version 3, as published by the Free Software Foundation.
See <http://www.gnu.org/licenses/> for details.

This reworked version contains substantial modifications by Jesaja Weintritt (2026) and has not been independently verified against the original. It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; use at your own risk and verify results independently.
*/


#ifndef TESolver_H
#define TESolver_H

#include "SolverBase.h"

class TESolver : public SolverBase {
public:
	/// Same constructor as SolverBase: a solver context for the given \c grating and math options \c mo. \c numThreads specifies how many threads to use for fine parallelization.
	using SolverBase::SolverBase;

protected:
	/// Weight of substrate order \c i in the debug transmission sum (see SolverBase::transmittedOrderWeight()): Re(beta1_[i]).
	double transmittedOrderWeight(int i) const override;

	const char* secondFieldName() const override { return "u'"; }
	bool secondFieldIsEpsScaled() const override { return false; }

	/// The TE equations for a block of trial columns: dU = W (= u'), dW = (diag(alpha_n^2) - toeplitz_) * U.
	void columnDerivatives(const Eigen::MatrixXcd& U, const Eigen::MatrixXcd& W, Eigen::MatrixXcd& dU, Eigen::MatrixXcd& dW) const override;
};

#endif