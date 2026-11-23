/*
Copyright (C) 2026 Jesaja Weintritt (jesaja.weintritt@stud.eah-jena.de) and 2012 Mark Boots (mark.boots@usask.ca).

This program was originally implemented as a part of the Parallel Efficiency of Gratings project PEG and got reworked in 2026. PEG is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License, version 3, as published by the Free Software Foundation.
See <http://www.gnu.org/licenses/> for details.

This reworked version contains substantial modifications by Jesaja Weintritt (2026) and has not been independently verified against the original. It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; use at your own risk and verify results independently.
*/


#include "TESolver.h"

double TESolver::transmittedOrderWeight(int i) const
{
	return beta1_[i].real();
}

void TESolver::columnDerivatives(const Eigen::MatrixXcd& U, const Eigen::MatrixXcd& W, Eigen::MatrixXcd& dU, Eigen::MatrixXcd& dW) const
{
	/// d(u)/dy = u', and u' is the second field itself.
	dU = W;

	/// d(u')/dy = (diag(alpha_n^2) - Toeplitz(k^2(y))) u for all columns of this block at once: one matrix product instead of one matrix-vector product per column. toeplitz_ was built once for this y by SolverBase::odeFunctionMatrix(). (The Jacobian of this system is, per column, [[0, I], [diag(alpha_n^2) - toeplitz_, 0]]; it is not needed by bulirsch_stoer. The old per-column odeJacobian() is in the git history.)
	dW.noalias() = toeplitz_ * U;
	for(int i=0; i<twoNp1_; ++i)
		dW.row(i) = (alpha_[i] * alpha_[i]) * U.row(i) - dW.row(i);
}