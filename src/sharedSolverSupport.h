/*
Copyright (C) 2026 Jesaja Weintritt (jesaja.weintritt@stud.eah-jena.de) and 2012 Mark Boots (mark.boots@usask.ca).

This program was originally implemented as a part of the Parallel Efficiency of Gratings project PEG and got reworked in 2026. PEG is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License, version 3, as published by the Free Software Foundation.
See <http://www.gnu.org/licenses/> for details.

This reworked version contains substantial modifications by Jesaja Weintritt (2026) and has not been independently verified against the original. It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; use at your own risk and verify results independently.
*/

// Provides for both solvers functions for solving, due their partial similar calculation steps

/*
Identical in TESolver and TMSolver:

- computeAlphaAndBeta, computeLayers and computeBMFromSMatrix
- both computeGratingExpansion overloads, and complex_sqrt_upperComplexPlane
- conditionNumber, which is never called
- GratingExpansionError, the matrix and vector allocations, and the T-matrix assembly loop
- most of getEffTE/getEffTM (99% the same; only the debug transmission line differs)
- the unused odeJacobian comment block
*/