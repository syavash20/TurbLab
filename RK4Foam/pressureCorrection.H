#include "continuityErrs.H"
U.correctBoundaryConditions();
solve(pEqn == fvc::div(U)/runTime.deltaT());
U = U - (fvc::grad(p) * runTime.deltaT());
U.correctBoundaryConditions();

