fvScalarMatrix pEqn
(
fvm::laplacian(p)
);
pEqn.setReference(pRefCell, pRefValue);
