// Solve the Momentum equation
if (name != "none")
{	

dimensionedScalar alpha = nu / Sc;
const volScalarField& invSct = mesh.objectRegistry::template lookupObject<volScalarField>("invScSgs");
const volScalarField& nut = mesh.objectRegistry::template lookupObject<volScalarField>("nut");

fvScalarMatrix SEqn	    
(
    fvm::ddt(S) + fvm::div(phi, S)
  - fvm::laplacian(alpha+(nut*invSct), S)
);

SEqn.relax();
SEqn.solve();

}
