Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field UMean\n" << endl;
volVectorField UMean
(
    IOobject
    (
        "UMean",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh
);

Info<< "Reading field pMean\n" << endl;
volScalarField pMean
(
    IOobject
    (
        "pMean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

volVectorField UPrime
(
    IOobject
    (
        "UPrime",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
           dimensionedVector
      (
           "UPrime", 
           dimVelocity, // Or Umean.dimensions()
           vector (0,0,0)
      )
);
//~ 
volScalarField pPrime
(
    IOobject
    (
        "pPrime",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedScalar
      (
           "pPrime", 
           dimVelocity*dimVelocity, // Or Umean.dimensions()
           scalar(0.0)
      )    
);
//~ 

volSymmTensorField B
(
    IOobject
    (
        "B",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
           dimensionedSymmTensor
      (
           "B", 
           dimVelocity*dimVelocity,
           symmTensor::zero
      )
);


volSymmTensorField SGSstrainTensor
(
    IOobject
    (
        "SGSstrainTensor",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
           dimensionedSymmTensor
      (
           "SGSstrainTensor", 
           dimVelocity/dimLength,
           symmTensor::zero
      )
);

volScalarField viscDiss
(
    IOobject
    (
        "viscDiss",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedScalar
      (
           "viscDiss", 
           dimVelocity*dimVelocity*dimVelocity/dimLength,
           0.0
      )
);

volScalarField pressDiff
(
    IOobject
    (
        "pressDiff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedScalar
      (
           "pressDiff", 
           dimVelocity*dimLength/(dimTime*dimTime),
           0.0
      )
);

volVectorField turbDiff
(
    IOobject
    (
        "turbDiff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedVector
      (
           "turbDiff", 
           dimVelocity*dimVelocity*dimVelocity, // Or Umean.dimensions()
           vector (0,0,0)
      )
);

volScalarField SGSDiss
(
    IOobject
    (
        "SGSDiss",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedScalar
      (
           "SGSDiss", 
           dimVelocity*dimVelocity*dimVelocity/dimLength,
           0.0
      )
);

volVectorField SGSDiff
(
    IOobject
    (
        "SGSDiff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    
           dimensionedVector
      (
           "SGSDiff", 
           dimVelocity*dimVelocity*dimVelocity, // Or Umean.dimensions()
           vector (0,0,0)
      )
);
 
 
#include "createPhi.H"


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, mesh.solutionDict().subDict("PIMPLE"), pRefCell, pRefValue);


singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

autoPtr<incompressible::LESModel> SGSModel
(
incompressible::LESModel::New(U, phi, laminarTransport)
);
