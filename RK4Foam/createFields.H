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

volVectorField Uold
(
IOobject
(
"Uold",
runTime.timeName(),
mesh,
IOobject::NO_READ,
IOobject::NO_WRITE
),
U
);

volVectorField Uc
(
IOobject
(
"Uc",
runTime.timeName(),
mesh,
IOobject::NO_READ,
IOobject::NO_WRITE
),
U
);

volVectorField dU
(
IOobject
(
"dU",
runTime.timeName(),
mesh,
IOobject::NO_READ,
IOobject::NO_WRITE
),
U
);

#include "createPhi.H"


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, mesh.solutionDict().subDict("RK4"), pRefCell, pRefValue);


singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );


Info<< "Reading explicit source from transportProperties\n" << endl;

   dimensionedScalar delPx
    (
       dimensionedScalar::lookupOrAddToDict("delPx",transportProperties,0)
    );

   dimensionedScalar delPy
    (
       dimensionedScalar::lookupOrAddToDict("delPy",transportProperties,0)
        
    );
   dimensionedScalar delPz
    (
       dimensionedScalar::lookupOrAddToDict("delPz",transportProperties,0)
        
    );

Info << "\n delPx = " <<delPx.value()<< "\n delPy = " <<delPy.value()<<"\n delPz = " <<delPz.value()<< endl;

volVectorField deltaP
(
    IOobject
    (
        "deltaP",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("deltaP", dimVelocity/dimTime, vector(delPx.value(),delPy.value(),delPz.value()))
);




