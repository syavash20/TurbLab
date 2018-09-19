/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    postLES

Description
    Calculates and writes some mean fields from an LES calculation. The
    following fields are created:
        - resLES: resolutness of an LES
        - TKEMean: mean turbulent kinectic energy
        - TKEMeanProd: production term of the mean turbulent kinectic energy
		- turbDiffusionMean: turbulent diffusion term of the mean turbulent kinectic energy
		- SGSDiffusionMean: subgrid scale diffusion term of the mean turbulent kinectic energy
		- viscDiffusionMean: viscous diffusion term of the mean turbulent kinectic energy
    Fields UMean, UPrime2Mean, turbDiffMean, SGSDiffMean, and kMean must exist.

\*---------------------------------------------------------------------------*/

#include "calc.H"
#include "fvc.H"
#include "fvCFD.H"
#include "incompressible/singlePhaseTransportModel/singlePhaseTransportModel.H"
#include "incompressible/turbulenceModel/turbulenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::calc
(
const argList& args, 
const Time& runTime, 
const fvMesh& mesh 
)


{
       
    IOobject UPrime2MeanHeader
    (
        "UPrime2Mean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );

    IOobject kMeanHeader
    (
        "kMean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );

    IOobject UMeanHeader
    (
        "UMean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );
          
    IOobject turbDiffMeanHeader
    (
        "turbDiffMean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );         
    
    IOobject SGSDiffMeanHeader
    (
        "SGSDiffMean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );     
    if
    (
        UPrime2MeanHeader.headerOk()
     && kMeanHeader.headerOk()
     && UMeanHeader.headerOk()
     && turbDiffMeanHeader.headerOk()
     && SGSDiffMeanHeader.headerOk()     
    )
    
    {
        Info<< "    Reading average field UPrime2Mean" << endl;
        const volSymmTensorField UPrime2Mean(UPrime2MeanHeader, mesh);

        Info<< "    Reading average field kMean" << endl;
        const volScalarField kMean(kMeanHeader, mesh);

        Info<< "    Reading average field UMean" << endl;
        const volVectorField UMean(UMeanHeader, mesh);                          
        
        Info<< "    Reading average field turbDiffMean" << endl;
        const volVectorField turbDiffMean(turbDiffMeanHeader, mesh);                                  
          
        Info<< "    Reading average field SGSDiffMean" << endl;
        const volVectorField SGSDiffMean(SGSDiffMeanHeader, mesh); 
        
        volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);                 


///////////////////////////////////////////////////////////////////////        
        
        volScalarField TKEMean
        (
            IOobject
            (
                "TKEMean",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("dimkMean",kMean.dimensions(),0)
        );

        volScalarField resLES
        (
            IOobject
            (
                "resLES",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("dimless",dimless,0)
        );

        volScalarField TKEMeanProd
        (
            IOobject
            (
                "TKEMeanProd",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("dimTKEMeanProd",kMean.dimensions()/dimTime,0)
        );
        
        volScalarField turbDiffusionMean
        (
            IOobject
            (
                "turbDiffusionMean",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("turbDiffusionMean",kMean.dimensions()/dimTime,0)
        );        
        
        volScalarField SGSDiffusionMean
        (
            IOobject
            (
                "SGSDiffusionMean",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("SGSDiffusionMean",kMean.dimensions()/dimTime,0)
        );
        
        volScalarField viscDiffusionMean
        (
            IOobject
            (
                "viscDiffusionMean",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("viscDiffusionMean",kMean.dimensions()/dimTime,0)
        ); 
        #include "createPhi.H"  
        //~ #include "createFields.H" 
        
        singlePhaseTransportModel laminarTransport(U, phi);

        autoPtr<incompressible::turbulenceModel> turbulence
        (
                incompressible::turbulenceModel::New(U, phi, laminarTransport)
        );
        
        const volScalarField nuLam(turbulence->nu());                             

        Info<< "    Calculating mean turbulent kinetic energy TKEMean" << endl;
        TKEMean = 0.5 * (
                            UPrime2Mean.component(0)
                          + UPrime2Mean.component(3)
                          + UPrime2Mean.component(5)
                        );

        Info<< "    Calculating Mean LES resolutness resLES" << endl;
        forAll(resLES, i)
        {
            resLES[i] = TKEMean[i] / (TKEMean[i] + kMean[i]);
        }

        Info<< "    Calculating Mean TKE production TKEMeanProd" << endl;
        TKEMeanProd = - UPrime2Mean && fvc::grad(UMean);
        
        Info<< "    Calculating turbDiffusionMean" << endl;
        turbDiffusionMean = fvc::div(turbDiffMean);

        Info<< "    Calculating SGSDiffusionMean" << endl;
        SGSDiffusionMean = fvc::div(SGSDiffMean);        
        
        Info<< "    Calculating viscDiffusionMean" << endl;
        viscDiffusionMean = nuLam*fvc::laplacian(TKEMean);         
        
        resLES.write();
        TKEMean.write();
        TKEMeanProd.write(); 
        turbDiffusionMean.write();        
        SGSDiffusionMean.write();
        viscDiffusionMean.write();                
               
    }
    else
    {
        Info << "    No UPrime2Mean and/or kMean and/or UMean." << endl;
    }
    
    Info<< "End" << endl;
}


// ************************************************************************* //
