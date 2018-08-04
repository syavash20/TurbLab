/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2013-12-28 Oscar Ochoa: Adapted OpenFOAM's 2.1.1 wallShearStress to LES.
2014-06-22 Bruno Santos: Adapted to OpenFOAM 2.2.x.
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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
    integralThickness

Description
    Calculates and reports boundary-layer integral thickness 
	including theta (momentun thickness) and deltaStar (displacement thickness).
	
	The streamwise and spanwise directions of the flow must be X and Z, respectively.
	A time-averaged velocity field (Uvec) is required for the code to work.	
	This code has been developed for a fully hexahedral mesh but it may work for a non-hex mesh either.
	

	Ehsan Asgari (syavash20), e.asgari@aut.ac.ir	
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "incompressible/singlePhaseTransportModel/singlePhaseTransportModel.H"
#include "incompressible/LES/LESModel/LESModel.H"
#include "incompressible/turbulenceModel/turbulenceModel.H"
#include "LES/LESdeltas/LESdelta/LESdelta.H"
#include "compressible/LES/LESModel/LESModel.H"
#include "nearWallDist.H"
#include "wallDist.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void calcIncompressible
(
    const fvMesh& mesh,
    const Time& runTime,
    const volVectorField& U,    
    volScalarField& theta,
    volScalarField& deltaStar
)
{
    #include "createPhi.H"

    singlePhaseTransportModel laminarTransport(U, phi);

    autoPtr<incompressible::LESModel> model
    (
        incompressible::LESModel::New(U, phi, laminarTransport)
    );

  

volVectorField cellDelta_
        (
            IOobject
            (
                "cellDelta",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedVector("cellDelta", dimless, vector::zero)
        );

      const faceList & ff = mesh.faces();
      const pointField & pp = mesh.points();
     
      forAll (cellDelta_, celli)
      {
       const cell & cc = mesh.cells()[celli];
       labelList pLabels(cc.labels(ff));
       pointField pLocal(pLabels.size(), vector::zero);

           forAll (pLabels, pointi)
           {
           pLocal[pointi] = pp[pLabels[pointi]];
           }
        cellDelta_[celli].component(vector::X) = Foam::max(pLocal & vector(1,0,0)) - Foam::min(pLocal & vector(1,0,0));
        cellDelta_[celli].component(vector::Y) = Foam::max(pLocal & vector(0,1,0)) - Foam::min(pLocal & vector(0,1,0));
        cellDelta_[celli].component(vector::Z) = Foam::max(pLocal & vector(0,0,1)) - Foam::min(pLocal & vector(0,0,1));
      }   
 


   const fvPatchList& patches = mesh.boundary();

scalar xCor = 0.0;
scalar zCor = 0.0;
label n = 0;
label flag = 0;

forAll(mesh.boundary(), patchID) 
{

        const fvPatch& currPatch = patches[patchID];
        scalarField xGlobal(currPatch.size(), 100);


            if (isA<wallFvPatch>(currPatch))
            {

              forAll (mesh.boundary()[patchID],facei) 
              {

               flag = 0;
               xCor = mesh.Cf().boundaryField()[patchID][facei].component(0);
               zCor = mesh.Cf().boundaryField()[patchID][facei].component(2);
               
               forAll (mesh.boundary()[patchID],faceii)
               {
               if (mag(xCor - xGlobal[faceii]) < 1e-9)
               {
               theta.boundaryField()[patchID][facei]=theta.boundaryField()[patchID][faceii];
               deltaStar.boundaryField()[patchID][facei]=deltaStar.boundaryField()[patchID][faceii];
               flag = 1;
// Info<< "xGlobal " << xGlobal[facei]
//     << " zCor " << zCor
//     << " facei " << faceii
//     << " xCor " << xCor << nl << endl;
              
               break;              
               }
               }


              if (flag == 0)
              {
              forAll (mesh.C(), celli)
              {            
              if (mag(xCor - mesh.C()[celli].component(0)) < 0.5*cellDelta_[celli].component(vector::X) && mag(zCor - mesh.C()[celli].component(2)) < 0.5*cellDelta_[celli].component(vector::Z))
              {
              
              theta.boundaryField()[patchID][facei]+=cellDelta_[celli].component(vector::Y)*(U[celli].component(0)*(1.0-U[celli].component(0)));
              deltaStar.boundaryField()[patchID][facei]+=cellDelta_[celli].component(vector::Y)*(1.0-U[celli].component(0));                     
              }             
              }
               xGlobal[facei] = xCor;       
              }

//if (facei == 63)
//{
//Info << " facei " << facei
//     << "flag" << flag
//     << " xCor " << xCor << nl << endl;
//}

              }
                                        
            }
                const scalarField& theta_wall = theta.boundaryField()[patchID];

                Info<< "Patch " << patchID
                    << " named " << currPatch.name()
                    << " theta : min: " << gMin(theta_wall) << " max: " << gMax(theta_wall)
                    << " average: " << gAverage(theta_wall) << nl << endl;
}
                    

}


int main(int argc, char *argv[])
{
    timeSelector::addOptions();

    #include "addRegionOption.H"

    argList::addBoolOption
    (
        "compressible",
        "calculate compressible wall shear stress"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createNamedMesh.H"

    const bool compressible = args.optionFound("compressible");

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;
        mesh.readUpdate();

       
        volScalarField theta
        (
            IOobject
            (
                "theta",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("theta", dimless, 0.0)
        );
        
        volScalarField deltaStar
        (
            IOobject
            (
                "deltaStar",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("deltaStar", dimless, 0.0)
        );

        IOobject UHeader
        (
            "Uvec",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        if (UHeader.headerOk())
        {
            Info<< "Reading field Uvec\n" << endl;
            volVectorField U(UHeader, mesh);

            if (compressible)
            {
 //               calcCompressible(mesh, runTime, U);
            }
            else
            {
                calcIncompressible(mesh, runTime, U, theta, deltaStar);
            }
        }
        else
        {
            Info<< "    no Uvec field" << endl;
        }

        Info<< "Writing values to field " << theta.name()
            << nl << endl;

      
        theta.write();
        deltaStar.write();

    }
    

    Info<< "End" << endl;

    return 0;
}


// ************************************************************************* //
