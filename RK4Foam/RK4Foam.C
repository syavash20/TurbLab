/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    RK4Foam

Description
    Explicit 4th-order Runge-Kutta solver for transient incompressible problems
    
    Implemented by Ehsan Asgari, PhD Candidate at Amirkabir University of 
    Technology (Tehran Polytechnic): eh.asgari@gmail.com

    Acknowledgment:
    
   1- Dr. Ville Vuorinen     
   2- Yu Cheng 
   3- Prof. Filippo Maria Denaro
 

    Sub-models include:
    - turbulence modelling, i.e. laminar, RAS or LES
    - with the capability of assigning a constant source vector

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"
#include "fvIOoptionList.H"
#include "IOporosityModelList.H"
#include "IOMRFZoneList.H"
#include "fixedFluxPressureFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;
    
    #include "createPoissonMatrix.H"
    #include "RK4Coeff.H"

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        runTime++;
        dimensionedScalar dt("dt", dimensionSet(0,0,1,0,0,0,0),SMALL);        
        Info<< "Time = " << runTime.timeName() << nl << endl;
        Uold = U; Uc = U; dt=runTime.deltaT();
        phi = (fvc::interpolate(U) & mesh.Sf());
        dU = dt*(fvc::laplacian(turbulence->nuEff(), U) - fvc::div(phi,U));
        Uc = Uc +b1*dU; U = Uold + 0.5*dU;
        #include "pressureCorrection.H"

        phi = (fvc::interpolate(U) & mesh.Sf());
        dU = dt*(fvc::laplacian(turbulence->nuEff(), U) - fvc::div(phi,U));
        Uc = Uc +b2*dU;
        U = Uold + 0.5*dU;
        #include "pressureCorrection.H"

        phi = (fvc::interpolate(U) & mesh.Sf());
        dU = dt*(fvc::laplacian(turbulence->nuEff(), U) - fvc::div(phi,U));
        Uc = Uc +b3*dU; U = Uold + dU;
        #include "pressureCorrection.H"


        phi = (fvc::interpolate(U) & mesh.Sf());
        dU = dt*(fvc::laplacian(turbulence->nuEff(), U) - fvc::div(phi,U));
        Uc = Uc +b4*dU; U = Uc;
        #include "pressureCorrection.H"
        U = U + deltaP * dt;
        U.correctBoundaryConditions();
        phi = (fvc::interpolate(U) & mesh.Sf());
        turbulence->correct();
  
        runTime.write();        

      
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
