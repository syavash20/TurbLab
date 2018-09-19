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
    LESpimpleFoam

Description
    
    A modified pimpleFoam solver for calculating of TKE budjets
    
Ref. for mathematics:
* 1) Saeedi M, Wang B-C. Large-eddy simulation of turbulent flow and dispersion over a
matrix of wall-mounted cubes. Phys. Fluids. 2015; 27:115104
 

* 2) Lectures in Turbulence for the 21st Century By William K. George     

Regards,
Ehsan Asgari, PhD Candidate @ Amirkabir University of Technology (Tehran Polytechnic)
e.asgari@aut.ac.ir
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "incompressible/turbulenceModel/turbulenceModel.H"
#include "incompressible/LES/LESModel/LESModel.H"
#include "pimpleControl.H"
#include "fvIOoptionList.H"
#include "IOporosityModelList.H"
#include "IOMRFZoneList.H"
#include "fixedFluxPressureFvPatchScalarField.H"
//~ #include "AverageIOField.H"
//~ #include "AverageIOField2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"

    
    pimpleControl pimple(mesh);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();               
            }
        }
        
    const volScalarField nuLam(turbulence->nu());        
                
    const volSymmTensorField G(turbulence->R());
    
    const volScalarField nu_SGS(SGSModel->nuSgs());

            
            UPrime = (U-UMean);         //Resolved Velocity Fluctuation Vector
            pPrime = (p-pMean);              //Resolved Pressure Fluctuation
            B = -2.0*nu_SGS*symm(fvc::grad(U));
            turbDiff = -0.5*(UPrime*magSqr(UPrime));      //Turbulent Diffusion Term--divergence operator will be applied afterwards
            pressDiff = -UPrime & fvc::grad(pPrime);           //Pressure Diffusion Term
            SGSstrainTensor = symm(fvc::grad(UPrime));       //Tensor of Strain Rate of Resolved Fluctuations
            viscDiss = -2*nuLam*(SGSstrainTensor && SGSstrainTensor); //Viscous Dissipation of Resolved Fluctuations
            SGSDiff = -UPrime & G;                        //SGS Diffusion Term--divergence operator will be applied afterwards
            SGSDiss = B && SGSstrainTensor;                  //SGS Dissipation Term
            runTime.write();
        
       
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
