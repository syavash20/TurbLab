/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

Class
    Foam::LESModels::DTCSM

Group
    grpLESTurbulence

Description
    Non-Boussinesq subgrid-scale model with dynamic procedure

    Reference:
    \verbatim
        Agrawal et al. (2022).
        Non-Boussinesq subgrid-scale model with dynamic tensorial coefficients
        Physical Review Fluids 7, 074602 (2022).
    \endverbatim
    Implemented by Ehsan Asgari, PhD
SourceFiles
    DTCSM.C

\*---------------------------------------------------------------------------*/

#ifndef DTCSM_H
#define DTCSM_H

#include "LESModel.H"
#include "LESeddyViscosity.H"
#include "ReynoldsStress.H"
#include "simpleFilter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                    Class DeardorffDiffStress Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class DTCSM
:
    public ReynoldsStress<LESModel<BasicTurbulenceModel>>
{
    // Private Member Functions

        DTCSM(const DTCSM&) = delete;

        //- No copy assignment
        void operator=(const DTCSM&) = delete;

protected:

    // Protected data

        // Model constants

            dimensionedScalar Ck_;
            dimensionedScalar Cm_;
            dimensionedScalar Ce_;
            dimensionedScalar Cs_;

	    volTensorField C_M_;
//            volScalarField C_12_;
//            volScalarField C_13_;
//            volScalarField C_23_;

            simpleFilter simpleFilter_;
            autoPtr<LESfilter> filterPtr_;
            LESfilter& filter_;


    // Protected Member Functions

        //- Update the eddy-viscosity
        virtual void correctNut();


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("DTCSM");


    // Constructors

        //- Constructor from components
        DTCSM
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~DTCSM() = default;
    


    // Member Functions

        //- Read model coefficients if they have changed
        virtual bool read();

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const;

        //- Return the specific dissipation rate
        virtual tmp<volScalarField> omega() const;

        //- Correct sub-grid stress, eddy-Viscosity and related properties
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "DTCSM.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
