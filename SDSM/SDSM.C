/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "SDSM.H"
//#include "fvOptions.H"
//#include "fvOptionList.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void SDSM<BasicTurbulenceModel>::correctNut
(
    const tmp<volTensorField>& gradU
)
{
    this->nut_ = max(Cs_*sqr(this->delta())*mag(dev(symm(gradU))),-1.0*this->nu());
    this->nut_.correctBoundaryConditions();
    //fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


template<class BasicTurbulenceModel>
void SDSM<BasicTurbulenceModel>::correctNut()
{
    correctNut(fvc::grad(this->U_));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
SDSM<BasicTurbulenceModel>::SDSM
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    LESeddyViscosity<BasicTurbulenceModel>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    k_
    (
        IOobject
        (
            IOobject::groupName("k", this->U_.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),

    Cs_
    (
        IOobject
        (
            IOobject::groupName("Cs", this->U_.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar ("Cs", dimless,SMALL)
    ),

    invScSgs_
    (
        IOobject
        (
            IOobject::groupName("invScSgs", this->U_.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar ("invScSgs", dimless,SMALL)
    ),

    simpleFilter_(U.mesh()),
    filterPtr_(LESfilter::New(U.mesh(), this->coeffDict())),
    filter_(filterPtr_())
{
    bound(k_, this->kMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool SDSM<BasicTurbulenceModel>::read()
{
    if (LESeddyViscosity<BasicTurbulenceModel>::read())
    {
        filter_.read(this->coeffDict());        

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
void SDSM<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
//    const surfaceScalarField& phi = this->phi_;
    const volVectorField& U = this->U_;
    //fv::options& fvOptions(fv::options::New(this->mesh_));

    LESeddyViscosity<BasicTurbulenceModel>::correct();

    tmp<volTensorField> tgradU(fvc::grad(U));
    const volTensorField& gradU = tgradU();

    volSymmTensorField D(dev(symm(gradU)));
    volScalarField magD(mag(D));

    volVectorField Uf(filter_(U));

    volSymmTensorField Df(filter_(D));  
//    volSymmTensorField Sf(dev(symm(fvc::grad(Uf))));
    
    volScalarField magDf(mag(Df));
          
    
    volSymmTensorField LL =
    (dev(filter_(sqr(U)) - (sqr(filter_(U)))));

    volSymmTensorField MM
    (
        sqr(this->delta())*(filter_(magD*D) - 4.0*magDf*Df)
    );
    
    volScalarField MMMM = fvc::average(magSqr(MM));
    MMMM.max(VSMALL);

    Cs_ = 0.5 * fvc::average(LL && MM)/MMMM;

    IOdictionary transportProperties
    (
        IOobject
	(
            "transportProperties",
            this->runTime_.constant(),
            this->mesh_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
     );

word name;
word tmp(transportProperties.lookupOrDefault<word>("scalar", "none"));
name = tmp;

if (name != "none")
  {
    Info<< "Reading field " << name << endl;
    const volScalarField& S = this->mesh_.objectRegistry::template lookupObject<volScalarField>(name);    

    volVectorField F = filter_(S*U) - (filter_(S)*filter_(U));
    
    volVectorField H = sqr(this->delta())*(-4.0*magDf*fvc::grad(filter_(S)) + filter_(magD*fvc::grad(S)));
    
    volScalarField HH = fvc::average(magSqr(H));
    HH.max(VSMALL);

    volScalarField C_S = fvc::average(F & H)/HH;

    invScSgs_ =  C_S / (2.0 * Cs_ + VSMALL);
  }  
else
  {
    invScSgs_ = 0;
  }  
    this->invScSgs_.correctBoundaryConditions();


    volScalarField KK =
    0.5*(filter_(magSqr(U)) - magSqr(filter_(U)));
    
    volScalarField mm
    (
        sqr(this->delta())*(4.0*sqr(mag(Df)) - filter_(sqr(magD)))
       
    );

    volScalarField mmmm = fvc::average(magSqr(mm));
    mmmm.max(VSMALL);

    k_ = fvc::average(KK*mm)/mmmm * sqr(this->delta())*magSqr(D);

    correctNut(gradU);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
