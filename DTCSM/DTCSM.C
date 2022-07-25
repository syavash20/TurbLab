/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "DTCSM.H"
#include "fvOptions.H"
#include "scalarMatrices.H"
#include "MatrixTools.H"
#include "LUscalarMatrix.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void DTCSM<BasicTurbulenceModel>::correctNut()
{
    this->nut_ = Ck_*sqrt(mag(this->k()))*this->delta();
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
DTCSM<BasicTurbulenceModel>::DTCSM
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
    ReynoldsStress<LESModel<BasicTurbulenceModel>>
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

    Ck_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Ck",
            this->coeffDict_,
            0.094
        )
    ),

    Cm_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cm",
            this->coeffDict_,
            4.13
        )
    ),
    Ce_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Ce",
            this->coeffDict_,
            1.05
        )
    ),
    Cs_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cs",
            this->coeffDict_,
            0.25
        )
    ),

    C_M_
    (
        IOobject
        (
            IOobject::groupName("C_M", this->U_.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedTensor ("C_M", dimless,tensor(0,0,0,0,0,0,0,0,0))
    ),

    simpleFilter_(U.mesh()),
    filterPtr_(LESfilter::New(U.mesh(), this->coeffDict())),
    filter_(filterPtr_())

{
    if (type == typeName)
    {
        this->printCoeffs(type);
        this->boundNormalStress(this->R_);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool DTCSM<BasicTurbulenceModel>::read()
{
    if (ReynoldsStress<LESModel<BasicTurbulenceModel>>::read())
    {
        Ck_.readIfPresent(this->coeffDict());
        Cm_.readIfPresent(this->coeffDict());
        Ce_.readIfPresent(this->coeffDict());
        Cs_.readIfPresent(this->coeffDict());
        filter_.read(this->coeffDict());

        return true;
    }

    return false;
}

template<class BasicTurbulenceModel>
tmp<volScalarField> DTCSM<BasicTurbulenceModel>::epsilon() const
{
    volScalarField k(this->k());

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("epsilon", this->alphaRhoPhi_.group()),
                this->runTime_.timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->Ce_*k*sqrt(k)/this->delta()
        )
    );
}

template<class BasicTurbulenceModel>
tmp<volScalarField> DTCSM<BasicTurbulenceModel>::omega() const
{
    volScalarField k(this->k());
    volScalarField epsilon(this->Ce_*k*sqrt(k)/this->delta());

    return tmp<volScalarField>::New
    (
        IOobject
        (
            IOobject::groupName("omega", this->alphaRhoPhi_.group()),
            this->runTime_.timeName(),
            this->mesh_
        ),
        epsilon/(0.09*k)
    );
}


template<class BasicTurbulenceModel>
void DTCSM<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volSymmTensorField& R = this->R_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    ReynoldsStress<LESModel<BasicTurbulenceModel>>::correct();

    tmp<volTensorField> tgradU(fvc::grad(U));
    const volTensorField& gradU = tgradU();

    volSymmTensorField D(symm(gradU));

    volSymmTensorField P(-twoSymm(R & gradU));

    volScalarField k(this->k());

    volSymmTensorField alpha2 = sqr(2.0*this->delta()/this->delta()) * sqrt(2.0) * mag(filter_(D)) * filter_(D);

    volSymmTensorField beta = sqrt(2.0) * mag(D) * D;

    volSymmTensorField L =
        filter_(sqr(U)) - sqr(filter_(U));

   volSymmTensorField M
    (
        alpha2 - filter_(beta)
    );

     forAll(M.ref(), celli)
     {

     scalarRectangularMatrix A(6, 4, Zero);

        A[0][0] = 2.0 * M[celli].component(symmTensor::XX);

        A[0][1] = 2.0 * M[celli].component(symmTensor::XY);

        A[0][2] = 2.0 * M[celli].component(symmTensor::XZ);

        A[0][3] = 0.0;

        A[1][0] = 2.0 * M[celli].component(symmTensor::YY);

        A[1][1] = -2.0 * M[celli].component(symmTensor::XY);

        A[1][2] = 0.0;

        A[1][3] = 2.0 * M[celli].component(symmTensor::YZ);

        A[2][0] = 2.0 * M[celli].component(symmTensor::ZZ);

        A[2][1] = 0.0;

        A[2][2] = -2.0 * M[celli].component(symmTensor::XZ);

        A[2][3] = -2.0 * M[celli].component(symmTensor::YZ);

        A[3][0] = 2.0 * M[celli].component(symmTensor::XY);

        A[3][1] = M[celli].component(symmTensor::YY) - M[celli].component(symmTensor::XX);

        A[3][2] = M[celli].component(symmTensor::YZ);

        A[3][3] = M[celli].component(symmTensor::XZ);

        A[4][0] = 2.0 * M[celli].component(symmTensor::XZ);

        A[4][1] = M[celli].component(symmTensor::YZ);

        A[4][2] = M[celli].component(symmTensor::ZZ) - M[celli].component(symmTensor::XX);

        A[4][3] = -M[celli].component(symmTensor::XY);

        A[5][0] = 2.0 * M[celli].component(symmTensor::YZ);

        A[5][1] = -M[celli].component(symmTensor::XZ);

        A[5][2] = -M[celli].component(symmTensor::XY);

        A[5][3] = M[celli].component(symmTensor::ZZ) - M[celli].component(symmTensor::YY);


     scalarRectangularMatrix Coeff(6, 1, 0);
        Coeff[0][0] = L[celli].component(symmTensor::XX);

        Coeff[1][0] = L[celli].component(symmTensor::YY);

        Coeff[2][0] = L[celli].component(symmTensor::ZZ);

        Coeff[3][0] = L[celli].component(symmTensor::XY);

        Coeff[4][0] = L[celli].component(symmTensor::XZ);

        Coeff[5][0] = L[celli].component(symmTensor::YZ);

     scalarRectangularMatrix LHS(4, 4, Zero);
     scalarRectangularMatrix RHS(4, 1, 0);

	LHS = A.T() * A;

     scalarSquareMatrix LHS2(LHS);
	RHS = A.T() * Coeff;
     scalarDiagonalMatrix RHS2(4, 0);
     RHS2[0] = RHS[0][0];
     RHS2[1] = RHS[1][0];
     RHS2[2] = RHS[2][0];
     RHS2[3] = RHS[3][0];

        LUsolve(LHS2, RHS2);

     C_M_.ref()[celli].component(tensor::XX) = RHS2[0];
     C_M_.ref()[celli].component(tensor::XY) = RHS2[1];
     C_M_.ref()[celli].component(tensor::XZ) = RHS2[2];
     C_M_.ref()[celli].component(tensor::YX) = -RHS2[1];
     C_M_.ref()[celli].component(tensor::YY) = RHS2[0];
     C_M_.ref()[celli].component(tensor::YZ) = RHS2[3];
     C_M_.ref()[celli].component(tensor::ZX) = -RHS2[2];
     C_M_.ref()[celli].component(tensor::ZY) = -RHS2[3];
     C_M_.ref()[celli].component(tensor::ZZ) = RHS2[0];


     R.ref()[celli] = -twoSymm((C_M_.ref()[celli] & D[celli])) * sqrt(2.0) * mag(D[celli]);

     }
    

//    fvOptions.correct(R);
//    this->boundNormalStress(R);

    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
