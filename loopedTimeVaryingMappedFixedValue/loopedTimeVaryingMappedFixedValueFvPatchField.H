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

Class
    Foam::loopedTimeVaryingMappedFixedValueFvPatchField

Group
    grpInletBoundaryConditions grpCoupledBoundaryConditions

Description
    This boundary conditions interpolates the values from a set of supplied
    points in space and time.  Supplied data should be specified in
    constant/boundaryData/\<patchname\> where:
    - points : pointField with locations
    - ddd    : supplied values at time ddd
    The points should be more or less on a plane since they get triangulated
    in 2-D.

    At startup, this condition generates the triangulation and performs a
    linear interpolation (triangle it is in and weights to the 3 vertices)
    for every face centre.

    Values are interpolated linearly between times.

    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        setAverage   | flag to activate setting of average value | yes |
        perturb      | perturb points for regular geometries | no | 1e-5
        fieldTableName | alternative field name to sample | no| this field name
    \endtable

    /verbatim
    myPatch
    {
        type            timeVaryingMappedFixedValue;
        setAverage      false;
        //perturb       0.0;
        //fieldTableName samples;
    }
    /endverbatim

Note
    Switch on debug flag to have it dump the triangulation (in transformed
    space) and transform face centres.

SeeAlso
    Foam::fixedValueFvPatchField

SourceFiles
    loopedTimeVaryingMappedFixedValueFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef loopedTimeVaryingMappedFixedValueFvPatchField_H
#define loopedTimeVaryingMappedFixedValueFvPatchField_H

#include "fixedValueFvPatchFields.H"
#include "FixedList.H"
#include "instantList.H"
#include "pointToPointPlanarInterpolation.H"
#include "DataEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class loopedTimeVaryingMappedFixedValueFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class loopedTimeVaryingMappedFixedValueFvPatchField
:
    public fixedValueFvPatchField<Type>
{
    // Private data

        //- Name of the field data table, defaults to the name of the field
        word fieldTableName_;

        //- If true adjust the mapped field to maintain average value
        bool setAverage_;

        //- Fraction of perturbation (fraction of bounding box) to add
        scalar perturb_;

        //- Time limit, the ending time stant of sampled data afterwhich restarting from the first sampled time instant        
        scalar timeLimit_;

        //- 2D interpolation
        autoPtr<pointToPointPlanarInterpolation> mapperPtr_;

        //- List of boundaryData time directories
        instantList sampleTimes_;

        //- Current starting index in sampleTimes
        label startSampleTime_;

        //- Interpolated values from startSampleTime
        Field<Type> startSampledValues_;

        //- If setAverage: starting average value
        Type startAverage_;

        //- Current end index in sampleTimes
        label endSampleTime_;

        //- Interpolated values from endSampleTime
        Field<Type> endSampledValues_;

        //- If setAverage: end average value
        Type endAverage_;

        //- Time varying offset values to interpolated data
        autoPtr<DataEntry<Type> > offset_;




public:

    //- Runtime type information
    TypeName("loopedTimeVaryingMappedFixedValue");


    // Constructors

        //- Construct from patch and internal field
        loopedTimeVaryingMappedFixedValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        loopedTimeVaryingMappedFixedValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given loopedTimeVaryingMappedFixedValueFvPatchField
        //  onto a new patch
        loopedTimeVaryingMappedFixedValueFvPatchField
        (
            const loopedTimeVaryingMappedFixedValueFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        loopedTimeVaryingMappedFixedValueFvPatchField
        (
            const loopedTimeVaryingMappedFixedValueFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
            (
                new loopedTimeVaryingMappedFixedValueFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        loopedTimeVaryingMappedFixedValueFvPatchField
        (
            const loopedTimeVaryingMappedFixedValueFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type> > clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type> >
            (
                new loopedTimeVaryingMappedFixedValueFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return startSampledValues
            const Field<Type> startSampledValues()
            {
                 return startSampledValues_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<Type>&,
                const labelList&
            );


       

            //- Find boundary data inbetween current time and interpolate
            void checkTable();


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "loopedTimeVaryingMappedFixedValueFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
