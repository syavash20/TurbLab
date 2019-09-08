/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::channelIndex

Description
    Does averaging of fields over layers of cells. Assumes layered mesh.

SourceFiles
    meshIndex.C

\*---------------------------------------------------------------------------*/

#ifndef meshIndex_H
#define meshIndex_H

#include "regionSplit.H"
#include "direction.H"
#include "scalarField.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                           Class channelIndex Declaration
\*---------------------------------------------------------------------------*/

class meshIndex
{

    // Private data

        static const NamedEnum<vector::components, 3> vectorComponentsNames_;

        //- direction to sort
        const direction dir_;

        //- directions for output
        const direction dirAlt1_;
        const direction dirAlt2_;

        //- Per cell the global region
        autoPtr<regionSplit> cellRegion_;

        // number of cells
        scalar cellCount_;       

        //- Per global region the number of cells (scalarField so we can use
        //  field algebra)
        scalarField regionCount_;

        //- From sorted region back to unsorted global region
        labelList sortMap_;

        //- Sorted components of cell centres
        scalarField x_;
        scalarField y_;
        scalarField z_;



    // Private Member Functions

        void walkOppositeFaces
        (
            const polyMesh& mesh,
            const labelList& startFaces,
            boolList& blockedFace
        );

        void calcLayeredRegions
        (
            const polyMesh& mesh,
            const labelList& startFaces
        );

        //- Disallow default bitwise copy construct and assignment
        meshIndex(const meshIndex&);
        void operator=(const meshIndex&);


public:

    // Constructors

        //- Construct from dictionary
        meshIndex(const polyMesh&, const dictionary&);

        //- Construct from supplied starting faces
        meshIndex
        (
            const polyMesh& mesh,
            const labelList& startFaces,
            const direction dir,
            const direction dirAlt1,
            const direction dirAlt2
        );


    // Member Functions

        // Access

            //- Sum field per region
            template<class T>
            Field<T> regionSum(const Field<T>& cellField) const;

            //- collapse a field to a line
            template<class T>
            Field<T> collapse
            (
                const Field<T>& vsf
            ) const;

            //- return the field of X,Y and Z locations from the cell centres
            const scalarField& x() const
            {
                return x_;
            }
            const scalarField& y() const
            {
                return y_;
            }
            const scalarField& z() const
            {
                return z_;
            }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "meshIndexTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
