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

\*---------------------------------------------------------------------------*/

#include "meshIndex.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
Foam::Field<T> Foam::meshIndex::regionSum(const Field<T>& cellField) const
{
    Field<T> regionField(cellRegion_().nRegions(), pTraits<T>::zero);

    forAll(cellRegion_(), cellI)
    {
        regionField[cellRegion_()[cellI]] += cellField[cellI];
    }

    // Global sum
    Pstream::listCombineGather(regionField, plusEqOp<T>());
    Pstream::listCombineScatter(regionField);

    return regionField;
}


template<class T>
Foam::Field<T> Foam::meshIndex::collapse
(
    const Field<T>& cellField
) const
{
    // Average and order
    const Field<T> summedField(regionSum(cellField));
    const Field<T> regionCountField = regionSum(scalarField(cellCount_, 1.0));

    Field<T> spanaveragedField(cellCount_, pTraits<T>::zero);

    forAll(cellRegion_(), cellI)
    {
        spanaveragedField[cellI] = summedField[cellRegion_()[cellI]]
                                 / regionCountField[cellRegion_()[cellI]];
    }

    return spanaveragedField;
}


// ************************************************************************* //
