// Copyright held by original authors

#include "IncompressibleTurbulenceModel.H"
#include "incompressible/transportModel/transportModel.H"
#include "addToRunTimeSelectionTable.H"
#include "makeTurbulenceModel.H"

#include "LESModel.H"

namespace Foam
{
    typedef IncompressibleTurbulenceModel<transportModel>
        transportModelIncompressibleTurbulenceModel;
    typedef LESModel<transportModelIncompressibleTurbulenceModel>
        LEStransportModelIncompressibleTurbulenceModel;
}

#include "dynamicSmagorinsky.H"
makeTemplatedTurbulenceModel
(
    transportModelIncompressibleTurbulenceModel,
    LES,
    dynamicSmagorinsky
);
