//     fileName path(UMean.rootPath()/UMean.caseName()/"span-averaged_field");
//     mkDir(path);

// VAR1



    volScalarField avgUx
    (
	IOobject
	(
	    "avgUx",
	    runTime.timeName(),
	    mesh,
	    IOobject::NO_READ,
	    IOobject::AUTO_WRITE
	),
        mesh,
        dimensionedScalar("avgUx", dimVelocity, 0.0)
    );
 
    volScalarField avgUy
    (
	IOobject
	(
	    "avgUy",
	    runTime.timeName(),
	    mesh,
	    IOobject::NO_READ,
	    IOobject::AUTO_WRITE
	),
        mesh,
        dimensionedScalar("avgUy", dimVelocity, 0.0)
    );

    volScalarField avgUz
    (
	IOobject
	(
	    "avgUz",
	    runTime.timeName(),
	    mesh,
	    IOobject::NO_READ,
	    IOobject::AUTO_WRITE
	),
        mesh,
        dimensionedScalar("avgUz", dimVelocity, 0.0)
    );

    volScalarField avgP
    (
	IOobject
	(
	    "avgP",
	    runTime.timeName(),
	    mesh,
	    IOobject::NO_READ,
	    IOobject::AUTO_WRITE
	),
        mesh,
        dimensionedScalar("avgP", dimPressure, 0.0)
    );

    volScalarField avgS
    (
        IOobject
        (
            "avgS",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("avgS", dimless, 0.0)
    );


    volScalarField avgK
    (
         IOobject
         (
            "avgK",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("avgK", dimVelocity*dimVelocity, 0.0)
    );	


    volScalarField avgNuSgs
    (
          IOobject
          (
            "avgNuSgs",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("avgNuSgs", dimVelocity*dimVelocity*dimTime, 0.0)
    );

    volScalarField pP2M
    (
	IOobject
	(
	    "pP2M",
	    runTime.timeName(),
	    mesh,
	    IOobject::NO_READ,
	    IOobject::AUTO_WRITE
	),
        mesh,
        dimensionedScalar("pP2M", dimPressure*dimPressure, 0.0)
    );

    volScalarField UPrime2MeanXX
    (
	IOobject
	(
	    "UPrime2MeanXX",
	    runTime.timeName(),
	    mesh,
	    IOobject::NO_READ,
	    IOobject::AUTO_WRITE
	),
        mesh,
        dimensionedScalar("UPrime2MeanXX", dimVelocity*dimVelocity, 0.0)
    );
    volScalarField UPrime2MeanYY
    (
	IOobject
	(
	    "UPrime2MeanYY",
	    runTime.timeName(),
	    mesh,
	    IOobject::NO_READ,
	    IOobject::AUTO_WRITE
	),
        mesh,
        dimensionedScalar("UPrime2MeanYY", dimVelocity*dimVelocity, 0.0)
    );
    
    volScalarField UPrime2MeanZZ
    (
	IOobject
	(
	    "UPrime2MeanZZ",
	    runTime.timeName(),
	    mesh,
	    IOobject::NO_READ,
	    IOobject::AUTO_WRITE
	),
        mesh,
        dimensionedScalar("UPrime2MeanZZ", dimVelocity*dimVelocity, 0.0)
    );
    volScalarField UPrime2MeanXY
    (
        IOobject
        (
            "UPrime2MeanXY",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("UPrime2MeanXY", dimVelocity*dimVelocity, 0.0)
    );
    volScalarField UPrime2MeanXZ
    (
        IOobject
        (
            "UPrime2MeanXZ",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("UPrime2MeanXZ", dimVelocity*dimVelocity, 0.0)
    );
    volScalarField UPrime2MeanYZ
    (
        IOobject
        (
            "UPrime2MeanYZ",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("UPrime2MeanYZ", dimVelocity*dimVelocity, 0.0)
    );
    
    volVectorField Uvec
    (
	IOobject
	(
	    "Uvec",
	    runTime.timeName(),
	    mesh,
	    IOobject::NO_READ,
	    IOobject::AUTO_WRITE
	),
        mesh,
        dimensionedVector
        (
                "Uvec",
                dimVelocity,
                vector::zero          
        )        
    );
    
    volSymmTensorField UP2M
    (
        IOobject
        (
            "UP2M",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor
        (
                "UP2M",
                dimVelocity*dimVelocity,
                symmTensor::zero
        )
    );

//VAR1
    if (UMeanHeader.headerOk())
    {
    volVectorField UMean(UMeanHeader, mesh);
    scalarField spanavgUMeanX
    (
        meshIndexing.collapse(UMean.component(vector::X)())
    );
    avgUx.internalField() = spanavgUMeanX;
    avgUx.boundaryField() = UMean.boundaryField().component(vector::X);
    
 //   avgUx.write();    
// VAR2
    scalarField spanavgUMeanY
    (
        meshIndexing.collapse(UMean.component(vector::Y)())
    );
    
    avgUy.internalField() = spanavgUMeanY;
    avgUy.boundaryField() = UMean.boundaryField().component(vector::Y);

//    avgUy.write();
// VAR3
    scalarField spanavgUMeanZ
    (
        meshIndexing.collapse(UMean.component(vector::Z)())
    );
    
    avgUz.internalField() = spanavgUMeanZ;
    avgUz.boundaryField() = UMean.boundaryField().component(vector::Z);
        Uvec.replace(0,avgUx); 
        Uvec.replace(1,avgUy);
        Uvec.replace(2,avgUz);
        Uvec.write();

    }
//    avgUz.write();
// VAR4

    if (pMeanHeader.headerOk())
    {
    volScalarField pMean(pMeanHeader, mesh);
    scalarField spanavgpMean
    (
        meshIndexing.collapse(pMean)
    );
    
    avgP.internalField() = spanavgpMean;
//    avgP.boundaryField() = pMean.boundaryField();
    avgP.write();

    }

    if (SMeanHeader.headerOk())
    {
     volScalarField SMean(SMeanHeader, mesh);
    scalarField spanavgSMean
    (
        meshIndexing.collapse(SMean)
    );

    avgS.internalField() = spanavgSMean;
//    avgP.boundaryField() = pMean.boundaryField();
    avgS.write();
    }

    if (kMeanHeader.headerOk())
    {
    volScalarField kMean(kMeanHeader, mesh);   
    scalarField spanavgkMean
    (
         meshIndexing.collapse(kMean)
    );
    avgK.internalField() = spanavgkMean;
//    avgP.boundaryField() = pMean.boundaryField();
    avgK.write();   
    }

    if (nuSgsMeanHeader.headerOk())
    {
    volScalarField nuSgsMean(nuSgsMeanHeader, mesh);
    scalarField spanavgnuSgsMean
    (
	 meshIndexing.collapse(nuSgsMean)
    );
    avgNuSgs.internalField() = spanavgnuSgsMean;
 //    avgP.boundaryField() = pMean.boundaryField();
    avgNuSgs.write();
    }

    if (pPrime2MeanHeader.headerOk())
    {
    volScalarField pPrime2Mean(pPrime2MeanHeader, mesh);	   
    scalarField spanavgpPrime2
    (
        meshIndexing.collapse(pPrime2Mean)
    );
    
    pP2M.internalField() = spanavgpPrime2;
//    pP2M.boundaryField() = pPrime2Mean.boundaryField();
    pP2M.write();
    }

    if (UPrime2MeanHeader.headerOk())
    {
    volSymmTensorField UPrime2Mean(UPrime2MeanHeader, mesh);
    scalarField spanavgUPrim2XX
    (
        meshIndexing.collapse(UPrime2Mean.component(symmTensor::XX)())
    );
    
    UPrime2MeanXX.internalField() = spanavgUPrim2XX;    
    UPrime2MeanXX.boundaryField() = UPrime2Mean.boundaryField().component(symmTensor::XX);
//    UPrime2MeanXX.write();

    
    scalarField spanavgUPrim2YY
    (
        meshIndexing.collapse(UPrime2Mean.component(symmTensor::YY)())
    );
    
    UPrime2MeanYY.internalField() = spanavgUPrim2YY;
    UPrime2MeanYY.boundaryField() = UPrime2Mean.boundaryField().component(symmTensor::YY);
//    UPrime2MeanYY.write();
    
    
    
    scalarField spanavgUPrim2ZZ
    (
        meshIndexing.collapse(UPrime2Mean.component(symmTensor::ZZ)())
    );
    
    UPrime2MeanZZ.internalField() = spanavgUPrim2ZZ;
    UPrime2MeanZZ.boundaryField() = UPrime2Mean.boundaryField().component(symmTensor::ZZ);
//    UPrime2MeanZZ.write();
    scalarField spanavgUPrim2XY
    (
        meshIndexing.collapse(UPrime2Mean.component(symmTensor::XY)())
    );

    UPrime2MeanXY.internalField() = spanavgUPrim2XY;
    UPrime2MeanXY.boundaryField() = UPrime2Mean.boundaryField().component(symmTensor::XY);

    scalarField spanavgUPrim2XZ
    (
        meshIndexing.collapse(UPrime2Mean.component(symmTensor::XZ)())
    );

    UPrime2MeanXZ.internalField() = spanavgUPrim2XZ;
    UPrime2MeanXZ.boundaryField() = UPrime2Mean.boundaryField().component(symmTensor::XZ);

    scalarField spanavgUPrim2YZ
    (
        meshIndexing.collapse(UPrime2Mean.component(symmTensor::YZ)())
    );

    UPrime2MeanYZ.internalField() = spanavgUPrim2YZ;
    UPrime2MeanYZ.boundaryField() = UPrime2Mean.boundaryField().component(symmTensor::YZ);


        UP2M.replace(0, UPrime2MeanXX);
        UP2M.replace(3, UPrime2MeanYY);
        UP2M.replace(5, UPrime2MeanZZ);
        UP2M.replace(1, UPrime2MeanXY);
        UP2M.replace(2, UPrime2MeanXZ);
        UP2M.replace(4, UPrime2MeanYZ);
        UP2M.write();     

        }    



 
