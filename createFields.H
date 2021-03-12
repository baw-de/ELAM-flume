/*---------------------------------------------------------------------------*\
Application
    createFields

Purpose
    Read OpenFOAM variable fields
    
Copyright (c) 2021 David Gisen, Bundesanstalt fuer Wasserbau
    See license file for information on usage

\*---------------------------------------------------------------------------*/

    cout << "Reading field alpha.water" << endl;
    volScalarField alpha
    (
        IOobject
        (
            "alpha.water",
            lastTime.name(),   // runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


//    cout << "Reading field p" << endl;
//    volScalarField p
//    (
//        IOobject
//        (
//            "p",
//            lastTime.name(),   // runTime.timeName(),
//            mesh,
//            IOobject::MUST_READ,
//            IOobject::AUTO_WRITE
//        ),
//        mesh
//    );


    cout << "Reading field UMean" << endl;
    volVectorField UMean
    (
        IOobject
        (
            "UMean",
            lastTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    cout << "Reading field acclMag" << endl;
    volScalarField acclMag
    (
        IOobject
        (
            "acclMag",
            lastTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    cout << "Reading field k" << endl;
    // TKE
    volScalarField k
    (
        IOobject
        (
            "k",
            lastTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    cout << endl;
