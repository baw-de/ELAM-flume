/*---------------------------------------------------------------------------*\
Application
    hydroInterpolation

Purpose
    Get CFD data from OpenFOAM at arbitrary Cartesian coordinates for fish
    location and all sensory points. Store in fishSensoryVelocity_NP and
    fishSensoryFieldVars_NP.
    Check for boundary collisions. Return patch type, number, and hitPoint.

Copyright (c) 2021 David Gisen, Bundesanstalt fuer Wasserbau
    See license file for information on usage

\*---------------------------------------------------------------------------*/

// Header guard
#ifndef hydroInterpolation_H
#define hydroInterpolation_H

#include "volFields.H"          // Diverse fields called by fvCFD.H
#include "fvCFD.H"              // OpenFOAM Finite Volume tools - location in Makefile
#include <string>               // Sequences of chars
#include "interpolationCellPoint.H"   // Interpolate within a cell using triangulation
#include <iostream>             // Input/Output, cout, cin
#include "treeDataFace.H"       // Search faces
#include "meshSearch.H"         // mesh search engine
//#include <ctime>                // Speed optimization

class HydroInterpolation
{

// Constructor, destructor and assignment operator are compiler-generated

public:

    void interpolateAtFSO
    (
        const fvMesh&       mesh,
        const meshSearch&   searchEngine1,
        const polyMesh::cellDecomposition& decompMode,
        int                 timeStep,
        const int*          nTecFieldVar, // still unsure why for const value * is required and not &
        const int*          FSOlimit,
        double              fishSensoryLocation_fishNumber[][3],
        double              fishSensoryVelocity_FSO[][3],
        double              fishSensoryFieldVars_FSO[][*nTecFieldVar],
        int                 SPfound_FSO[*FSOlimit],
        const interpolationCellPoint<scalar>& triangulateCellsAlpha,
        // const interpolationCellPoint<scalar>& triangulateCellsP,
        const interpolationCellPoint<scalar>& triangulateCellsAcclMag,
        // const interpolationCellPoint<scalar>& triangulateCellsGradUsum,
        const interpolationCellPoint<scalar>& triangulateCellsK,
        const interpolationCellPoint<vector>& triangulateCellsUMean,
        bool&               accelOn,
        bool&               TKEconstOn,
        bool&               extraDiagnostics
    );


    void getWallDistances
    (
        const fvMesh&       mesh,
        const meshSearch&   searchEngine1,
        double&             xFish,  // double& takes 'reference' = value itself.
        double&             yFish,  // double* takes 'pointer', gets value from address.
        double&             zFish,
        double&             visualRange,
        const int*          FSOlimit,
        double              SP_XYZ_relative[][3],       // 2D array
        double              wallDistances[*FSOlimit]    // 1D array, best way (clearest)
    );


    std::string hitBoundary_type
    (
        const fvMesh&       mesh,
        const meshSearch&   searchEngine1,
        double&             xStart,
        double&             yStart,
        double&             zStart,
        double&             xEnd,
        double&             yEnd,
        double&             zEnd,
        int&                patchNumber
    );


    bool hitInternalPatch
    (
        const fvMesh&       mesh,
        const meshSearch&   searchEngine1,
        Foam::point&        pStart,
        Foam::point&        pEnd,
        Foam::point&        pHit
    );


    Foam::point resetFSOoutside
    (
        const fvMesh&       mesh,
        const meshSearch&   searchEngine1,
        Foam::point&        insidePoint,
        Foam::point&        outsidePoint,
        Foam::label&        cellID
    );

};

#endif   // hydroInterpolation_H
