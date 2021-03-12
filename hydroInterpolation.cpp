/*---------------------------------------------------------------------------*\
hydroInterpolation.cpp
    See header file for information on purpose

Copyright (c) 2021 David Gisen, Bundesanstalt fuer Wasserbau
    See license file for information on usage
\*---------------------------------------------------------------------------*/

#include "hydroInterpolation.h"

using namespace Foam;
using std::cout;
using std::endl;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void HydroInterpolation::interpolateAtFSO
(
    const       fvMesh& mesh,
    const       meshSearch& searchEngine1,
    const       polyMesh::cellDecomposition& decompMode,
    int         timeStep,
    const int*  nTecFieldVar,
    const int*  FSOlimit,
    double      fishSensoryLocation_FSO[][3], //_fishNumber statt FSO
    double      fishSensoryVelocity_FSO[][3],
    double      fishSensoryFieldVars_FSO[][*nTecFieldVar],
    int         SPfound_FSO[*FSOlimit],
    const       interpolationCellPoint<scalar>& triangulateCellsAlpha,
    // const       interpolationCellPoint<scalar>& triangulateCellsP,
    const       interpolationCellPoint<scalar>& triangulateCellsAcclMag,
    // const       interpolationCellPoint<scalar>& triangulateCellsGradUsum,
    const       interpolationCellPoint<scalar>& triangulateCellsK,
    const       interpolationCellPoint<vector>& triangulateCellsUMean,
    bool&       accelOn,
    bool&       TKEconstOn,
    bool&       extraDiagnostics
)
{
    double      xFSO, yFSO, zFSO,
                alphaAtSP, pAtSP, acclMagAtSP,
                TKEatSP, THSatSP;

    Foam::point     FSO_position;
    Foam::vector    UatSP, in_out;
    Foam::label     FSO_cell(-1);

    xFSO = yFSO = zFSO = 0.0;

    Foam::point FSO_center = Foam::point(
                                fishSensoryLocation_FSO[0][0],
                                fishSensoryLocation_FSO[0][1],
                                fishSensoryLocation_FSO[0][2]);

    // Begin ovoid loop
    // FSO 1..7 in Fortran.
    // Executed just once (for center) if TS = 0, see break below
    for (int FSO = 0; FSO < *FSOlimit; FSO++)
    {
        xFSO = fishSensoryLocation_FSO[FSO][0];
        yFSO = fishSensoryLocation_FSO[FSO][1];
        zFSO = fishSensoryLocation_FSO[FSO][2];

        FSO_position = Foam::point(xFSO, yFSO, zFSO);

        // Will change to 0 if not in mesh
        SPfound_FSO[FSO] = 1;

        // Get cell ID, -1 if outside
        FSO_cell = searchEngine1.findCell(FSO_position,-1,true); // useTreeSearch

        // check could be omitted for FSO=0 (center), but FSO_cell is required anyway
        if (FSO_cell == -1)
        {
            if (extraDiagnostics)
            {
                cout << "    Resetting FSO position #" << FSO+1 << endl;
            }

            // get new FSO_position inside mesh
            FSO_position = this->resetFSOoutside
            (
                mesh,
                searchEngine1,
                FSO_center,     // Fish center / must be an inside point
                FSO_position,   // Current (old) Sensory Ovoid Point Position
                FSO_cell        // new FSO_cell
            );

            // Mark as "Not in mesh before reset"
            SPfound_FSO[FSO] = 0;

        } // end "if FSO not in mesh"

        // (Re-)Initialize variables
        UatSP       = Foam::vector::zero;
        alphaAtSP   =  0.0;
        pAtSP       = -1.0; // unused
        acclMagAtSP =  0.0;
        TKEatSP     =  0.0;
        THSatSP     = -1.0; // unused

        // Check for ovoid points in air
        in_out = Foam::vector::zero;
        in_out = FSO_position - FSO_center;

        alphaAtSP = triangulateCellsAlpha.interpolate(FSO_position, FSO_cell);
        while (alphaAtSP < 0.5)
        {
            if (extraDiagnostics)
            {
                cout << "    hydroInterp: FSO #" << FSO+1 << " in air "
                "- move by 25%" << endl;
            }

            // Move FSO_position inwards to get alpha > 0.5
            FSO_position -= 0.25 * in_out;
            FSO_cell = mesh.findCell(FSO_position, decompMode);

            if (extraDiagnostics)
            {
                Foam::Info << "    New FSO: " << FSO_position << Foam::endl;
            }

            if (FSO_cell == -1) continue; // Prevent crash at enclaves

            alphaAtSP =
                     triangulateCellsAlpha.interpolate(FSO_position, FSO_cell);
        }

        // Interpolate
        UatSP = triangulateCellsUMean.interpolate(FSO_position, FSO_cell);
        if (accelOn)
        {
            acclMagAtSP =
                triangulateCellsAcclMag.interpolate(FSO_position, FSO_cell);
        }
        if (TKEconstOn)
        {
            TKEatSP = triangulateCellsK.interpolate(FSO_position, FSO_cell);
        }
        // pAtSP = triangulateCellsP.interpolate(FSO_position, FSO_cell);
        //THSatSP =
        //    triangulateCellsGradUsum.interpolate(FSO_position, FSO_cell);

        // Set u, v, w, p, AcclM, TKE, THS
        fishSensoryVelocity_FSO[FSO][0]  = UatSP[0];    // u at sensory point
        fishSensoryVelocity_FSO[FSO][1]  = UatSP[1];    // v at sensory point
        fishSensoryVelocity_FSO[FSO][2]  = UatSP[2];    // w at sensory point

        fishSensoryFieldVars_FSO[FSO][0] = pAtSP;       // Pressure at sensory point
        fishSensoryFieldVars_FSO[FSO][1] = TKEatSP;     // Turbulence kinetic energy at sensory point
        fishSensoryFieldVars_FSO[FSO][2] = acclMagAtSP; // Acceleration magnitude at sensory point
        fishSensoryFieldVars_FSO[FSO][3] = THSatSP;     // Total hydraulic strain mag at sensory point

        // Omit all SP except center for TS = 0 - as in StrELAM_main.f:1945
        // reason probably efficiency, not needed for TS = 0 - DG.
        if (timeStep <= 0)
        {
            break;
        }
    }   // end for(FSO < *FSOlimit)

    return;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void HydroInterpolation::getWallDistances
// Compute distances from fish to wall
// SP2)     front (here i_SP==1, array indexing)
// SP4+5)   left and right
(
    const fvMesh& mesh,
    const meshSearch& searchEngine1,
    double& xFish,
    double& yFish,
    double& zFish,
    double& visualRange,
    const int* FSOlimit,
    double SP_XYZ_relative[][3],
    double wallDistances[*FSOlimit]
)
{
    Foam::point fishPoint(xFish, yFish, zFish);
    const Foam::fvPatchList& patches = mesh.boundary();
    int patchNumber;
    int numberOfHitpoints = 0;

    // initialize to infinity range - there is no wall detectable
    for ( int i_SP = 0; i_SP < *FSOlimit; i_SP++ )
    {
        wallDistances[i_SP] = VGREAT;
    }

    // check all sensory point directions for hits
    for ( int i_SP = 0; i_SP < *FSOlimit; i_SP++ )
    {
        // skip for all but front, left & right
        if (!((i_SP == 1) || (i_SP == 3) || (i_SP == 4))) continue;

        Foam::vector in_out(SP_XYZ_relative[i_SP][0],
                            SP_XYZ_relative[i_SP][1],
                            SP_XYZ_relative[i_SP][2]);

        Foam::point maxVisionP = fishPoint + in_out/mag(in_out) * visualRange;

        // Search line direction in -> out to catch internal faces
        pointIndexHit pHit = searchEngine1.intersection(fishPoint, maxVisionP);

        if (pHit.hit())
        {
            // get coordinates and compute distance
            Foam::point hitPoint = pHit.hitPoint();
            wallDistances[i_SP]  = mag(hitPoint-fishPoint);
        }
        // catch OF 4.1 bug#1: no hit was detected, but point is outside
        // catch OF 4.1 bug#2: point is inside, but not detected as such
        else if (searchEngine1.findCell(maxVisionP,-1,true) == -1)
        {
            // Simple solution. Theory: bug#1 happens in wall area where distance ~ visualRange
            //   bug#2 is rare and can be ignored
            wallDistances[i_SP] = visualRange;
            Info << "    getWallDistances: OF bug bypassed @ " << maxVisionP << endl;
        }
        else
        {
            // do nothing, wall distance was initialized & is kept at VGREAT
        }
    } // end: for i_SP

    return;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


std::string HydroInterpolation::hitBoundary_type
// "Ray shooting" from StrELAM_main.f
// Returns just boundary patch type
// Assumes outsideP to be outside the mesh
(
    const   fvMesh& mesh,
    const   meshSearch& searchEngine1,
    double& xStart,
    double& yStart,
    double& zStart,
    double& xEnd,
    double& yEnd,
    double& zEnd,
    int&    patchID
)
{
    Foam::point insideP(xStart, yStart, zStart);
    Foam::point outsideP(xEnd, yEnd, zEnd);
    std::string patchType = "noHits";

    pointIndexHit pHit = searchEngine1.intersection(insideP, outsideP);

    if (!pHit.hit())
    {
        // not hit despite updateFishLocation found it outside - bug in OF4.1
        // mimic behavior of ELAM-diss (reset on faulty mesh spots):
        Info << "    hitBoundary_type: OF bug bypassed." << endl;
        patchType = "wall";
        return patchType;
    }

    int faceID              = pHit.index();
    patchID                 = mesh.boundaryMesh().whichPatch(faceID);
    // constant reference to pp1 - reason: no copy to return = more efficient?
    const polyPatch& pp1    = mesh.boundaryMesh()[patchID];
    patchType               = pp1.type();

    return patchType;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


bool HydroInterpolation::hitInternalPatch
// Checks if fish trajectory collided with internal faces, such as the slot
(
    const fvMesh&       mesh,
    const meshSearch&   searchEngine1,
    Foam::point&        pStart,
    Foam::point&        pEnd,
    Foam::point&        pHit
)
{
    pointIndexHit pIndexHit = searchEngine1.intersection(pStart, pEnd);

    int faceID              = pIndexHit.index();
    int patchID             = mesh.boundaryMesh().whichPatch(faceID);
    // constant reference to pp1 - reason: no copy to return = more efficient
    const polyPatch& pp1    = mesh.boundaryMesh()[patchID];

    if ((pIndexHit.hit()) && pp1.inGroup("internal"))
    {
        pHit = pIndexHit.hitPoint();
        return true;
    }
    else // no hits found in "internal" patches
        return false;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


Foam::point HydroInterpolation::resetFSOoutside
(
    const fvMesh&       mesh,
    const meshSearch&   searchEngine1,
    Foam::point&        insideP,
    Foam::point&        outsideP,
    Foam::label&        cellID
)
{
    Foam::point newOutsideP = outsideP;
    Foam::vector in_out;
    int counter;
    int nResets;

    in_out = Foam::vector::zero;
    in_out = outsideP - insideP;
    cellID = -1;
    counter = 0;
    nResets = 3; // 2 not so much faster

    if (Foam::mag(in_out) == 0)
    {
        cout << "ABORT: resetFSOoutside: Zero movement, error happened before" << endl;
        exit(EXIT_FAILURE);
    }

//    Foam::Info << insideP << " " << outsideP << Foam::endl;

    do
    {
        counter++;

        if (counter > nResets)
        {
            cout << "ABORT: resetFSOoutside: Can't find point inside" << endl;
            Info << newOutsideP << Foam::endl;
            exit(EXIT_FAILURE);
        }

        newOutsideP -= 1.0/nResets*in_out;

        // Bugfix because 3*1/3 was not guaranteed to return the same insidePoint
        if (counter == nResets)
        {
            newOutsideP = insideP;
        }
        cellID = searchEngine1.findCell(newOutsideP);
    }
    while (cellID == -1);

    return newOutsideP;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
