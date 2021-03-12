/*---------------------------------------------------------------------------*\
Application
    updateFishLocation

Purpose
    Update particle position according to flow vector (passive) or swim vector
    (active). Check and correct out-of-bounds and note normal exit or reset
    position, if necessary. Skip (and transfer position to next timestep) if
    particle has left the domain.
    
Copyright (c) 2021 David Gisen, Bundesanstalt fuer Wasserbau.
    See license file for information on usage

\*---------------------------------------------------------------------------*/

#include "volFields.H"              // For diverse fields called by fvCFD.H
#include "fvCFD.H"                  // OpenFOAM Finite Volume tools
#include "hydroInterpolation.h"     // Get CFD data
#include "cellSet.H"                // cellSet
#include <iostream>                 // Input/Output, cout, cin
#include <string>                   // Sequences of chars

using std::cout;
using std::endl;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void updateFishLocation
(
    const fvMesh&   mesh,
    const meshSearch& searchEngine1,
    const polyMesh::cellDecomposition decompMode,
    int&            fishNumber,
    bool&           passiveTransport,
    double&         dt,
    double&         xFish_NPp1,
    double&         yFish_NPp1,
    double&         zFish_NPp1,
    double&         xFish_NP,
    double&         yFish_NP,
    double&         zFish_NP,
    double&         u,
    double&         v,
    double&         w,
    double&         uFish,
    double&         vFish,
    double&         wFish,
    int&            fishAttribute_NP_fish,
    int&            fishAttribute_NPp1_fish,
    int             boundaryHitCount_Type[][2],
    bool&           extraDiagnostics,
    const           interpolationCellPoint<scalar>& triangulateCellsAlpha
)
{
    double              xFish, yFish, zFish;
    HydroInterpolation  OFinterp;
    Foam::label         fishCell(-1);
    Foam::point         fishPosition_NP(xFish_NP, yFish_NP, zFish_NP);
    Foam::point         pHit(-1,-1,-1);

    // Skip update if "particle" has left
    if (fishAttribute_NP_fish != 0)
    {
        xFish_NPp1 = xFish_NP;
        yFish_NPp1 = yFish_NP;
        zFish_NPp1 = zFish_NP;
        fishAttribute_NPp1_fish = fishAttribute_NP_fish;
        return;
    }

    // Update x,y,z positions:
    if (passiveTransport)
    {
        xFish_NPp1 = xFish_NP + u * dt;
        yFish_NPp1 = yFish_NP + v * dt;
        zFish_NPp1 = zFish_NP + w * dt;
    }
    else
    {
        xFish_NPp1 = xFish_NP + (uFish + u) * dt;
        yFish_NPp1 = yFish_NP + (vFish + v) * dt;
        zFish_NPp1 = zFish_NP + (wFish + w) * dt;
    }

    Foam::point fishPosition(xFish_NPp1, yFish_NPp1, zFish_NPp1);

    fishCell = searchEngine1.findCell(fishPosition,-1,true); // useTreeSearch

    // Check if individual should be removed from simulation
    // (because it was captured):
    if (fishCell == -1) // fish outside domain
    {
        double xHit, yHit, zHit;
        int patchNumber;
        std::string patchType;

        patchNumber = -1;
        patchType = "noHits";

        if (extraDiagnostics)
        {
            Info << "    New pos. outside: " << fishPosition << endl;
        }

        patchType = OFinterp.hitBoundary_type
        (
            mesh,
            searchEngine1,
            xFish_NP, yFish_NP, zFish_NP,
            xFish_NPp1, yFish_NPp1, zFish_NPp1,
            patchNumber
        );

        if ((patchType == "wall") || (patchType == "genericPatch"))
        // Warning: genericPatch catches "screen", but also wrongly named patches.
        {
            if (extraDiagnostics)
            {
                cout << "    Resetting fish position NP+1 @ " << patchType << endl;
            }

           // Reset position to old position - unconditionally stable
            xFish_NPp1 = xFish_NP;
            yFish_NPp1 = yFish_NP;
            zFish_NPp1 = zFish_NP;

            return;
        }
        else if (patchType == "patch")  // Fish leaves simulation
        {
            cout << "    updateFishLocation: Exit of fish #" << fishNumber+1
                << endl;

            fishAttribute_NPp1_fish = -1;
            boundaryHitCount_Type[patchNumber][0]++;

            return;
        }
        else
        {
            cout << "ABORT: updateFishLocation: Unknown boundary patch "
                "type: " << patchType << endl;
            exit(EXIT_FAILURE);
        }

    }   // END if (fishPosition outside mesh)


    // check if slot wall was crossed (fish still inside):
     // xxx magic number 1.0 depends on actual domain!
    if (fishPosition_NP.x() > -1.0 && fishPosition_NP.x() < 1.0)
    {
        if (OFinterp.hitInternalPatch
        (
            mesh,
            searchEngine1,
            fishPosition_NP,
            fishPosition,
            pHit
        ))
        {
            // Reset position to old position
              xFish_NPp1 = xFish_NP;
              yFish_NPp1 = yFish_NP;
              zFish_NPp1 = zFish_NP;
              Info << "    updateFishLocation: Hit generic patch "
                    "inside domain @ " << pHit << ". Reset." << endl;
              return;
        }
    }

    // check + move down if in air:
    // xxx magic number "0.1 m" depends on actual domain!
    while (triangulateCellsAlpha.interpolate(fishPosition, fishCell) < 0.5)
    {
        if (extraDiagnostics)
        {
            cout << "    New fish center in air "
            "- move by -10 cm in z" << endl;
        }

        fishPosition.z() -= 0.1; // 0.1 m fixed, xxx magic number
        fishCell = mesh.findCell(fishPosition, decompMode);

        if (fishCell == -1)
        {
            // Reset position to old position - safest, but trapping danger
            xFish_NPp1 = xFish_NP;
            yFish_NPp1 = yFish_NP;
            zFish_NPp1 = zFish_NP;

            Info << "    updateFishLocation: Reset: alpha request "
                "outside domain @ " << fishPosition << endl;

            break;
        }

        zFish_NPp1 = fishPosition.z();
    }

    return;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
