/*---------------------------------------------------------------------------*\
sensoryPointCreate.cpp
    See header file for information on purpose

Copyright (c) 2004 R. Andrew Goodwin.
    See license file for information on usage
\*---------------------------------------------------------------------------*/

#include "sensoryPointCreate.h"

using std::cout;
using std::endl;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void SensoryOvoid::createSensoryPoints
(
    const int*  FSOlimit,
    double      fishSensoryFieldVars_acclMagAtFish,
    const int*  nCoeff,
    double      coefficients[*nCoeff],
    double      fishBodyLengths_fishNumber,
    int&        validSVorient_XYZ_NPm1,
    double&     SV_angleOff_CFD_XYZ_NPm1,
    double&     ovoidLength,
    double      fishSP_angleOff_SV[*FSOlimit-2],
    double      SP_relative_XYZ[][3]
)
{
    double      sensoryOvoidX[*FSOlimit],
                sensoryOvoidY[*FSOlimit],
                sensoryOvoidZ[*FSOlimit];

    double      SPminDistance,
                SP_XY_resultant,
                SPangleOffCFD,
                SPangleOffSV,
                SPangleOffSVradian,
                lengthFactor,
                widenessFactor,
                heigthFactor;

    // First sensory point is the fish itself and, therefore, there is no
    // displacement distance defined for this sensory point:
    sensoryOvoidX[0] = 0.0;
    sensoryOvoidY[0] = 0.0;
    sensoryOvoidZ[0] = 0.0;

/*  Other sensory points:
    positive (+) "x" is in the direction the head of the fish is pointing
    negative (-) "x" is in the direction of the fish's tail
    positive (+) "y" is to the left of the fish
    negative (-) "y" is to the right of the fish
    positive (+) "z" is in the direction above the fish (adverse to gravity)
    negative (-) "z" is in the direction below the fish (direction of gravity) */

    SPminDistance = fishBodyLengths_fishNumber; // Body length

    lengthFactor    = coefficients[0];
    widenessFactor  = coefficients[1];
    heigthFactor    = coefficients[2];

    // Sensory points are defined as a displacement distance of
    // (x,y,z) = ( SENOVOIDX(fso),SENOVOIDY(fso),SENOVOIDZ(fso) )
    // from the location of the fish

    sensoryOvoidX[1] =  SPminDistance * lengthFactor;
    sensoryOvoidY[1] =  0;
    sensoryOvoidZ[1] =  0;

    sensoryOvoidX[2] = -SPminDistance * lengthFactor;
    sensoryOvoidY[2] =  0;
    sensoryOvoidZ[2] =  0;

    sensoryOvoidX[3] =  0;
    sensoryOvoidY[3] =  SPminDistance * widenessFactor;
    sensoryOvoidZ[3] =  0;

    sensoryOvoidX[4] =  0;
    sensoryOvoidY[4] = -SPminDistance * widenessFactor;
    sensoryOvoidZ[4] =  0;

    sensoryOvoidX[5] =  0;
    sensoryOvoidY[5] =  0;
    sensoryOvoidZ[5] =  SPminDistance * heigthFactor;

    sensoryOvoidX[6] =  0;
    sensoryOvoidY[6] =  0;
    sensoryOvoidZ[6] = -SPminDistance * heigthFactor;


    // Set default values (fish facing in CFD x direction) ====================
    ovoidLength = std::sqrt
                  (
                        (SPminDistance * lengthFactor)
                      * (SPminDistance * lengthFactor)
                      + (SPminDistance * widenessFactor)
                      * (SPminDistance * widenessFactor)
                      + (SPminDistance * heigthFactor)
                      * (SPminDistance * heigthFactor)
                  );

    // Unused in BehaviorRule.f90 - requires transformation to fish coordinates
    fishSP_angleOff_SV[0] =   0.0;
    fishSP_angleOff_SV[1] =   0.0;
    fishSP_angleOff_SV[2] = 180.0;
    fishSP_angleOff_SV[3] =  90.0;
    fishSP_angleOff_SV[4] = -90.0;

    for (int FSO = 1; FSO < *FSOlimit; FSO++)
    {
        SP_relative_XYZ[FSO][0] = sensoryOvoidX[FSO];   // Distance (x-component in fish coordinates) from fish center to sensory point
        SP_relative_XYZ[FSO][1] = sensoryOvoidY[FSO];   // Distance (y-component in fish coordinates) from fish center to sensory point
        SP_relative_XYZ[FSO][2] = sensoryOvoidZ[FSO];   // Distance (z-component in fish coordinates) from fish center to sensory point
    }


    // Calculate Sensory Point Distances ======================================
    // (i.e., that override the default distances) for actual fish orientation
    // and/or if the SP were varied before. Introduces deviation in the
    // 7th decimal place, probably due to approximation of pi, ignored.
    if (validSVorient_XYZ_NPm1 == 1)   // Valid xy-plane orientation of fish axis can be determined
    {
        for (int FSO = 1; FSO < *FSOlimit-2; FSO++)   // Just for x/y-plane
        {
            SP_XY_resultant = std::sqrt
                              (
                                    sensoryOvoidX[FSO]*sensoryOvoidX[FSO]
                                   +sensoryOvoidY[FSO]*sensoryOvoidY[FSO]
                              );

            // atan2 computes the arc tangent of y/x using the signs of arguments to
            // determine the correct quadrant.
            SPangleOffCFD = std::atan2
                            (
                                sensoryOvoidY[FSO],
                                sensoryOvoidX[FSO]
                            )
                            *(180.0/3.14159);

            // Make ovoid angle relative to previous swim vector angle
            SPangleOffSV = SPangleOffCFD + SV_angleOff_CFD_XYZ_NPm1;

            // Perform bounds check
            if (SPangleOffSV < 0.0)
            {
                SPangleOffSV += 360.0;
            }
            else if (SPangleOffSV > 360.0)
            {
                SPangleOffSV -= 360.0;
            }

            // Ensure that angles go from 0.0 to 180.0 and 0.0 to -180.0
            if (SPangleOffSV > 180.0)
            {
                SPangleOffSV -= 360.0;
            }

            if      (SPangleOffSV >= 0.0 &&  SPangleOffSV <= 90.0)
            {
                SPangleOffSVradian =         SPangleOffSV * (3.14159/180.0);
                SP_relative_XYZ[FSO][0] = SP_XY_resultant * std::cos(SPangleOffSVradian);
                SP_relative_XYZ[FSO][1] = SP_XY_resultant * std::sin(SPangleOffSVradian);
            }
            else if (SPangleOffSV > 90.0 &&  SPangleOffSV <= 180.0)
            {
                SPangleOffSVradian =  (180.0-SPangleOffSV) * (3.14159/180.0);
                SP_relative_XYZ[FSO][0] = -SP_XY_resultant * std::cos(SPangleOffSVradian);
                SP_relative_XYZ[FSO][1] =  SP_XY_resultant * std::sin(SPangleOffSVradian);
            }
            else if (SPangleOffSV <= 0.0 &&  SPangleOffSV >= -90.0)
            {
                SPangleOffSVradian =        -SPangleOffSV  * (3.14159/180.0);
                SP_relative_XYZ[FSO][0] =  SP_XY_resultant * std::cos(SPangleOffSVradian);
                SP_relative_XYZ[FSO][1] = -SP_XY_resultant * std::sin(SPangleOffSVradian);
            }
            else if (SPangleOffSV < -90.0 && SPangleOffSV >= -180.0)
            {
                SPangleOffSVradian =  (180.0+SPangleOffSV) * (3.14159/180.0);
                SP_relative_XYZ[FSO][0] = -SP_XY_resultant * std::cos(SPangleOffSVradian);
                SP_relative_XYZ[FSO][1] = -SP_XY_resultant * std::sin(SPangleOffSVradian);
            }
            else
            {
                cout << "    SPcreate.cpp: Improper SPangleOffSV value."
                    << endl;
                return;
            }

        }   // End for FSO
    }   // End if validSVorient_XYZ_NPm1
    else
    {
        cout << "    SPcreate.cpp: Invalid xy orientation - falling back "
            "to default SP relative distance" << endl;
    }

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
