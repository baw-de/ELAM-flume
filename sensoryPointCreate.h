/*---------------------------------------------------------------------------*\
Application
    sensoryPointCreate

Purpose
    To create (x,y,z) coordinates for all sensory points of a given fish center
    with respect to the direction the fish is facing
    
Copyright (c) 2021 David Gisen, Bundesanstalt fuer Wasserbau
    See license file for information on usage

\*---------------------------------------------------------------------------*/

// Header guard
#ifndef sensoryPointCreate_H
#define sensoryPointCreate_H

#include <cmath>        // log10
#include <iostream>     // Input/Output, cout, cin

class SensoryOvoid
{

// Constructor, destructor and assignment operator are compiler-generated

public:

    void createSensoryPoints
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
        double      SP_XYZ_relative[][3]
    );

};

#endif   // sensoryPointCreate_H
