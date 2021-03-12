/*---------------------------------------------------------------------------*\
Application
    ELAM main class

Purpose
    To control time and fish loops of BehaviorRule.f90 and to control
    OpenFOAM 4.1 data acquisition on an unstructured mesh

Copyright (c) 2021 David Gisen, Bundesanstalt fuer Wasserbau
    See license file for information on usage

\*---------------------------------------------------------------------------*/


// === Includes ===============================================================
// Quotes = locally, Angle brackets = System directory
#include "volFields.H"              // Diverse fields called by fvCFD.H
#include "fvCFD.H"                  // OpenFOAM Finite Volume tools
#include <fstream>                  // read files
#include "hydroInterpolation.h"     // get CFD data
#include "sensoryPointCreate.h"     // Create Sensory Points (SP)
#include <iostream>                 // Input/Output, cout, cin
#include <iomanip>                  // Format I/O
#include <sstream>                  // str2double
#include <string>                   // Sequences of chars
//#include <ctime>                    // Speed optimization


using std::cout;
using std::endl;


// === Function definitions ===================================================


void readDoubleLinedInputFile
(
    std::string* fileName,
    const int* arrayLength,
    double (*array_)
)
/* Skips every second line as there are comments, puts read double values into
   array until file end */
{
    std::ifstream inputFile;
    std::string line;
    int position;
    double inputValue;

    inputFile.open(*fileName, std::ios::in);
    if (inputFile.is_open())
    {
        position = 0;
        while ( getline(inputFile, line) )
        {
            if (position > *arrayLength-1)
            {
                cout << "ERROR in readDoubleLinedInputFile(): Input file "
                        "too large." << endl;
                exit(EXIT_FAILURE);
            }
            getline(inputFile, line);
            std::istringstream i(line);
            i >> inputValue;
            *(array_+position) = inputValue;
            position++;
        }
        inputFile.close();
    }
    else
    {
        cout << "ERROR in readDoubleLinedInputFile(): Unable to open file."
             << endl;
        exit(EXIT_FAILURE);
    }

    return;
}


bool checkInputValues
(
    const int* valuesLength,
    double values_[*valuesLength]
)
{
    if ((values_[0] != 0) && (values_[0] != 1))
    {
        cout << "Input 1 flawed." << endl;
        return false;
    }
    if (values_[1] <= 0.0)
    {
        cout << "Input 2 flawed." << endl;
        return false;
    }
    if (values_[2] <= 0.0)
    {
        cout << "Input 3 flawed." << endl;
        return false;
    }
    if (values_[3] < 1.0)
    {
        cout << "Input 4 flawed." << endl;
        return false;
    }
    if (values_[4] < 0.0 || values_[4] > 2.0)
    {
        cout << "Input 5 flawed." << endl;
        return false;
    }
    // values_[5] not needed
    if (values_[6] < 0 || values_[6] > 1 || (values_[6] != (int)values_[6]))
    {
        cout << "Input 7 flawed." << endl;
        return false;
    }
    if (values_[7] != -1.0)
    {
        cout << "Input 8 flawed." << endl;
        return false;
    }
    if ((values_[8] < 0.0) || (values_[8] > 3.0))
    {
        cout << "Input 9 flawed." << endl;
        return false;
    }
    if ((values_[9] <= 0.0) || (values_[9] >= 1000.0))
    {
        cout << "Input 10 flawed." << endl;
        return false;
    }
    if ((values_[10] < 1) || (values_[10] != (int)values_[10]))
    {
        cout << "Input 11 flawed." << endl;
        return false;
    }
    if (values_[11] != (bool)values_[11])
    {
        cout << "Input 12 flawed." << endl;
        return false;
    }
    if (values_[12] != (bool)values_[12])
    {
        cout << "Input 13 flawed." << endl;
        return false;
    }
    if (values_[13] != (bool)values_[13])
    {
        cout << "Input 14 flawed." << endl;
        return false;
    }
    if (values_[14] != (bool)values_[14])
    {
        cout << "Input 15 flawed." << endl;
        return false;
    }
    return true;
}


void readFishPositions
(
    std::string& fileName,
    double fishLocations[][3],
    int    nFish,       // pass by value, recommended for unmodified int
    double fishBodyLengths[nFish],
    const int* nCoeff,  // only required to determine array length
    double coefficients[*nCoeff],
    double (*k_M_FN),   // motivation coefficients
    double (*k_F_FN)    // fatigue coefficients
)
{
    int fish_i;
    std::ifstream inputFile(fileName); // ifstream: read-only; fstream: both
    double x, y, z;
    double bodyLength;
    double k_M_i, k_F_i;
    std::string condition;

    fish_i = 0;

    if (inputFile.is_open())
    {
        while(inputFile >> x >> y >> z >> bodyLength >> condition)
        {
            fishLocations[fish_i][0] = x;
            fishLocations[fish_i][1] = y;
            fishLocations[fish_i][2] = z;

            // memory corruption if fish_i > nFish - catch below
            fishBodyLengths[fish_i] = bodyLength;

            if (condition == "fast")
            {
                *(k_M_FN+fish_i) = coefficients[57];
                *(k_F_FN+fish_i) = coefficients[59];
            }
            else if (condition == "slow")
            {
                *(k_M_FN+fish_i) = coefficients[58];
                *(k_F_FN+fish_i) = coefficients[60];
            }
            else
            {
                cout << "ERROR in readFishPositions(): fish must be 'fast' or "
                        " 'slow'."
                     << endl;
                exit(EXIT_FAILURE);
            }

            fish_i++;
        }
    }
    else
    {
        cout << "ERROR in readFishPositions(): Unable to open file." << endl;
        exit(EXIT_FAILURE);
    }

    if (fish_i != nFish)
    {
        cout << "ERROR: Released fish count incorrect!" << endl;
        exit(EXIT_FAILURE);
    }

    return;
}


template <typename T1> void initializeArray3D
(
    int idim,
    int jdim,
    int kdim,
    T1 (*arrayToOutput)
)
{
    for (int i = 0; i < idim; i++)
    {
        for (int j = 0; j < jdim; j++)
        {
            for (int k = 0; k < kdim; k++)
            {
                (*(arrayToOutput+i*jdim*kdim+j*kdim+k)) = 0;
            }
        }
    }
    return;
}


int sumArray2Dconditionally
(
    int iDim,
    int (*array2D)
)
// Sums up the exits through patches
{
    int sum = 0;

    for (int i = 0; i < iDim; i++)
    {
        if ( (*(array2D+i*2+1)) == 1 )
        {
            sum += (*(array2D+i*2+0));
        }
    }
    return sum;
}


// Fortran wrapper
extern "C"
{
    // Forward declare function prototypes
    void behaviorrule_   // lowercase and underscore forced by Fortran 95
    (
        int*        fishNumber,
        int*        nFish,
        int*        timeStep,
        double*     dt,
        const int*  nAgents,
        const int*  nTecFieldVar,
        const int*  FSOlimit,
        int         (*SPfound_NP_FN_FSO),   // Pointer at last entry array
        double      (*fishSensoryLocation_FN_FSO)[3],
//        double      (*fishLocations_NPm1_FN)[3],
        double      (*fishSensoryVelocity_FN_FSO)[3],
        double      (*fishSensoryFieldVars_FN_FSO)[*nTecFieldVar],
        int*        seed,
        double*     fishSpeedResultant_NPm1_FN,
        double*     fishSpeedResultant_NP_FN,
        double      (*SV_velocityOff_CFD_XYZ_NPm1_FN_3),   // points to array[3]
        double      (*SV_velocityOff_CFD_XYZ_FN_3),
        double      (*SV_angleOff_CFD_XYZ_NPm1_FN)[2],
        double      (*SV_angleOff_CFD_XYZ_NP_FN_2),
        double*     SV_ao_SV_XY_NPm1_FN,
        double*     SV_ao_SV_XY_NP_FN,
        double      (*agtProb_NPm1_FN_Agent),
        double      (*agtProb_NP_FN_Agent),
        double      (*agtDetctMetrcAmb_NP_FN_Agent),
        double      (*agtDetctMetrcAmb_NPp1_FN_Agent),
        int         (*agtDecision_nTimeSteps)[*nFish][*nAgents],
        int         (*numDecisions_FN),
        int         (*validSVorient_XYZ_NP_FN_2),
        bool*       extraDiagnostics,
        bool*       writeToDebugFile,
        const int*  nCoeff,
        double      (*coefficients),
        double      (*SVvo_FV_XYZ_NP_FN_3),
        double      (*SVvo_SV_XYZ_NP_FN_3),
        double*     SVao_FV_XY_NP_FN,
        double      (*agtDetectMetric_NP_FN_Agent),
        double      (*agtUtility_NP_FN_Agent),
        double      (*fishSP_angleOff_SV_FSO),
        int         (*fishAttribute_NPm1_FN)[1],
        int*        nTimeSteps,
        double      (*sAvg_NPm1_FN_3),
        double      (*t_same_spot),
        double*     ovoidLength,
        double*     agtCostAvg_FN,
        double      (*k_M_FN),
        double      (*k_F_FN),
        double      (*wallDistances),
        bool*       accelOn,
        bool*       visualOn,
        bool*       lowVelOn,
        bool*       TKEconstOn,
        double*     fishBodyLengths_FN
    );

    void openoutputfiles_
    (
        bool*       writeToDebugFile,
        int&        longResults
    );

    void writesensorylocations_
    (
        int*        fishNumber,
        int*        timeStep,
        const int*  FSOlimit,
        double*     ovoidLength,
        double      (*SP_XYZ_relative_FSOlimit)[3]
    );

    void outputfishdata_zonesaretime_
    (
        int*        fishNumber,
        int*        nFish,
        int*        timeStep,
        int*        nTimeSteps,
        double*     dt,
        const int*  FSOlimit,
        const int*  nAgents,
        double      (*fishLocations_nTimeSteps)[*nFish][3],
        double      (*SV_velocityOff_CFD_XYZ_nTimeSteps)[*nFish][3],
        double      (*agtUtility_nTimeSteps)[*nFish][*nAgents],
        double      (*fishSpeedResultant_nTimeSteps)[*nFish],
        int         (*fishAttribute_nTimeSteps)[*nFish][1],
        int*        totalExitTally,
        double      (*fishSensoryVelocity_nFish)[*FSOlimit][3]
    );

    void outputfishpassageanddecisions_
    (
        int*        nFish,
        int*        timeStep,
        double*     dt,
        int*        numExitRoutes,
        int*        numBoundaries,
        int         (*boundaryHitCount_Type)[2]
    );

    void endofrunoutput_
    (
        int*        nFish,
        int*        timeStep,
        int*        nTimeSteps,
        double*     dt,
        int         (*numDecisions_nFish),
        int*        outIntV,
        double      (*fishLocations_nTimeSteps)[*nFish][3],
        double      (*SV_velocityOff_CFD_XYZ_nTimeSteps)[*nFish][3],
        double      (*fishSpeedResultant_nTimeSteps)[*nFish],
        int         (*fishAttribute_nTimeSteps)[*nFish][1],
        double      (*SV_angleOff_CFD_XYZ_nTimeSteps)[*nFish][2],
        int&        longResults,
        double      (*statesStatistics_nFish)[8],
        int&        numDecisionsSum
    );
}


// Forward declare function prototype
void updateFishLocation
(
    const fvMesh& mesh,
    const meshSearch& searchEngine1,
    const polyMesh::cellDecomposition decompMode,
    int&    fishNumber,
    bool&   passiveTransport,
    double& dt,
    double& xFish_NP,
    double& yFish_NP,
    double& zFish_NP,
    double& xFish_NPm1,
    double& yFish_NPm1,
    double& zFish_NPm1,
    double& u,
    double& v,
    double& w,
    double& uFish,
    double& vFish,
    double& wFish,
    int&    fishAttribute_NPm1_fish,
    int&    fishAttribute_NP_fish,
    int     boundaryHitCount_Type[][2],
    bool&   extraDiagnostics,
    const interpolationCellPoint<scalar>& triangulateCellsAlpha
);


// Main program: read argument count and arg. vector
int main(int argc, char* argv[])
{
// === Fixed size variable declarations =======================================
    // Scalar variables declaration
    int     fishNumber,
            nFish,
            timeStep,
            hydroValuesLength,
            FSO,            // fish sensory ovoid, sensory point
            nTimeSteps,     // begin console input vars
            outIntV,
            outputVfiles,
            //diskOrRAM,    // not needed
            longResults,
            vFishInpNum,
            CDNvFishInp,
            seed,           // end console input vars
            NiterSinceOut,  // Number of iterations since last output "ITBOUT"
            numExitRoutes,
            totalExitTally;

    // Initialization
            totalExitTally = 0;

    // Provides const int FSOlimit, nFhAttrb, nCoeff, nTecFieldVar, nAgents,
    //                    consoleInputLength
    #include "input/rules.inp"

    double  dt,
            xFish,
            yFish,
            zFish,
            ovoidLength;

    // Initialization
            ovoidLength = 0.0;

    bool    passiveTransport = false,
            accelOn          = false,
            visualOn         = false,
            lowVelOn         = false,
            TKEconstOn       = false,
            extraDiagnostics = false,
            writeToDebugFile = false;

//    // For duration analyses:
//     struct timespec tStart, tEnd;
//     double elapsed;

    // Fixed size array declarations
    double  consoleInput[consoleInputLength],
            coefficients[nCoeff] = {0};

    // Further type variables declaration + initialization
    Foam::argList args(argc, argv, true, true, false);   // false = prevent initialization
    Foam::label fishCell(-1);

// === Read basic input =======================================================

    // "Console" input - from file
    std::string fileName = "input/simSettings.inp";
    readDoubleLinedInputFile(&fileName, &consoleInputLength, consoleInput);
    if (!checkInputValues(&consoleInputLength, consoleInput))
    {
        cout << "ERROR in main(): Input values incorrect." << endl;
        exit(EXIT_FAILURE);
    }

    // Set checked values from read
    passiveTransport = !consoleInput[0];// console: volitional transport
    dt              = consoleInput[1];
    nTimeSteps      = consoleInput[2];
    outIntV         = consoleInput[3];  // outputInterval for V files
    outputVfiles    = consoleInput[4];
    // diskOrRam not needed
    longResults     = consoleInput[6];
    vFishInpNum     = consoleInput[7];  // changed
    CDNvFishInp     = consoleInput[8];  // unused
    seed            = consoleInput[9];
    nFish           = consoleInput[10];
    accelOn         = consoleInput[11];
    visualOn        = consoleInput[12];
    lowVelOn        = consoleInput[13];
    TKEconstOn      = consoleInput[14];

    fileName = "input/agentBehaviorCoefficients.inp";
    readDoubleLinedInputFile(&fileName, &nCoeff, coefficients);

    NiterSinceOut  = 0; // Start output w/ first timestep

// === Variable size array declarations =======================================

    // If not directly as here: always declare array size using variables -
    // good programming practice.

    int     // 4th dimension omitted in SPfound compared to StrELAM - zone# not needed
            SPfound[nTimeSteps][nFish][FSOlimit],
            agtDecision[nTimeSteps][nFish][nAgents],
            validSVorient_XYZ[nTimeSteps][nFish][2],
            fishAttribute[nTimeSteps+1][nFish][1], // 0=Fish ; 1=Invertebrate ; -1=To be removed from simulation
            numDecisions[nFish]; // Limit output to the point the fish left the domain

    double  fishSensoryLocations[nFish][FSOlimit][3],
            fishLocations[nTimeSteps+1][nFish][3],
            fishSensoryVelocity[nFish][FSOlimit][3],
            fishSensoryFieldVars[nFish][FSOlimit][nTecFieldVar],
            fishSpeedResultant[nTimeSteps][nFish],  // TS0 ini in BR
            agtProb[nTimeSteps][nFish][nAgents],
            agtDetctMetrcAmb[nTimeSteps+1][nFish][nAgents],
            SV_velocityOff_CFD_XYZ[nTimeSteps][nFish][3], // TS0 ini in BR
            SV_angleOff_CFD_XYZ[nTimeSteps][nFish][2],
            SV_angleOff_SV_XY[nTimeSteps][nFish],   // TS0 ini in BR
            SV_vo_FV_XYZ[nTimeSteps][nFish][3],     // off Flow Vector, allTS ini in BR
            SV_vo_SV_XYZ[nTimeSteps][nFish][3],     // off Swim Vector
            SV_ao_FV_XY[nTimeSteps][nFish],
            agtDetectMetric[nTimeSteps][nFish][nAgents],
            agtUtility[nTimeSteps][nFish][nAgents],
            fishSP_angleOff_SV[FSOlimit-2],
            SP_XYZ_relative[FSOlimit][3],
            sAvg_NPm1[nFish][3],
            t_same_spot[nFish],
            agtCostAvg[nFish],
            k_M[nFish],
            k_F[nFish],
            statesStatistics[nFish][8],
            wallDistances[FSOlimit],
            motivationSum[nFish],
            motivationMean[nFish],
            fatigueSum[nFish],
            fatigueMean[nFish],
            fishBodyLengths[nFish];

// === Open data files for output - Fortran format ============================

    if (outputVfiles >= 1.0) { extraDiagnostics = true; }
    if (outputVfiles >= 2.0) { writeToDebugFile = true; }
    openoutputfiles_(&writeToDebugFile, longResults);


// === First timestep: Initialization =========================================

    timeStep = 0;

    // Environmental mesh field variable values object
    HydroInterpolation OFinterp;
    // Sensory points object
    SensoryOvoid ovoid1;

    if (!args.checkRootCase())   // Mimics "setRootCase.H"
    {
        Foam::FatalError.exit();
    }

    // "Lazy" headers, used by OF to outsource code
    #include "createTime.H"     // Solution control
    #include "createMesh.H"     // Create "mesh" object

    // Create object for tree-based search(es)
    // could be created with smaller boundBox, maybe faster for large mesh (5000000+ ?)
    // polyMesh.H:
    // FACE_PLANES,      //- Faces considered as planes
    // FACE_CENTRE_TRIS, //- Faces decomposed into triangles using face-centre
    // FACE_DIAG_TRIS,   //- Faces decomposed into triangles diagonally
    // CELL_TETS         //- Cell decomposed into tets
    const polyMesh::cellDecomposition decompMode = polyMesh::FACE_CENTRE_TRIS;
    const Foam::meshSearch searchEngine1(mesh, decompMode);

    // Latest timestep: Read selected field variables
    const instantList& times = runTime.times();
    const instant& lastTime = times[times.size()-1];
    #include "createFields.H"

    // Interpolate fields in advance:
    cout << "Interpolating CFD model fields\n" << endl;

    interpolationCellPoint<scalar> triangulateCellsAlpha(alpha);
    // interpolationCellPoint<scalar> triangulateCellsP(p);
    interpolationCellPoint<scalar> triangulateCellsAcclMag(acclMag);
    // interpolationCellPoint<scalar> triangulateCellsGradUsum(gradUsum);
    interpolationCellPoint<scalar> triangulateCellsK(k);
    interpolationCellPoint<Foam::vector> triangulateCellsUMean(UMean);

    cout << "Reading fish positions\n" << endl;
    fileName = "input/fishPositions.inp";
    readFishPositions(
        fileName,
        fishLocations[timeStep],
        nFish, // pass by value
        fishBodyLengths, // decay to pointer
        &nCoeff,
        coefficients,
        k_M,
        k_F
    );

    // no movement from TS=0 to TS=1, so keep initial position
    for (int i = 0; i < nFish*3; i++)
    {
        *(&fishLocations[0][0][0]+nFish*3+i) = *(&fishLocations[0][0][0]+i);
    }

    // Initialize arrays to 0
    initializeArray3D<double>(nTimeSteps, nFish, 2, &SV_angleOff_CFD_XYZ[0][0][0]);
    initializeArray3D<double>(nFish,FSOlimit,3,&fishSensoryVelocity[0][0][0]);
    initializeArray3D<double>(nFish, FSOlimit, nTecFieldVar, &fishSensoryFieldVars[0][0][0]);
    initializeArray3D<double>(nTimeSteps, nFish, nAgents, &agtProb[0][0][0]);
    initializeArray3D<double>(nTimeSteps+1, nFish, nAgents, &agtDetctMetrcAmb[0][0][0]); // values set by BR.f90
    initializeArray3D<double>(FSOlimit,3,1,&SP_XYZ_relative[0][0]);
    initializeArray3D<double>(nFish,3,1,&sAvg_NPm1[0][0]);
    initializeArray3D<double>(nFish,1,1,&t_same_spot[0]);
    initializeArray3D<double>(nFish,1,1,&agtCostAvg[0]);
    initializeArray3D<double>(nFish,7,1,&statesStatistics[0][0]);
    initializeArray3D<double>(FSOlimit,1,1,&wallDistances[0]); // 1D
    initializeArray3D<int>(nTimeSteps,nFish,FSOlimit,&SPfound[0][0][0]);
    initializeArray3D<int>(nTimeSteps,nFish,nAgents,&agtDecision[0][0][0]);
    initializeArray3D<int>(nTimeSteps,nFish,2,&validSVorient_XYZ[0][0][0]);
    initializeArray3D<int>(nTimeSteps, nFish, 1, &fishAttribute[0][0][0]);
    initializeArray3D<int>(nFish,1,1,&numDecisions[0]);


// === Initialize single fish values ==========================================
    cout << "Initializing time step #" << timeStep+1 << "...\n" << endl;

    for (fishNumber = 0; fishNumber < nFish; fishNumber++)
    {
        if (extraDiagnostics)
        {
            cout << "Fish #" << fishNumber+1 << endl;
        }

        xFish = fishLocations[timeStep][fishNumber][0];
        yFish = fishLocations[timeStep][fishNumber][1];
        zFish = fishLocations[timeStep][fishNumber][2];

        Foam::point fishPosition = Foam::point(xFish, yFish, zFish);
        if (!searchEngine1.isInside(fishPosition)) // about the same speed as findCell()
        {
            cout << "ERROR: Out-of-Bound Release: Fish #" << fishNumber+1 << endl;
            cout << "Position: " << xFish << " " << yFish << " " << zFish << endl;
            exit(EXIT_FAILURE);
        }

        // FSO(center) is sufficient for TS=0, no FSO loop
        fishSensoryLocations[fishNumber][0][0] = xFish;
        fishSensoryLocations[fishNumber][0][1] = yFish;
        fishSensoryLocations[fishNumber][0][2] = zFish;

        sAvg_NPm1[fishNumber][0] = xFish;
        sAvg_NPm1[fishNumber][1] = yFish;
        sAvg_NPm1[fishNumber][2] = zFish;

        // Downstream = 0.0; Upstream = 180.0;
        SV_angleOff_CFD_XYZ[timeStep][fishNumber][0] = 180.0;

        OFinterp.interpolateAtFSO
        (
            mesh,
            searchEngine1,
            decompMode,
            timeStep,
            &nTecFieldVar,
            &FSOlimit,
            fishSensoryLocations[fishNumber],
            fishSensoryVelocity[fishNumber],
            fishSensoryFieldVars[fishNumber],
            SPfound[timeStep][fishNumber],
            triangulateCellsAlpha,
//            triangulateCellsP,
            triangulateCellsAcclMag,
//            triangulateCellsGradUsum,
            triangulateCellsK,
            triangulateCellsUMean,
            accelOn,
            TKEconstOn,
            extraDiagnostics
        );

        statesStatistics[fishNumber][0] = fishNumber+1;     // conversion int>double
        statesStatistics[fishNumber][2] = coefficients[38]; // initial motivation (min)
        motivationSum[fishNumber] = coefficients[38];       // initial motivation
        fatigueSum[fishNumber]    = agtCostAvg[fishNumber]; // initial fatigue

        // Initialize variables
        for (int agent = 0; agent < nAgents; agent++)
        {
            agtDecision[timeStep][fishNumber][agent] = -999;   // Allows subroutine BehaviorRule to initialize needed variables
        }
        validSVorient_XYZ[timeStep][fishNumber][0] = 1;   // Allows subroutine SensoryPtCreate to initialize sensory points

    }   // end for fishNumber

// === Prepare wall hits and exits count array ================================

    const Foam::fvPatchList& patches = mesh.boundary();
    int numBoundaries = patches.size();             // count defaultFaces
    int boundaryHitCount_Type[numBoundaries][2];    // Type: 0 = wall, 1 = exit
    numExitRoutes = 0;

    initializeArray3D<int>(numBoundaries, 2, 1, &boundaryHitCount_Type[0][0]);

    forAll(patches, patchI)
    {
        if (patchI == 0) // skip defaultFaces
        {
            continue;
        }
        const polyPatch& polyPatch1 = patches[patchI].patch();
        if (polyPatch1.type() == "patch")
        {
            boundaryHitCount_Type[patchI][1] = 1;
            numExitRoutes++;
        }
    }



// === Further timesteps (> 0) ================================================
    cout << "Starting time loop...\n" << endl;

    // C++ index count: 0..(length-1)
    for (timeStep = 1; timeStep < nTimeSteps; timeStep++)
    {
        if (extraDiagnostics)
        {
            cout << "Time step #" << timeStep+1 << " is starting, t = "
                << timeStep*dt << " s." << endl;
        }

        NiterSinceOut++;

        for (fishNumber = 0; fishNumber < nFish; fishNumber++)
        {
            if (fishAttribute[timeStep][fishNumber][0] == 0) // fish not left
            {

            if (extraDiagnostics)
            {
                cout << "Fish #" << fishNumber+1 << " is starting." << endl;
                cout << "  Calling createSensoryPoints()" << endl;
            }

//            clock_gettime(CLOCK_MONOTONIC, &tStart);

            ovoid1.createSensoryPoints
            (
                &FSOlimit,
                // next line in Fortran: [FN][FSO=0=fish center][acclMag]
                fishSensoryFieldVars[fishNumber][0][2],
                &nCoeff,
                coefficients,
                fishBodyLengths[fishNumber],
                validSVorient_XYZ[timeStep-1][fishNumber][0],
                SV_angleOff_CFD_XYZ[timeStep-1][fishNumber][0],
                ovoidLength,
                fishSP_angleOff_SV,
                SP_XYZ_relative //[FSOlimit][3]
            );

            if (writeToDebugFile)
            {
                writesensorylocations_
                (
                    &fishNumber,
                    &timeStep,
                    &FSOlimit,
                    &ovoidLength,
                    SP_XYZ_relative
                );
            }

            xFish = fishLocations[timeStep][fishNumber][0];
            yFish = fishLocations[timeStep][fishNumber][1];
            zFish = fishLocations[timeStep][fishNumber][2];

            // Compute FSO locations
            for (FSO = 0; FSO < FSOlimit; FSO++)
            {
                fishSensoryLocations[fishNumber][FSO][0] =
                    xFish + SP_XYZ_relative[FSO][0];
                fishSensoryLocations[fishNumber][FSO][1] =
                    yFish + SP_XYZ_relative[FSO][1];
                fishSensoryLocations[fishNumber][FSO][2] =
                    zFish + SP_XYZ_relative[FSO][2];
            } // end for FSO

//            clock_gettime(CLOCK_MONOTONIC, &tEnd);
//            elapsed = (tEnd.tv_sec - tStart.tv_sec);
//            elapsed += (tEnd.tv_nsec - tStart.tv_nsec) / 1000000000.0;
//            cout << elapsed << " s " << endl;
//            clock_gettime(CLOCK_MONOTONIC, &tStart);

            if (extraDiagnostics)
            {
                cout << "  Calling interpolateAtFSO()" << endl;
            }

            OFinterp.interpolateAtFSO
            (
                mesh,
                searchEngine1,
                decompMode,
                timeStep,
                &nTecFieldVar,
                &FSOlimit,
                fishSensoryLocations[fishNumber],
                fishSensoryVelocity[fishNumber],
                fishSensoryFieldVars[fishNumber],
                SPfound[timeStep][fishNumber],
                triangulateCellsAlpha,
//                triangulateCellsP,
                triangulateCellsAcclMag,
//                triangulateCellsGradUsum,
                triangulateCellsK,
                triangulateCellsUMean,
                accelOn,
                TKEconstOn,
                extraDiagnostics
            );

            if (extraDiagnostics)
            {
                cout << "  Calling getWallDistances()" << endl;
            }

            OFinterp.getWallDistances
            (
                mesh,
                searchEngine1,
                xFish, yFish, zFish,
                coefficients[13],   // wallDetectionRange / visualRange
                &FSOlimit,
                SP_XYZ_relative,    // [FSOlimit][3]
                wallDistances       // [FSOlimit]
            );

//            clock_gettime(CLOCK_MONOTONIC, &tEnd);
//            elapsed = (tEnd.tv_sec - tStart.tv_sec);
//            elapsed += (tEnd.tv_nsec - tStart.tv_nsec) / 1000000000.0;
//            cout << elapsed << " s " << endl;
//            clock_gettime(CLOCK_MONOTONIC, &tStart);

            if (extraDiagnostics)
            {
                cout << "  Calling behaviorRule()" << endl;
            }

            behaviorrule_
            (
            // All arguments in and to FORTRAN are passed by reference.
            // Integers, doubles, and bools, need the address-of operator "&", arrays dont
                &fishNumber,
                &nFish,
                &timeStep,
                &dt,
                &nAgents,
                &nTecFieldVar,
                &FSOlimit,
                SPfound[timeStep][fishNumber],                  // send partial arrays
                fishSensoryLocations[fishNumber],               // pointer at 2D array within a 3D
//                fishLocations[timeStep-1],
                fishSensoryVelocity[fishNumber],
                fishSensoryFieldVars[fishNumber],
                &seed,
                &fishSpeedResultant[timeStep-1][fishNumber],    // scalar access, NPm1
                &fishSpeedResultant[timeStep][fishNumber],
                SV_velocityOff_CFD_XYZ[timeStep-1][fishNumber],
                SV_velocityOff_CFD_XYZ[timeStep][fishNumber],
                SV_angleOff_CFD_XYZ[timeStep-1],                // pointer at 2D array
                SV_angleOff_CFD_XYZ[timeStep][fishNumber],
                &SV_angleOff_SV_XY[timeStep-1][fishNumber],     // scalar, last TS
                &SV_angleOff_SV_XY[timeStep][fishNumber],
                agtProb[timeStep-1][fishNumber],                // pointer at 1D array within a 3D
                agtProb[timeStep][fishNumber],                  // Motivation M
                agtDetctMetrcAmb[timeStep][fishNumber],
                agtDetctMetrcAmb[timeStep+1][fishNumber],
                agtDecision,
                &numDecisions[fishNumber],
                validSVorient_XYZ[timeStep][fishNumber],
                &extraDiagnostics,
                &writeToDebugFile,
                &nCoeff,
                coefficients,
                SV_vo_FV_XYZ[timeStep][fishNumber],
                SV_vo_SV_XYZ[timeStep][fishNumber],
                &SV_ao_FV_XY[timeStep][fishNumber],             // scalar within a 3D array
                agtDetectMetric[timeStep][fishNumber],
                agtUtility[timeStep][fishNumber],
                fishSP_angleOff_SV,
                fishAttribute[timeStep],
                &nTimeSteps,
                sAvg_NPm1[fishNumber],
                t_same_spot,
                &ovoidLength,
                &agtCostAvg[fishNumber],                        // Fatigue F
                &k_M[fishNumber],
                &k_F[fishNumber],
                wallDistances,
                &accelOn,
                &visualOn,
                &lowVelOn,
                &TKEconstOn,
                &fishBodyLengths[fishNumber]
            );

//            clock_gettime(CLOCK_MONOTONIC, &tEnd);
//            elapsed = (tEnd.tv_sec - tStart.tv_sec);
//            elapsed += (tEnd.tv_nsec - tStart.tv_nsec) / 1000000000.0;
//            cout << elapsed << " s " << endl;

            // Gather statistics
            // 0=FN, 1=exit time, 2=min M, 3=max M, 4=mean M, 5=min F, 6=max F, 7=mean F
            statesStatistics[fishNumber][1] = timeStep * dt;
            // M
            if (*agtProb[timeStep][fishNumber] < statesStatistics[fishNumber][2])
               {statesStatistics[fishNumber][2] = *agtProb[timeStep][fishNumber];}
            if (*agtProb[timeStep][fishNumber] > statesStatistics[fishNumber][3])
               {statesStatistics[fishNumber][3] = *agtProb[timeStep][fishNumber];}
            motivationSum[fishNumber] += *agtProb[timeStep][fishNumber];
            motivationMean[fishNumber] = motivationSum[fishNumber]/(timeStep+1);
            statesStatistics[fishNumber][4] = motivationMean[fishNumber];
            // F
            if (agtCostAvg[fishNumber] < statesStatistics[fishNumber][5])
               {statesStatistics[fishNumber][5] = agtCostAvg[fishNumber];}
            if (agtCostAvg[fishNumber] > statesStatistics[fishNumber][6])
               {statesStatistics[fishNumber][6] = agtCostAvg[fishNumber];}
            fatigueSum[fishNumber] += agtCostAvg[fishNumber];
            fatigueMean[fishNumber] = fatigueSum[fishNumber]/(timeStep+1);
            statesStatistics[fishNumber][7] = fatigueMean[fishNumber];
            // Gather statistics end

            }
            else   // fish left; invertebrates ignored just like in ELAM-2014
            {
                if (extraDiagnostics)
                {
                    cout << "Fish #" << fishNumber+1 << " skipped." << endl;
                }

                // Prevent undefined behavior - put 0.0 to left fishes
                for (int i = 0; i < 3; i++)
                {
                    SV_velocityOff_CFD_XYZ[timeStep][fishNumber][i] = 0.0;
                    SV_vo_FV_XYZ[timeStep][fishNumber][i] = 0.0;
                    SV_vo_SV_XYZ[timeStep][fishNumber][i] = 0.0;
                }
                SV_ao_FV_XY[timeStep][fishNumber]       = 0.0;
                SV_angleOff_SV_XY[timeStep][fishNumber] = 0.0;
            } // end if fishAttribute == 0

            if (extraDiagnostics)
            {
                cout << "  Calling updateFishLocation()" << endl;
            }

//            clock_gettime(CLOCK_MONOTONIC, &tStart);

            // updateFishLocation is outside the fishAttribute "if" brackets
            // because it is needed to "transport" (frozen) coordinates of left
            // fishes through time
            updateFishLocation
            (
                mesh,
                searchEngine1,
                decompMode,
                fishNumber,
                passiveTransport,
                dt,
                fishLocations[timeStep+1][fishNumber][0],   // xFish new
                fishLocations[timeStep+1][fishNumber][1],   // yFish new
                fishLocations[timeStep+1][fishNumber][2],   // zFish new
                fishLocations[timeStep][fishNumber][0],     // xFish
                fishLocations[timeStep][fishNumber][1],     // yFish
                fishLocations[timeStep][fishNumber][2],     // zFish
                fishSensoryVelocity[fishNumber][0][0],      // u at FSO center
                fishSensoryVelocity[fishNumber][0][1],      // v at FSO center
                fishSensoryVelocity[fishNumber][0][2],      // w at FSO center
                SV_velocityOff_CFD_XYZ[timeStep][fishNumber][0],   // uFish
                SV_velocityOff_CFD_XYZ[timeStep][fishNumber][1],   // vFish
                SV_velocityOff_CFD_XYZ[timeStep][fishNumber][2],   // wFish
                fishAttribute[timeStep][fishNumber][0],
                fishAttribute[timeStep+1][fishNumber][0],
                boundaryHitCount_Type,
                extraDiagnostics,
                triangulateCellsAlpha
            );

//            clock_gettime(CLOCK_MONOTONIC, &tEnd);
//            elapsed = (tEnd.tv_sec - tStart.tv_sec);
//            elapsed += (tEnd.tv_nsec - tStart.tv_nsec) / 1000000000.0;
//            cout << elapsed << " s " << endl;

            totalExitTally = sumArray2Dconditionally
            (
                numBoundaries, &boundaryHitCount_Type[0][0]
            );

            if ((NiterSinceOut == outIntV) && (longResults == 1))
            {
                outputfishdata_zonesaretime_
                (
                    &fishNumber,
                    &nFish,
                    &timeStep,
                    &nTimeSteps,
                    &dt,
                    &FSOlimit,
                    &nAgents,
                    fishLocations,
                    SV_velocityOff_CFD_XYZ,
                    agtUtility,
                    fishSpeedResultant,
                    fishAttribute,
                    &totalExitTally,
                    fishSensoryVelocity
                );
            }
        } // end for fishNumber (FN)

        if ((NiterSinceOut == outIntV) || (nFish-totalExitTally <= 0))
        {
            outputfishpassageanddecisions_
            (
                &nFish,
                &timeStep,
                &dt,
                &numExitRoutes,
                &numBoundaries,
                boundaryHitCount_Type
            );
        }

        // Omit outIntV timesteps from output
        if (NiterSinceOut == outIntV)
        {
            NiterSinceOut = 0;
        }

        // Exit TS loop if all fish have left the system
        if (nFish-totalExitTally <= 0)
        {
            break;
        }

        if (extraDiagnostics)
        {
            cout << endl;
        }

    } // end for timeStep (TS)

    int numDecisionsSum = 0;
    for(auto& num_i : numDecisions)
        numDecisionsSum += num_i;

    endofrunoutput_
    (
        &nFish,
        &timeStep,
        &nTimeSteps,
        &dt,
        numDecisions,
        &outIntV,
        fishLocations,
        SV_velocityOff_CFD_XYZ,
        fishSpeedResultant,
        fishAttribute,
        SV_angleOff_CFD_XYZ,
        longResults,
        statesStatistics,
        numDecisionsSum
    );

    cout << endl;
    cout << "End of run at t = " << timeStep*dt << " s." << endl;
    cout << endl;

    exit(EXIT_SUCCESS);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
