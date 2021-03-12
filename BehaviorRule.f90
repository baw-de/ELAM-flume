!/*--------------------------------------------------------------------*\
! Application
!   BehaviorRule
!
! Purpose
!   Determine volitional swim direction, horizontally and vertical, and
!   swim speed.
!
! Portions of BehaviorRule are
!   Copyright (c) 2021 David Gisen, Bundesanstalt fuer Wasserbau
!
! Further portions of BehaviorRule are
!   provided courtesy of the United States Army Corps of Engineers,
!   Engineer Research and Development Center, Environmental Laboratory
!
! Portions are indicated in the code by respective headings. See license
!   file for information on usage.
!\*--------------------------------------------------------------------*/

!/*--------------------------------------------------------------------*\
! The following portions are provided courtesy of the United States Army
!   Corps of Engineers, Engineer Research and Development Center,
!   Environmental Laboratory. See license file for information on usage.
!\*--------------------------------------------------------------------*/

                                                                         !                          Equivalence in DRIVER code:
      Subroutine BehaviorRule(
     I                        FishNumber,                                ! Local scalar             Fish number (FN)
     I                        NFish,                                     ! Local scalar
     I                        TimeStep,                                  ! Local scalar             Time step
     I                        dt,                                        ! Local scalar             TimeStepLength
     I                        NAgents,                                   ! Local scalar
     I                        NTecFieldVar,                              ! Local scalar
     I                        FSOLIMIT,                                  ! Global parameter
     I                        SPFound_NP,                                ! Local array              SPFound(1,1:FSOLIMIT,FN,NP)
     I                        FishSensoryLocation_NP,                    ! Local array              FishSensoryLocation_NP(1:3,1:FSOLIMIT,FN)
!     I                        FishLocation_NPm1,                         ! Local array              FishLocation(1:3,1:NFish,NPorTS-1)
!     I                        FishLocation_NP,                           ! Local array              FishLocation(1:3,1:NFish,NPorTS)
     I                        FishSensoryVelocity_NP,                    ! Local array              FishSensoryVelocity_NP(1:3,1:FSOLIMIT,FN)
     I                        FishSensoryFieldVars_NP,                   ! Local array              FishSensoryFieldVars_NP(1:NTecFieldVar,1:FSOLIMIT,FN)
     I                        SEED,                                      ! Local scalar
     I                        FhSpdRes_NPm1,                             ! Local scalar             FhSpdRes(FN,NP-1)
     O                        FhSpdRes_NP,                               ! Local scalar             FhSpdRes(FN,NP)
     I                        SVvoCFDXYZ_NPm1,                           ! Local array              1st position => VelFishCFD(1,FN,NP-1)
                                                                         !                          2nd position => VelFishCFD(2,FN,NP-1)
                                                                         !                          3rd position => VelFishCFD(3,FN,NP-1)
     O                        SVvoCFDXYZ_NP,                             ! Local array              1st position => VelFishCFD(1,FN,NP)
                                                                         !                          2nd position => VelFishCFD(2,FN,NP)
                                                                         !                          3rd position => VelFishCFD(3,FN,NP)
     I                        SVaoCFDXYZ_NPm1,                           ! Local array              1st position => SVAOCFDXYZ(1,1:NFish,NP-1)
                                                                         !                          2nd position => SVAOCFDXYZ(2,1:NFish,NP-1)
     O                        SVaoCFDXYZ_NP,                             ! Local array              1st position => SVAOCFDXYZ(1,FN,NP)
                                                                         !                          2nd position => SVAOCFDXYZ(2,FN,NP)
     I                        SVaoSVXY_NPm1,                             ! Local scalar             SVOSVXY(FN,NP-1)
     O                        SVaoSVXY_NP,                               ! Local scalar             SVOSVXY(FN,NP)
     I                        AgtProb_NPm1,                              ! Local array              AgtProb(1:NAgents,FN,NP-1)
     O                        AgtProb_NP,                                ! Local array              AgtProb(1:NAgents,FN,NP)
     I                        AgtDetctMetrcAmb_NP,                       ! Local array              AgtDetctMetrcAmb(1:NAgents,FN,NP)
     O                        AgtDetctMetrcAmb_NPp1,                     ! Local array              AgtDetctMetrcAmb(1:NAgents,FN,NP+1)
     I                        AgtDecision,                               ! Local array              AgtDecision(1:NAgents,1:NFish,1:TotTimeSteps)
     O                        NumDecisions_NP,                           ! Local scalar             NumDecisions(FN)
     O                        VldSVOrientXYZ_NP,                         ! Local array              VldSVOrientXYZ(1:2,FN,NP)
     I                        ExtraDiagnostics,                          ! Local logical            ExtraDiagnostics
     I                        WriteToDebugFile,                          ! Local logical            ExtraDiagnostics
     I                        nCoeff,
     I                        Coefficients,                              ! Local array              Coefficients(1:nCoeff)
     O                        SVvoFVXYZ_NP,                              ! Local array              VelFishPTV(1:3,FN,NP)
     O                        SVvoSVXYZ_NP,                              ! Local array              VelFishPSV(1:3,FN,NP)
     O                        SVaoFVXY_NP,                               ! Local scalar             SVOFVXY(FN,NP)
!     O                        AgtDetctThrshld,                           ! Local array              AgtDetctThrshld(1:NAgents)
     O                        AgtDetctMetrc_NP,                          ! Local array              AgtDetctMetrc(1:NAgents,FN,NP)
     O                        AgtUtil_NP,                                ! Local array              AgtUtil(1:NAgents,FN,NP)
     I                        FhSPaoSV,                                  ! Local array              FhSPAOSV(1:5)
!     I                        FhAttrb_NPm1,                              ! Local array              1st position => FhAttrb(1:nFhAttrb,1:NFish,NP-1)
     O                        FhAttrb_NP,                                !                          1st position => FhAttrb(1:nFhAttrb,1:NFish,NP)
     I                        TotTimeSteps,                              ! Local scalar
     I                        sAvg_NPm1,                                 ! Local array (new)
     I                        t_same_spot,                               ! Local array
     I                        ovoidLength,                               ! Local scalar
     I                        AgtCostAvg_NPm1,
     I                        k_M,
     I                        k_F,
     I                        wallDistances,
     I                        accelOn,
     I                        visualOn,
     I                        lowVelOn,
     I                        TKEconstOn,
     I                        fishBodyLength
     .                        )

      implicit none

      integer :: NAgents, NFish                                          ! Scalars (for dimensioning arrays)
      integer :: TotTimeSteps                                            ! Scalars (for dimensioning arrays)
      integer, intent(in) :: FSOLIMIT, nCoeff                            ! Scalars (const)

      integer*4 :: SPFound_NP(FSOLIMIT),                                 ! Arrays
     .        VldSVOrientXYZ_NP(2),
     .        VldSVOrientXYZ_NPm1(2),
     .        VldFVOrientXYZ_NP(2),
     .        AgtDecision(NAgents,NFish,TotTimeSteps),
     .        AgtEvent_NP(NAgents),
     .        NAgtEvents_NP(NAgents),
     .        AgentChosen(2),
     .        VldRVOrientXYZ(2),
     .        VldMVOrientXYZ(2),
     .        FhAttrb_NPm1(1,NFish),
     .        FhAttrb_NP(1,NFish),
     .        FishNumber,                                                ! Scalars
     .        TimeStep,
     .        NTecFieldVar,
     .        NumDecisions_NP,
     .        FSO,
     .        SEED,
     .        VectorType,
     .        dim

      real*8 :: FishSensoryVelocity_NP(3,FSOLIMIT),                      ! Arrays
     .        FishSensoryFieldVars_NP(NTecFieldVar,FSOLIMIT),
     .        SVvoCFDXYZ_NPm1(3),
     .        SVvoCFDXYZ_NP(3),
     .        SVaoCFDXYZ_NPm1(2,NFish),
     .        SVaoCFDXYZ_NP(2),
     .        AgtProb_NPm1(NAgents),
     .        AgtProb_NP(NAgents),
     .        AgtDetctThrshld(NAgents),
     .        AgtDetctMetrcAmbMem(NAgents),
     .        AgtDetctMetrc_NP(NAgents),
     .        AgtDetctMetrcAmb_NP(NAgents),
     .        AgtDetctMetrcAmb_NPp1(NAgents),
     .        AgtCost_NP(NAgents),
     .        AgtUtil_NP(NAgents),
     .        AgtMem(NAgents),
!     .        AgtIntUtil(NAgents),
     .        SVvoFVXYZ_NP(3),
     .        SVvoSVXYZ_NP(3),
     .        FVaoSVXY_NP(FSOLIMIT),
     .        FVaoCFDXYZ_NP(2,FSOLIMIT),
     .        FVvoSVXYZ_NP(3,FSOLIMIT),
     .        CondVelM(FSOLIMIT),
     .        CondAcclM(FSOLIMIT),
     .        Coefficients(nCoeff),
     .        FhSPaoSV(5),                                                ! unused
     .        MVComponentsCFD(3),
     .        RVaoCFDXYZ(2),
     .        MVaoCFDXYZ(2),
     .        RVvoCFDXYZ(3),
     .        MVvoRVXYZ(3),
     .        FishSensoryLocation_NP(3,FSOLIMIT),
!     .          FishLocation_NPm1(3,NFish),
!     .          FishLocation_NP(3,NFish),
     .        t_same_spot(NFish),
     .        sAvg(3),sAvg_NPm1(3),
     .        TKE(FSOLIMIT),
     .        FishEleva_NP,                                               ! Scalars
     .        RRR,
     .        FhSpdRes_NPm1,
     .        FhSpdRes_NP,
     .        SVaoSVXY_NPm1, SVaoSVXY_NP,
     .        SVaoFVXY_NP,SVaoFVXYTemp,SVaoFVXYRad,
     .        MVaoRVXY,
     .        SVaoCFDXYTemp,SVaoCFDXYRad,
     .        STRBASELOG,
     .        dt,
     .        RAND1,
     .        SVaoCFDZRad,FhSpdResXY,SVaoSVXYRad,
     .        FHBODYLENGTH,
     .        dim_diff_mag,
     .        ovoidLength,
     .        acclMagTolMigr, velMagTol,
     .        deltaHolding,
     .        dRmag, m_s, m_F, rSameSpot,
     .        AgtCostAvg, AgtCostAvg_NPm1,
     .        k_M, k_F,
     .        wallDistances(FSOLIMIT),
     .        TKEdiffLeft, TKEdiffMid, TKEdiffRight,
     .        driftStraightProb,
     .        fishBodyLength

      logical*1 :: WriteToDebugFile, ExtraDiagnostics,
     .        accelOn, visualOn, lowVelOn, TKEconstOn



!***************************************************************************************************
!***************************************************************************************************
!***************************************    SECTION 1.1    *****************************************
!***************************************************************************************************
!******************************** Do Not Modify Between Lines **************************************

! Initialize variables used in subroutine BehaviorRule.

! Account for different array indexing in Fortran/C++ - to be reversed at file end
      FishNumber = FishNumber+1
      TimeStep = TimeStep+1

!      FishEleva_NPm1 = REAL(FishLocation_NPm1(3,FishNumber))            ! unused
      FishEleva_NP   = real(FishSensoryLocation_NP(3,1))                 ! Fish center z value

!      WriteToDebugFile = .TRUE.
      IF (INT(FhAttrb_NP(1,FishNumber)) /= 0) THEN                       ! Attribute #1: 0=Fish ; 1=Invertebrate ; -1=To be removed from simulation
        WriteToDebugFile = .FALSE.                                       ! Do not write for fish which have exited
      END IF

      DO FSO=1,FSOLIMIT
        IF (FishSensoryFieldVars_NP(3,FSO) < 2E-6) THEN                  ! 2E-6 is arbitrary but needed value to prevent negative #s
          FishSensoryFieldVars_NP(3,FSO) = 2E-6                          ! Treatment just for acclMag - log scale
        END IF
      END DO

      IF (AgtDecision(1,FishNumber,TimeStep-1) == -999) THEN
        AgtDecision(1:NAgents,FishNumber,TimeStep-1) = 0
        SVaoSVXY_NPm1                                = 0.0
        SVaoCFDXYZ_NPm1(1:2,FishNumber)              = 0.0
        SVvoCFDXYZ_NPm1(1:3)                         = 0.0
        FhSpdRes_NPm1                                = 0.0
        AgtDetctMetrcAmb_NP(4)                       = FishEleva_NP
        AgtProb_NPm1(1:NAgents)                      = Coefficients(39)   ! Initial motivation
      END IF

      SVaoFVXY_NP                                = 0.0
      SVvoFVXYZ_NP(1:3)                          = 0.0
      SVvoSVXYZ_NP(1:3)                          = 0.0
      VldSVOrientXYZ_NPm1(1:2)                   = 0
      VldSVOrientXYZ_NP(1:2)                     = 0
      VldFVOrientXYZ_NP(1:2)                     = 0

      ! 3 VectorRelation variables
      FVaoSVXY_NP(1:FSOLIMIT)                    = 0.0
      FVaoCFDXYZ_NP(1:2,1:FSOLIMIT)              = 0.0
      FVvoSVXYZ_NP(1:3,1:FSOLIMIT)               = 0.0

      AgtDecision(1:NAgents,FishNumber,TimeStep) = 0
      AgtEvent_NP(1:NAgents)                     = 0
      NAgtEvents_NP(1:NAgents)                   = 1
      AgtDetctMetrc_NP(1:NAgents)                = 0.0
      AgtCost_NP(1:NAgents)                      = 0.0
      AgtUtil_NP(1:NAgents)                      = 0.0
      AgtMem(1:NAgents)                          = 0.0
!      AgtIntUtil(1:NAgents)                      = 0.0
      AgtDetctThrshld(1:NAgents)                 = 0.0
      AgtDetctMetrcAmbMem(1:NAgents)             = 0.0
      NumDecisions_NP                            = NumDecisions_NP + 1
      dim_diff_mag                               = 0.0
      AgentChosen(1:2)                           = 0


!********************************* Do Not Modify Between Lines **************************************
!****************************************************************************************************
!****************************************    SECTION 1.2    *****************************************
!****************************************************************************************************
!****************************************************************************************************

!Rules Developed Assume 7 Sensory Points: One at the Fish Centroid and in Each Coordinate Direction.

!Agent Parameters.

!     1 = Migrating - Swim against flow vector (default)
!     2 = Holding - Stand against flow vector
!     3 = Drifting - Drift with flow vector
!     4 = Swim toward acclimatized depth


! Memory coefficients:
      AgtMem(1) = Coefficients(5)                                        ! Memory coefficient for Agt 1
      AgtMem(2) = Coefficients(6)                                        ! Memory coefficient for Agt 2
      AgtMem(3) = Coefficients(7)                                        ! Memory coefficient for Agt 3
      AgtMem(4) = Coefficients(8)                                        ! Memory coefficient for Agt 4

! Intrinsic utilities:
!      AgtIntUtil(1) = Coefficients(13)                                   ! Intrinsic utility for Agt 1
!      AgtIntUtil(2) = Coefficients(14)                                   ! Intrinsic utility for Agt 2
!      AgtIntUtil(3) = Coefficients(15)                                   ! Intrinsic utility for Agt 3
!      AgtIntUtil(4) = Coefficients(16)                                   ! Intrinsic utility for Agt 4

! Stimulus thresholds:
!      AgtDetctThrshld(1) = Coefficients(21)                              !
!      AgtDetctThrshld(2) = Coefficients(22)                              !
!      AgtDetctThrshld(3) = Coefficients(23)                              !
      AgtDetctThrshld(4) = Coefficients(24)                              ! Threshold of perceived change in AgtDetctMetrc_NP (Pressure) kB{4}

! Memory coefficients for acclimatized values:
!      AgtDetctMetrcAmbMem(1) = Coefficients(29)                          !
!      AgtDetctMetrcAmbMem(2) = Coefficients(30)                          !
!      AgtDetctMetrcAmbMem(3) = Coefficients(31)                          !
      AgtDetctMetrcAmbMem(4) = Coefficients(32)                          ! Memory coefficient for updating acclimatization to Pressure (Agt 4)

! Set Fish Size
      FHBODYLENGTH      = fishBodyLength                                  ! Length of fish (meters)

! Other Parameters
      STRBASELOG = Coefficients(55)                                      ! Base reference value for calculation of AcclM in decibel scale. In theory,
                                                                         !   this would be analogous to the minimum "threshold of detection" of AcclM in water.


!****************************************************************************************************
!****************************************************************************************************
!****************************************    SECTION 2.1    *****************************************
!****************************************************************************************************
!********************************* Do Not Modify Between Lines **************************************

!Determine Flow Vector Velocity Relative to Swim Vector Velocity.

      IF ((SVvoCFDXYZ_NPm1(1) == 0.0) .AND.
     .    (SVvoCFDXYZ_NPm1(2) == 0.0)) THEN
        VldSVOrientXYZ_NPm1(1) = 0                                       ! No valid xy-plane orientation of fish axis => no fish movement from t-1 to t
      ELSE
        VldSVOrientXYZ_NPm1(1) = 1                                       ! Valid xy-plane orientation of fish axis can be determined for t-1 to t
      END IF

      DO FSO=1,1
        IF (SPFound_NP(FSO) == 0) GOTO 2000                              ! out-of-bounds?
        MVComponentsCFD = FishSensoryVelocity_NP(:,FSO)
        VectorType      = 1                                              ! PositionVector=0 ; VelocityVector=1
        VldRVOrientXYZ  = VldSVOrientXYZ_NPm1
        RVaoCFDXYZ      = SVaoCFDXYZ_NPm1(:,FishNumber)                  ! t-1 (to t) is the most up-to-date info on where fish is pointing
        RVvoCFDXYZ      = SVvoCFDXYZ_NPm1
        CALL VectorRelation(
     I                      VectorType,                                  ! INTEGER => PositionVector=0 ; VelocityVector=1
     I                      MVComponentsCFD,                             ! REAL
     I                      VldRVOrientXYZ,                              ! INTEGER
     O                      VldMVOrientXYZ,                              ! INTEGER
     I                      RVaoCFDXYZ,                                  ! REAL
     I                      RVvoCFDXYZ,                                  ! REAL
     O                      MVaoRVXY,                                    ! REAL    => (M)aster (V)ector (a)ngle (o)ff (R)eference (V)ector's xy-axis
     O                      MVaoCFDXYZ,                                  ! REAL
     O                      MVvoRVXYZ                                    ! REAL    => (M)aster (V)ector (v)elocity (o)ff (R)eference (V)ector's xy-axis
     .                      )
        VldFVOrientXYZ_NP    = VldMVOrientXYZ
        FVaoSVXY_NP(FSO)     = MVaoRVXY                                  ! Flow velocity Vector xy-angle off (previous) Swim Vector axis (i.e., flow velocity xy-angle relative to swim direction) (-180..+180)
        FVaoCFDXYZ_NP(:,FSO) = MVaoCFDXYZ                                ! Flow velocity Vector angle off xz-plane (horizontal, 0..360) and xy-plane (vertical, -90..+90)
        FVvoSVXYZ_NP(:,FSO)  = MVvoRVXYZ                                 ! Flow velocity components off moving Swim Vector axis (i.e., velocity components as if fish is stationary)
 2000   CONTINUE
      END DO   ! FSO


!********************************* Do Not Modify Between Lines **************************************
!****************************************************************************************************
!****************************************    SECTION 2.2    *****************************************
!****************************************************************************************************
!****************************************************************************************************

! Format definitions
  221 format (' ', A, F8.3) ! (empty) control character, character, real # - Positions up to 9999.999 m
  223 format (' ', A, I1)   ! space, character, integer #, up to 99

      IF (ExtraDiagnostics .AND. WriteToDebugFile) THEN
!        WRITE(72,*) ' '
!        WRITE(72,*) '  FN  = ',FishNumber
!        WRITE(72,*) '  TS  = ',TimeStep
        WRITE(72,223) '    SPFound_NP(1) = ',SPFound_NP(1)
        WRITE(72,221) '    x = ',FishSensoryLocation_NP(1,1)
        WRITE(72,221) '    y = ',FishSensoryLocation_NP(2,1)
        WRITE(72,221) '    z = ',FishSensoryLocation_NP(3,1)
        WRITE(72,221) '      SVaoCFDXYZ_NPm1(1,FishNumber) = ',
     .                     SVaoCFDXYZ_NPm1(1,FishNumber)
        WRITE(72,223) '       SPFound_NP(2) = ',SPFound_NP(2)
        WRITE(72,221) '        SP2x = ',FishSensoryLocation_NP(1,2)
        WRITE(72,221) '        SP2y = ',FishSensoryLocation_NP(2,2)
        WRITE(72,221) '        SP2z = ',FishSensoryLocation_NP(3,2)
        WRITE(72,223) '       SPFound_NP(3) = ',SPFound_NP(3)
        WRITE(72,221) '        SP3x = ',FishSensoryLocation_NP(1,3)
        WRITE(72,221) '        SP3y = ',FishSensoryLocation_NP(2,3)
        WRITE(72,221) '        SP3z = ',FishSensoryLocation_NP(3,3)
        WRITE(72,223) '       SPFound_NP(4) = ',SPFound_NP(4)
        WRITE(72,221) '        SP4x = ',FishSensoryLocation_NP(1,4)
        WRITE(72,221) '        SP4y = ',FishSensoryLocation_NP(2,4)
        WRITE(72,221) '        SP4z = ',FishSensoryLocation_NP(3,4)
        WRITE(72,223) '       SPFound_NP(5) = ',SPFound_NP(5)
        WRITE(72,221) '        SP5x = ',FishSensoryLocation_NP(1,5)
        WRITE(72,221) '        SP5y = ',FishSensoryLocation_NP(2,5)
        WRITE(72,221) '        SP5z = ',FishSensoryLocation_NP(3,5)
        WRITE(72,223) '       SPFound_NP(6) = ',SPFound_NP(6)
        WRITE(72,221) '        SP6x = ',FishSensoryLocation_NP(1,6)
        WRITE(72,221) '        SP6y = ',FishSensoryLocation_NP(2,6)
        WRITE(72,221) '        SP6z = ',FishSensoryLocation_NP(3,6)
        WRITE(72,223) '       SPFound_NP(7) = ',SPFound_NP(7)
        WRITE(72,221) '        SP7x = ',FishSensoryLocation_NP(1,7)
        WRITE(72,221) '        SP7y = ',FishSensoryLocation_NP(2,7)
        WRITE(72,221) '        SP7z = ',FishSensoryLocation_NP(3,7)
      END IF


!Set Default Agent Behavior (swim with flow vector).

      AgtEvent_NP(1:NAgents) = 0
      AgtEvent_NP(1)         = 1


! Calculate [VelM], [AcclM] and TKE:
      DO FSO=1,FSOLIMIT
          CondVelM(FSO) = SQRT( FishSensoryVelocity_NP(1,FSO) *
     .                          FishSensoryVelocity_NP(1,FSO) +
     .                          FishSensoryVelocity_NP(2,FSO) *
     .                          FishSensoryVelocity_NP(2,FSO) +
     .                          FishSensoryVelocity_NP(3,FSO) *
     .                          FishSensoryVelocity_NP(3,FSO) )
          CondAcclM(FSO) = FishSensoryFieldVars_NP(3,FSO)
          TKE(FSO)       = FishSensoryFieldVars_NP(2,FSO)
      END DO

!****************************************************************************************************
!****************************************************************************************************
!****************************************    SECTION 2.3    *****************************************
!****************************************************************************************************
! The following portions are
! Copyright (c) 2021 David Gisen, Bundesanstalt fuer Wasserbau
!   See license file for information on usage
!****************************************************************************************************

! Auxiliary variables for Motivation calculation

! DG - Format
!  222 format (' ', A, ES10.3)     ! character, scientific #


!     Calculate time interval being close to the same spot
      m_s = AgtMem(4)
      sAvg = (1-m_s) * FishSensoryLocation_NP(:,1)                      ! sAvg = position average
     .+ m_s * sAvg_NPm1

      sAvg_NPm1 = sAvg
      dRmag = norm2(FishSensoryLocation_NP(:,1) - sAvg)


      rSameSpot = Coefficients(13)

      if ( dRmag >= rSameSpot*FHBODYLENGTH .or.
     .  t_same_spot(FishNumber) > Coefficients(57)+dt ) then
        t_same_spot(FishNumber) = 0.0
      else
        t_same_spot(FishNumber) = t_same_spot(FishNumber)+dt
      end if


!****************************************************************************************************
!****************************************************************************************************
!****************************************    SECTION 2.4    *****************************************
!****************************************************************************************************
!****************************************************************************************************

! Does change in [Elevation/Pressure] exceed threshold:

      AgtDetctMetrc_NP(4) = FishEleva_NP - AgtDetctMetrcAmb_NP(4)
      if (AgtDetctMetrc_NP(4) > AgtDetctThrshld(4)*dt) then
        AgentChosen(2) = 4   ! get down
      elseif (AgtDetctMetrc_NP(4) < -AgtDetctThrshld(4)*dt) then
        AgentChosen(2) = 5   ! get up
      end if

!****************************************************************************************************
!****************************************************************************************************
!****************************************    SECTION 2.5    *****************************************
!****************************************************************************************************
!****************************************************************************************************

!Update Acclimatized Value of [Elevation/Pressure].

      AgtDetctMetrcAmb_NPp1(4) = (1.0 - AgtDetctMetrcAmbMem(4)) *
     .                           FishEleva_NP +
     .                           AgtDetctMetrcAmbMem(4) *
     .                           AgtDetctMetrcAmb_NP(4)


!****************************************************************************************************
!****************************************************************************************************
!****************************************    SECTION 2.6    *****************************************
!****************************************************************************************************
!****************************************************************************************************

!Calculate Bioenergetic Costs Associated with Respond to Each 'Agent' Regardless of Success.

!     Fatigue F
      AgtCost_NP(1) = 1.0/k_F * FhSpdRes_NPm1
     .    / FHBODYLENGTH
!     Cost must be within the range of 0.0 to 1.0:
      if (AgtCost_NP(1) > 1) AgtCost_NP(1) = 1.0


!     Fatigue memory
      if (AgtCost_NP(1) <= AgtCostAvg_NPm1) then
        m_F = AgtMem(2)
      else
        m_F = AgtMem(3)
      endif

      AgtCostAvg = (1.0-m_F) * AgtCost_NP(1)
     .                + m_F * AgtCostAvg_NPm1

      AgtCostAvg_NPm1 = AgtCostAvg


!****************************************************************************************************
!****************************************************************************************************
!****************************************    SECTION 2.7    *****************************************
!****************************************************************************************************
!****************************************************************************************************

!Calculate the Probability and Utility (Motivation) Associated with Each 'Agent'.


!     Motivation M
      AgtUtil_NP(1) = 1.0/k_M * t_same_spot(FishNumber)

!     M must be within the range of 0.0 to 1.0:
      if (AgtUtil_NP(1) > 1.0) AgtUtil_NP(1) = 1.0


!     Motivation memory
      AgtUtil_NP(1) = (1.0-AgtMem(1)) * AgtUtil_NP(1)
     .                + AgtMem(1) * AgtProb_NPm1(1)

      AgtProb_NP(1) = AgtUtil_NP(1)



!     chooseBehavior(): Motivation vs. fatigue
      deltaHolding = Coefficients(56)
      if      (AgtUtil_NP(1) > AgtCostAvg + deltaHolding) then
        !Migrate
        AgentChosen(1) = 1
      else if (AgtUtil_NP(1) < AgtCostAvg - deltaHolding) then
        !Drift
        AgentChosen(1) = 3
      else
        !Hold
        AgentChosen(1) = 2
      end if



      IF (ExtraDiagnostics .AND. WriteToDebugFile) THEN
        WRITE(72,*) ' '
        WRITE(72,*) '   => FhAttrb_NP                         = ',
     .                     FhAttrb_NP(1,FishNumber)
        WRITE(72,*) '   => AgentChosen(1)                     = ',
     .                     AgentChosen(1)
        WRITE(72,*) '      AgentChosen(2) [Vertical Override] = ',
     .                     AgentChosen(2)
        WRITE(72,*) '         AgtEvent_NP(1) = ',
     .                        AgtEvent_NP(1)
        WRITE(72,*) '         AgtEvent_NP(2) = ',
     .                        AgtEvent_NP(2)
        WRITE(72,*) '         AgtEvent_NP(3) = ',
     .                        AgtEvent_NP(3)
        WRITE(72,*) '         AgtEvent_NP(4) = ',
     .                        AgtEvent_NP(4)

        WRITE(72,*) ' '
        WRITE(72,*) '          AgtUtil_NP(1) = ',
     .                         AgtUtil_NP(1)
        WRITE(72,*) '          AgtUtil_NP(2) = ',
     .                         AgtUtil_NP(2)
        WRITE(72,*) '          AgtUtil_NP(3) = ',
     .                         AgtUtil_NP(3)
        WRITE(72,*) '          AgtUtil_NP(4) = ',
     .                         AgtUtil_NP(4)
        WRITE(72,*) ' '
!        WRITE(72,*) ' '
        WRITE(72,*) '          AgtCost_NP(1) = ',
     .                         AgtCostAvg
        WRITE(72,*) '          AgtDetctMetrcAmb_NP(4) = ',
     .                         AgtDetctMetrcAmb_NP(4)
        write(72,*) '          t_same_spot =',
     .                         t_same_spot(FishNumber)
      END IF


!****************************************************************************************************
!****************************************************************************************************
!****************************************    SECTION 3.1    *****************************************
!****************************************************************************************************
!****************************************************************************************************

!At this point in the behavior rule code you must have:
!  - decided which behavior to implement
!
!You will now lay out the code to implement each behavior.
!There are 2 general categories of behavior rules:
!  (1) using a vector as the basis for a behavior response. For instance, behavior
!      could be moving at an angle relative to the direction of flow (e.g., swim with flow vector).
!  (2) using the gradient of a value for a behavior response. For instance, behavior
!      could be moving at an angle relative to the gradient in velocity magnitude,
!      acceleration, temperature, dissolved oxygen, etc.

!****************************************************************************************************

! Behavior implemented here: swim against flow vector (Migrating, default)

      if (AgentChosen(1) == 1) then

!     Assign swim angle:
        ! Horizontal, against flow vector, relative to swim vector:
        SVaoCFDXYTemp = FVaoCFDXYZ_NP(1,1) - 180.0

        ! 0-360 check, robust version
        do while (SVaoCFDXYTemp < 0.0 )
            SVaoCFDXYTemp = SVaoCFDXYTemp + 360.0
        end do
        do while (SVaoCFDXYTemp >= 360.0)
            SVaoCFDXYTemp = SVaoCFDXYTemp - 360.0
        end do

        SVaoSVXY_NP = SVaoCFDXYTemp - SVaoCFDXYZ_NPm1(1,FishNumber)

        ! Towards lower velM; SP 4=left, 5=right
        if (lowVelOn) then
            velMagTol = Coefficients(40)
            if      ( CondVelM(4) < CondVelM(5)/velMagTol ) then
                SVaoSVXY_NP = SVaoSVXY_NP + Coefficients(43)*dt ! go left
            else if ( CondVelM(4) > CondVelM(5)*velMagTol ) then
                SVaoSVXY_NP = SVaoSVXY_NP - Coefficients(43)*dt ! go right
            else
                !stay
            end if
        end if

        ! Towards higher acclmag
        if (accelOn) then
            acclMagTolMigr = Coefficients(40)
            if      ( CondAcclM(4) > CondAcclM(5)*acclMagTolMigr ) then
                SVaoSVXY_NP = SVaoSVXY_NP + Coefficients(45)*dt ! go left
            else if ( CondAcclM(4) < CondAcclM(5)/acclMagTolMigr ) then
                SVaoSVXY_NP = SVaoSVXY_NP - Coefficients(45)*dt ! go right
            else
                !stay
            end if
        end if

        ! Towards constant TKE
        if (TKEconstOn) then
            TKEdiffLeft  = abs(TKE(4)-TKE(1))
            TKEdiffMid   = abs(TKE(2)-TKE(1))
            TKEdiffRight = abs(TKE(5)-TKE(1))
            if ((TKEdiffLeft  > TKEdiffMid) .and.
     .          (TKEdiffRight > TKEdiffMid)) then
                if (TKEdiffLeft > TKEdiffRight) then
                    SVaoSVXY_NP = SVaoSVXY_NP - Coefficients(45)*dt ! go right
                else
                    SVaoSVXY_NP = SVaoSVXY_NP + Coefficients(45)*dt ! go left
                end if
            end if
        end if

        ! Towards closer wall
        if (visualOn) then
            if      (wallDistances(4) > wallDistances(5)) then
                SVaoSVXY_NP = SVaoSVXY_NP - Coefficients(46)*dt ! go right
            else if (wallDistances(4) < wallDistances(5)) then
                SVaoSVXY_NP = SVaoSVXY_NP + Coefficients(46)*dt ! go left
            else
                ! stay
            end if
        end if

        ! Wall collision avoidance - sides: orient towards FV-180 (against FV)
        if ((SPFound_NP(4) == 0) .or.
     .      (SPFound_NP(5) == 0)) then
            SVaoSVXY_NP = ( FVaoSVXY_NP(1) - 180.0 )
        end if

        ! Random overlay
        call RandomFromSeed(RRR,SEED)
        if ((SPFound_NP(4) == 1) .and. (SPFound_NP(5) == 1)) then
            ! open water: full angle range
            rand1 = (RRR*2.0 - 1.0) ! rand1 = -1..+1
            SVaoSVXY_NP = SVaoSVXY_NP + rand1*Coefficients(44)*dt
        else
            ! tactile contact to wall: that side impossible
            rand1 = RRR             ! rand1 =  0..1
            if      ((SPFound_NP(4) == 0) .and. (SPFound_NP(5) /= 0))
     .      then
                SVaoSVXY_NP = SVaoSVXY_NP - rand1*Coefficients(44)*dt ! go right
            else if ((SPFound_NP(5) == 0) .and. (SPFound_NP(4) /= 0))
     .      then
                SVaoSVXY_NP = SVaoSVXY_NP + rand1*Coefficients(44)*dt ! go left
            else
                ! both sides confined, no random overlay
            end if
        end if

        ! Vertical angle
        SVaoCFDXYZ_NP(2) = -FVaoCFDXYZ_NP(2,1)

        ! anti-stuck
        call RandomFromSeed(RRR,SEED)
        if (t_same_spot(FishNumber) > Coefficients(57)) then
            SVaoSVXY_NP = RRR * 359.9
            ! unnecessary limit? Could probably be omitted in future work
            if (RRR > 0.5) then
                SVaoCFDXYZ_NP(2) = 0.0
            end if
        end if

!     Assign swim speed (as in still water):
        FhSpdRes_NP = CondVelM(1) + FHBODYLENGTH * Coefficients(52)

!     Check swim angle limits:
        ! 0-360 check, robust version
        do while (SVaoSVXY_NP < 0.0)
            SVaoSVXY_NP = SVaoSVXY_NP + 360.0
        end do
        do while (SVaoSVXY_NP >= 360.0)
            SVaoSVXY_NP = SVaoSVXY_NP - 360.0
        end do

        ! Make angles go from 0.0 to 180.0 and 0.0 to -180.0
        if (SVaoSVXY_NP > 180.0) then
          SVaoSVXY_NP = SVaoSVXY_NP - 360.0
        end if

        ! -90..+90 force check
        if (SVaoCFDXYZ_NP(2) >  90.0) SVaoCFDXYZ_NP(2) =  90.0
        if (SVaoCFDXYZ_NP(2) < -90.0) SVaoCFDXYZ_NP(2) = -90.0


!****************************************************************************************************
!****************************************************************************************************
!****************************************    SECTION 3.2    *****************************************
!****************************************************************************************************
!****************************************************************************************************

! Behavior implemented here: swim with flow vector (Holding)

      else if (AgentChosen(1) == 2) then

        ! same amount (2) of RRR iterations per behavior for repeatability
        call RandomFromSeed(RRR,SEED)
        call RandomFromSeed(RRR,SEED)

!     Assign swim angle:
        ! Horizontal, against flow vector, relative to swim vector:
        SVaoCFDXYTemp = FVaoCFDXYZ_NP(1,1) + 180.0

        ! 0-360 check, robust version
        do while (SVaoCFDXYTemp < 0.0 )
            SVaoCFDXYTemp = SVaoCFDXYTemp + 360.0
        end do
        do while (SVaoCFDXYTemp >= 360.0)
            SVaoCFDXYTemp = SVaoCFDXYTemp - 360.0
        end do

        SVaoSVXY_NP = SVaoCFDXYTemp - SVaoCFDXYZ_NPm1(1,FishNumber)

        ! Vertical, absolute:
        SVaoCFDXYZ_NP(2) = -FVaoCFDXYZ_NP(2,1)

!     Assign swim speed (as in still water):
        FhSpdRes_NP = CondVelM(1)

!     Check swim angle limits:
        ! 0-360 check, robust version
        do while (SVaoSVXY_NP < 0.0)
            SVaoSVXY_NP = SVaoSVXY_NP + 360.0
        end do
        do while (SVaoSVXY_NP >= 360.0)
            SVaoSVXY_NP = SVaoSVXY_NP - 360.0
        end do

        ! Make angles go from 0.0 to 180.0 and 0.0 to -180.0
        if (SVaoSVXY_NP > 180.0) then
          SVaoSVXY_NP = SVaoSVXY_NP - 360.0
        end if

        ! -90..+90 force check
        if (SVaoCFDXYZ_NP(2) >  90.0) SVaoCFDXYZ_NP(2) =  90.0
        if (SVaoCFDXYZ_NP(2) < -90.0) SVaoCFDXYZ_NP(2) = -90.0


!****************************************************************************************************
!****************************************************************************************************
!****************************************    SECTION 3.3    *****************************************
!****************************************************************************************************
!****************************************************************************************************

! Behavior implemented here: drift with flow vector (Drifting)

      else if (AgentChosen(1) == 3) then

!     Assign swim angle:
        ! Horizontal, against flow vector, relative to swim vector:
        SVaoCFDXYTemp = FVaoCFDXYZ_NP(1,1) + 180.0

        ! 0-360 check, robust version
        do while (SVaoCFDXYTemp < 0.0 )
            SVaoCFDXYTemp = SVaoCFDXYTemp + 360.0
        end do
        do while (SVaoCFDXYTemp >= 360.0)
            SVaoCFDXYTemp = SVaoCFDXYTemp - 360.0
        end do

        SVaoSVXY_NP = SVaoCFDXYTemp - SVaoCFDXYZ_NPm1(1,FishNumber)

        ! Towards random side, in (1-driftStraightProb)*100 percent
        call RandomFromSeed(RRR,SEED)
        driftStraightProb = Coefficients(41)
        if ( RRR<=driftStraightProb ) then
            ! drift straight
        else
            ! drift sideways, 50:50 for side
            if (RRR-driftStraightProb > (1.0-driftStraightProb)/2) then
                SVaoSVXY_NP = SVaoSVXY_NP + Coefficients(18)
            else
                SVaoSVXY_NP = SVaoSVXY_NP - Coefficients(18)
            end if
        end if

!     Assign swim speed:
        call RandomFromSeed(RRR,SEED)
        FhSpdRes_NP = CondVelM(1) * RRR

!     Check swim angle limits:
        ! 0-360 check, robust version
        do while (SVaoSVXY_NP < 0.0)
            SVaoSVXY_NP = SVaoSVXY_NP + 360.0
        end do
        do while (SVaoSVXY_NP >= 360.0)
            SVaoSVXY_NP = SVaoSVXY_NP - 360.0
        end do

        ! Make angles go from 0.0 to 180.0 and 0.0 to -180.0
        if (SVaoSVXY_NP > 180.0) then
          SVaoSVXY_NP = SVaoSVXY_NP - 360.0
        end if

!       ! -90..+90 force check not required


!****************************************************************************************************
!****************************************************************************************************
!****************************************    SECTION 3.4    *****************************************
!****************************************************************************************************
!********************************* Do Not Modify Between Lines **************************************

      else

        WRITE(*,*) 'Error: AgentChosen(1) calculation:'
        WRITE(*,*) '  AgentChosen(1) = ', AgentChosen(1)
        STOP

      end if   ! (AgentChosen(1) == ?)


!****************************************************************************************************
!****************************************************************************************************
!****************************************    SECTION 5.1    *****************************************
!****************************************************************************************************
!****************************************************************************************************

!Vertical Override Agent Behavior.

! Behavior implemented here: swim toward acclimatized elevation/pressure:
!  Involves implementing the behavior selected for xy-plane movement, but
!  vertical component is overridden by swimming toward acclimatized elevation/pressure:

      if (AgentChosen(2) == 4) then                                     ! Pressure (swim toward acclimatized elevation)

        ! Assign swim angle:
        SVaoCFDXYZ_NP(2) = -Coefficients(48)*dt

      elseif (AgentChosen(2) == 5) then

        ! Assign swim angle:
        SVaoCFDXYZ_NP(2) = Coefficients(48)*dt

      end if

! Check swim angle limits:
      if (SVaoCFDXYZ_NP(2) >  90.0) SVaoCFDXYZ_NP(2) =  90.0
      if (SVaoCFDXYZ_NP(2) < -90.0) SVaoCFDXYZ_NP(2) = -90.0

! Prevent swimming through bottom if close to it
      if ( (SPFound_NP(7) == 0) .and.
     .     (SVaoCFDXYZ_NP(2) < 5.0 ) ) then
        SVaoCFDXYZ_NP(2) = 5.0
      end if


!****************************************************************************************************
!****************************************************************************************************
!****************************************    SECTION 6.1    *****************************************
!****************************************************************************************************
! The following portions are provided courtesy of the United States Army
!   Corps of Engineers, Engineer Research and Development Center,
!   Environmental Laboratory. See license file for information on usage.
!****************************** Do Not Modify Code Past This Line ***********************************

      IF (ExtraDiagnostics .AND. WriteToDebugFile) THEN
        WRITE(72,*) ' '
        WRITE(72,*) '    AgtDecision(1,FishNumber,TimeStep)      = ',
     .                   AgtDecision(1,FishNumber,TimeStep)
        WRITE(72,*) '    AgtDecision(2,FishNumber,TimeStep)      = ',
     .                   AgtDecision(2,FishNumber,TimeStep)
        WRITE(72,*) '    AgtDecision(3,FishNumber,TimeStep)      = ',
     .                   AgtDecision(3,FishNumber,TimeStep)
        WRITE(72,*) '    AgtDecision(4,FishNumber,TimeStep)      = ',
     .                   AgtDecision(4,FishNumber,TimeStep)
        WRITE(72,*) ' '
        WRITE(72,*) '    SVaoSVXY_NP      After  SP        Check = ',
     .                   SVaoSVXY_NP
        WRITE(72,*) '    SVaoCFDXYZ_NP(2) After  SP & Agt4 Check = ',
     .                   SVaoCFDXYZ_NP(2)
        WRITE(72,*) '      FVaoCFDXYZ_NP(2,1) = ',FVaoCFDXYZ_NP(2,1)
        WRITE(72,*) ' '
        WRITE(72,*) '    FhSpdRes_NP          = ',FhSpdRes_NP
        WRITE(72,*) '      FhBodyLngthVel     = ',
     .                     FhSpdRes_NP / FHBODYLENGTH
        WRITE(72,*) '      FhAttrb_NP         = ',
     .                     FhAttrb_NP(1,FishNumber)
      END IF


  230 format (F7.1, I6, 2(F6.2), I3)

      if (ExtraDiagnostics .and. WriteToDebugFile) then
        write(110,230) TimeStep*dt, FishNumber,
     .  AgtUtil_NP(1), AgtCostAvg, AgentChosen(1)
      end if


!End of Behavior Calculations.
!  Purpose of the behavior rule code above is for you to have at this point the following 3 variables:
!    SVaoSVXY_NP      = the swim vector angle off of (relative to) the previous swim vector in the
!                       horizontal plane. This value must range from -180.0 to 180.0 and will be
!                       used in the "Do Not Modify" section below to calculate SVaoCFDXYZ_NP(1), which is
!                       the swim vector angle off of (relative to) the CFD model mesh xy-plane axes.
!    SVaoCFDXYZ_NP(2) = the vertical swim vector angle off of (relative to) the CFD model mesh xy-plane; in other
!                       words, the vertical angle off horizontal. This value must range from -90.0 to 90.0.
!    FhSpdRes_NP      = fish velocity (resultant value).


!****************************** Do Not Modify Code Past This Line ***********************************
!****************************************************************************************************
!****************************************    SECTION 6.2    *****************************************
!****************************************************************************************************
!****************************** Do Not Modify Code Past This Line ***********************************

!Check SVaoSVXY_NP and SVaoCFDXYZ_NP(2) Values.

      IF ((SVaoSVXY_NP < -180.0) .OR.
     .    (SVaoSVXY_NP >  180.0)) THEN
        WRITE(*,*) 'Error: SVaoSVXY_NP calculation'
        WRITE(*,*) '  SVaoSVXY_NP = ',SVaoSVXY_NP
        STOP
      END IF
      IF ((SVaoCFDXYZ_NP(2) < -90.0) .OR.
     .    (SVaoCFDXYZ_NP(2) >  90.0)) THEN
        WRITE(*,*) 'Error: SVaoCFDXYZ_NP(2) calculation'
        WRITE(*,*) '  SVaoCFDXYZ_NP(2) = ',SVaoCFDXYZ_NP(2)
        STOP
      END IF


!****************************** Do Not Modify Code Past This Line ***********************************
!****************************************************************************************************
!****************************************    SECTION 6.3    *****************************************
!****************************************************************************************************
!****************************** Do Not Modify Code Past This Line ***********************************

!Convert Swim Angle/Velocity to Vectors.

      SVaoCFDZRad = SVaoCFDXYZ_NP(2) * (3.14159 / 180.0)
      FhSpdResXY  = FhSpdRes_NP * COS(SVaoCFDZRad)
      IF      ((SVaoSVXY_NP >=  0.0) .AND.
     .         (SVaoSVXY_NP <= 90.0)) THEN
        SVaoSVXYRad     =          SVaoSVXY_NP  * (3.14159 / 180.0)
        SVvoSVXYZ_NP(1) =  FhSpdResXY * COS(SVaoSVXYRad)                 ! (S)wim (V)ector (v)elocity (o)ff the previous (S)wim (V)ector's xy-axis
        SVvoSVXYZ_NP(2) =  FhSpdResXY * SIN(SVaoSVXYRad)                 ! (S)wim (V)ector (v)elocity (o)ff the previous (S)wim (V)ector's xy-axis
      ELSE IF ((SVaoSVXY_NP >=  90.0) .AND.
     .         (SVaoSVXY_NP <= 180.0)) THEN
        SVaoSVXYRad     = (180.0 - SVaoSVXY_NP) * (3.14159 / 180.0)
        SVvoSVXYZ_NP(1) = -FhSpdResXY * COS(SVaoSVXYRad)
        SVvoSVXYZ_NP(2) =  FhSpdResXY * SIN(SVaoSVXYRad)
      ELSE IF ((SVaoSVXY_NP <=   0.0) .AND.
     .         (SVaoSVXY_NP >= -90.0)) THEN
        SVaoSVXYRad     =         -SVaoSVXY_NP  * (3.14159 / 180.0)
        SVvoSVXYZ_NP(1) =  FhSpdResXY * COS(SVaoSVXYRad)
        SVvoSVXYZ_NP(2) = -FhSpdResXY * SIN(SVaoSVXYRad)
      ELSE IF ((SVaoSVXY_NP <=  -90.0) .AND.
     .         (SVaoSVXY_NP >= -180.0)) THEN
        SVaoSVXYRad     = (180.0 + SVaoSVXY_NP) * (3.14159 / 180.0)
        SVvoSVXYZ_NP(1) = -FhSpdResXY * COS(SVaoSVXYRad)
        SVvoSVXYZ_NP(2) = -FhSpdResXY * SIN(SVaoSVXYRad)
      ELSE
        WRITE(*,*) 'Error: Invalid SVaoSVXY_NP value'
        STOP
      END IF

! Determine angle of the swim vector relative to the
! CFD xy-axis:

      SVaoCFDXYZ_NP(1) = SVaoCFDXYZ_NPm1(1,FishNumber) + SVaoSVXY_NP     ! (S)wim (V)ector (a)ngle (o)ff CFD xy-axis orientation
      IF      (SVaoCFDXYZ_NP(1) < 0.0) THEN
        SVaoCFDXYZ_NP(1) = SVaoCFDXYZ_NP(1) + 360.0
      ELSE IF (SVaoCFDXYZ_NP(1) > 360.0) THEN
        SVaoCFDXYZ_NP(1) = SVaoCFDXYZ_NP(1) - 360.0
      END IF
      SVaoCFDXYTemp = SVaoCFDXYZ_NP(1)
      IF (SVaoCFDXYTemp > 180.0) THEN                                    ! Ensure that angles go from 0.0 to 180.0 and 0.0 to -180.0
        SVaoCFDXYTemp = SVaoCFDXYTemp - 360.0
      END IF
      IF      ((SVaoCFDXYTemp >=  0.0) .AND.
     .         (SVaoCFDXYTemp <= 90.0)) THEN
        SVaoCFDXYRad     =          SVaoCFDXYTemp  * (3.14159 / 180.0)
        SVvoCFDXYZ_NP(1) =  FhSpdResXY * COS(SVaoCFDXYRad)               ! (S)wim (V)ector (v)elocity (o)ff CFD xy-axis orientation
        SVvoCFDXYZ_NP(2) =  FhSpdResXY * SIN(SVaoCFDXYRad)               ! (S)wim (V)ector (v)elocity (o)ff CFD xy-axis orientation
      ELSE IF ((SVaoCFDXYTemp >=  90.0) .AND.
     .         (SVaoCFDXYTemp <= 180.0)) THEN
        SVaoCFDXYRad     = (180.0 - SVaoCFDXYTemp) * (3.14159 / 180.0)
        SVvoCFDXYZ_NP(1) = -FhSpdResXY * COS(SVaoCFDXYRad)
        SVvoCFDXYZ_NP(2) =  FhSpdResXY * SIN(SVaoCFDXYRad)
      ELSE IF ((SVaoCFDXYTemp <=   0.0) .AND.
     .         (SVaoCFDXYTemp >= -90.0)) THEN
        SVaoCFDXYRad     =         -SVaoCFDXYTemp  * (3.14159 / 180.0)
        SVvoCFDXYZ_NP(1) =  FhSpdResXY * COS(SVaoCFDXYRad)
        SVvoCFDXYZ_NP(2) = -FhSpdResXY * SIN(SVaoCFDXYRad)
      ELSE IF ((SVaoCFDXYTemp <=  -90.0) .AND.
     .         (SVaoCFDXYTemp >= -180.0)) THEN
        SVaoCFDXYRad     = (180.0 + SVaoCFDXYTemp) * (3.14159 / 180.0)
        SVvoCFDXYZ_NP(1) = -FhSpdResXY * COS(SVaoCFDXYRad)
        SVvoCFDXYZ_NP(2) = -FhSpdResXY * SIN(SVaoCFDXYRad)
      ELSE
        WRITE(*,*) 'Error: Invalid SVaoCFDXYZ_NP(1) value'
        STOP
      END IF
      SVvoCFDXYZ_NP(3) = FhSpdRes_NP * SIN(SVaoCFDZRad)
      IF (ExtraDiagnostics .AND. WriteToDebugFile) THEN
        WRITE(72,*) ' '
        WRITE(72,*) '    SVaoCFDXYZ_NP(1)     = ',SVaoCFDXYZ_NP(1)
        WRITE(72,*) '      SVaoCFDXYZ_NP_180  = ',SVaoCFDXYTemp
        WRITE(72,*) '      FVaoCFDXYZ_NP(1,1) = ',FVaoCFDXYZ_NP(1,1)
      END IF
      IF ((SVvoCFDXYZ_NP(1) == 0.0) .AND.
     .    (SVvoCFDXYZ_NP(2) == 0.0)) THEN
        VldSVOrientXYZ_NP(1) = 0                                         ! No valid xy-plane orientation of fish axis => no fish movement
      ELSE
        VldSVOrientXYZ_NP(1) = 1                                         ! Valid xy-plane orientation of fish axis can be determined
      END IF

! Determine angle of the swim vector relative to the
! flow vector xy-axis:

      IF ((VldSVOrientXYZ_NP(1) == 1) .AND.
     .    (VldFVOrientXYZ_NP(1) == 1)) THEN
        SVaoFVXY_NP = SVaoCFDXYZ_NP(1) - FVaoCFDXYZ_NP(1,1)              ! (S)wim (V)ector (a)ngle (o)ff (F)low (V)ector's xy-axis
        IF      (SVaoFVXY_NP < 0.0) THEN
          SVaoFVXY_NP = SVaoFVXY_NP + 360.0
        ELSE IF (SVaoFVXY_NP > 360.0) THEN
          SVaoFVXY_NP = SVaoFVXY_NP - 360.0
        END IF
        IF (SVaoFVXY_NP > 180.0) SVaoFVXY_NP = SVaoFVXY_NP - 360.0       ! Make it so that angles go from 0.0 to 180.0 and 0.0 to -180.0
      ELSE
        SVaoFVXY_NP = 0.0                                                ! (S)wim (V)ector (a)ngle (o)ff (F)low (V)ector's xy-axis
      END IF

! Determine velocity of the swim vector relative to the
! flow vector xyz-axis velocity:

      IF      ((SVaoFVXY_NP >=  0.0) .AND.
     .         (SVaoFVXY_NP <= 90.0)) THEN
        SVaoFVXYRad     =          SVaoFVXY_NP  * (3.14159 / 180.0)
        SVvoFVXYZ_NP(1) =  FhSpdResXY * COS(SVaoFVXYRad)                 ! (S)wim (V)ector (v)elocity (o)ff (F)low (V)ector's xy-axis
        SVvoFVXYZ_NP(2) =  FhSpdResXY * SIN(SVaoFVXYRad)                 ! (S)wim (V)ector (v)elocity (o)ff (F)low (V)ector's xy-axis
      ELSE IF ((SVaoFVXY_NP >=  90.0) .AND.
     .         (SVaoFVXY_NP <= 180.0)) THEN
        SVaoFVXYRad     = (180.0 - SVaoFVXY_NP) * (3.14159 / 180.0)
        SVvoFVXYZ_NP(1) = -FhSpdResXY * COS(SVaoFVXYRad)
        SVvoFVXYZ_NP(2) =  FhSpdResXY * SIN(SVaoFVXYRad)
      ELSE IF ((SVaoFVXY_NP <=   0.0) .AND.
     .         (SVaoFVXY_NP >= -90.0)) THEN
        SVaoFVXYRad     =         -SVaoFVXY_NP  * (3.14159 / 180.0)
        SVvoFVXYZ_NP(1) =  FhSpdResXY * COS(SVaoFVXYRad)
        SVvoFVXYZ_NP(2) = -FhSpdResXY * SIN(SVaoFVXYRad)
      ELSE IF ((SVaoFVXY_NP <=  -90.0) .AND.
     .         (SVaoFVXY_NP >= -180.0)) THEN
        SVaoFVXYRad     = (180.0 + SVaoFVXY_NP) * (3.14159 / 180.0)
        SVvoFVXYZ_NP(1) = -FhSpdResXY * COS(SVaoFVXYRad)
        SVvoFVXYZ_NP(2) = -FhSpdResXY * SIN(SVaoFVXYRad)
      ELSE
        WRITE(*,*) 'Error: Invalid SVaoFVXY_NP value'
        STOP
      END IF
      SVvoFVXYZ_NP(3) = SVvoCFDXYZ_NP(3) - FishSensoryVelocity_NP(3,1)


      IF (ExtraDiagnostics .AND. WriteToDebugFile) THEN
        WRITE(72,*) ' '
        write(72,*) '    acclMagLeft  = ',CondAcclM(4)
        write(72,*) '    acclMagRight = ',CondAcclM(5)
        WRITE(72,*) ' '
        write(72,*) '    TKECenter    = ',TKE(1)
        write(72,*) '    TKEFront     = ',TKE(2)
        write(72,*) '    TKELeft      = ',TKE(4)
        write(72,*) '    TKERight     = ',TKE(5)
        write(72,*) ' '
        WRITE(72,*) '    VelMfish = ',
     .                   SQRT( SVvoCFDXYZ_NP(1) *
     .                         SVvoCFDXYZ_NP(1) +
     .                         SVvoCFDXYZ_NP(2) *
     .                         SVvoCFDXYZ_NP(2) +
     .                         SVvoCFDXYZ_NP(3) *
     .                         SVvoCFDXYZ_NP(3) )
        WRITE(72,*) '      Ufish = ',SVvoCFDXYZ_NP(1)
        WRITE(72,*) '      Vfish = ',SVvoCFDXYZ_NP(2)
        WRITE(72,*) '      Wfish = ',SVvoCFDXYZ_NP(3)
        WRITE(72,*) ' '
        WRITE(72,*) '    VelM = ',CondVelM(1)
        WRITE(72,*) '      U = ',FishSensoryVelocity_NP(1,1)
        WRITE(72,*) '      V = ',FishSensoryVelocity_NP(2,1)
        WRITE(72,*) '      W = ',FishSensoryVelocity_NP(3,1)
        WRITE(72,*) ' '
        WRITE(72,*) '+---------------------------------------------------
     .--------------------------+'
        WRITE(72,*) ' '
      END IF


! Account for different array numbering - to be set at file head
      FishNumber = FishNumber-1
      TimeStep = TimeStep-1


      return
      end
