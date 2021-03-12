!-----------------------------------------------------------------------
! Application
!   writeOutput
!
! Purpose
!   Output the results to files, for Tecplot, ParaView and Matlab
!
! Portions of writeOutput are
!   Copyright (c) 2021 David Gisen, Bundesanstalt fuer Wasserbau
!
! Further portions of BehaviorRule are
!   provided courtesy of the United States Army Corps of Engineers,
!   Engineer Research and Development Center, Environmental Laboratory
!
! Portions are indicated in the code by respective headings. See license
!   file for information on usage.
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! The following portions are
! Copyright (c) 2021 David Gisen, Bundesanstalt fuer Wasserbau
!   See license file for information on usage
!-----------------------------------------------------------------------

      Subroutine openOutputFiles
     .(
     I  writeToDebugFile,
     I  longResults
     .)

! Open Fish Data Files. CLOSE'd automatically after program end

      implicit none

      integer*4 :: ierror
      integer   :: longResults
      logical*1 :: writeToDebugFile ! "*" notation not in Fortran standard - provides transfer

      ierror = 0

      if (longResults == 1) then
        OPEN(UNIT=20001,FILE="output/v_TecTrack_ZonesAreTime.dat",
     .       STATUS='REPLACE',IOSTAT=ierror)                             ! Opening a file to output simulated fish tracks for Tecplot

        ! rely on user to create subfolder
        if (ierror /= 0) then
          stop 'ABORT: Could not find folder "output".'
        end if

        WRITE(20001,9997)                                                ! Creating an output file header for Tecplot
 9997   FORMAT('TITLE = "Virtual Fish Tracks - Organized by Time"')
        WRITE(20001,9972)
 9972   FORMAT('VARIABLES = "x_fish", "y_fish", "z_fish", "WSC", "Time",
     . "Node BC", "uFishCFD", "vFishCFD", "wFishCFD",
     . "AgtUtil1",
     . "FID", "FhSpdRes",
     . "FhAttrb1", "u_flow", "v_flow", "w_flow"')
      end if

      OPEN(UNIT=73,FILE="output/v_TecPassageResults.dat",
     .  STATUS='REPLACE')
      WRITE(73,*) 'TITLE = "Passage Results (%)"'
      WRITE(73,*) 'VARIABLES = "Route #", "Percent of Fish"'

      if (writeToDebugFile) then
        OPEN(UNIT=110,FILE="output/internalStatesDebug.txt",
     .    STATUS='REPLACE')

        WRITE(110,*) '"Time" "FishNumber" "Motivation" "Fatigue"
     .    "Behavior"'

        OPEN(UNIT=72,FILE="output/v_Debug.txt",STATUS='REPLACE')
      end if

      return
      end


!***********************************************************************


      Subroutine writeSensoryLocations
     .(
     I  FN,
     I  TimeStep,
     I  FSOlimit,
     I  OVOIDLENGTH,
     I  SenPointXYZ
     .)

      implicit none

      logical*1 ::  extraDiagnostics

      integer ::    FN, TimeStep, FSOlimit

      real*8 ::     OVOIDLENGTH,
     .              SenPointXYZ(3,FSOlimit)

      ! Account for different array indexing - to be reversed at end
      FN = FN+1
      TimeStep = TimeStep+1;

              WRITE(72,*) 'TimeStep = ',TimeStep
              WRITE(72,*) 'Fish #   = ',FN
              write(72,*) ' '
              WRITE(72,*) '    OVOID avg. radius = ',OVOIDLENGTH
              WRITE(72,*) ' '

      ! Account for different array numbering - to be set at begin
      FN = FN-1
      TimeStep = TimeStep-1;

      return
      end


!***********************************************************************

!-----------------------------------------------------------------------
! The following portions are provided courtesy of the United States Army
!   Corps of Engineers, Engineer Research and Development Center,
!   Environmental Laboratory. See license file for information on usage.
!-----------------------------------------------------------------------


      SUBROUTINE OutputFishData_ZonesAreTime
     .(
     I  FN,
     I  nFish,
     I  TimeStep,
     I  nTimeSteps,
     I  DT,
     I  FSOlimit,
     I  nAgents,
     I  FishLocation,
     I  VelFishCFD,
     I  AgtUtil,
     I  FhSpdRes,
     I  FhAttrb,
     I  totalExitTally,
     I  fishSensoryVelocity
     .)

      IMPLICIT NONE


      integer ::    FN,
     .              nFish,
     .              TimeStep,
     .              nTimeSteps,
     .              NPorTS,
     .              FhAttrb(1,nFish,nTimeSteps),
     .              totalExitTally,
     .              FSOlimit,
     .              nAgents

      real*8 ::     DT,
     .              FishLocation(3,nFish,nTimeSteps),
     .              VelFishCFD(3,nFish,nTimeSteps),
     .              AgtUtil(3,nFish,nTimeSteps),
     .              FhSpdRes(nFish,nTimeSteps),
     .              fishSensoryVelocity(3,FSOlimit,nFish)

      ! Account for different array indexing - to be reversed at end
      FN = FN+1
      TimeStep = TimeStep+1

      NPorTS = TimeStep

 9900 format(A,F8.3,A,I6,A,F8.3,A)

      if (FN == 1) then
        IF (TimeStep > 2) WRITE(20001,*) ' '
        WRITE(20001,9900) 'Zone T="Time',REAL(TimeStep)*DT,
     .      '", I=',nFish,       !-totalExitTally, ! only subtract if line output reduces with leaves
     .      ', SOLUTIONTIME=',REAL(TimeStep)*DT,
     .      ', STRANDID=1, F=POINT'
      end if

!Output Fish Data for Each Virtual Fish.

      WRITE(20001,9973) FishLocation(1,FN,NPorTS),                       ! X-position (m) of virtual fish at time step
     .                  FishLocation(2,FN,NPorTS),                       ! Y-position (m) of virtual fish at time step
     .                  FishLocation(3,FN,NPorTS),                       ! Z-position (m) of virtual fish at time step
     .                  0,                                               ! Must be zero so virtual fish is not "value-blanked" in Tecplot
     .                  REAL(TimeStep)*DT,                               ! Time (sec)
     .                  0,                                               ! Must be zero so virtual fish is not "value-blanked" in Tecplot
     .                  VelFishCFD(1,FN,NPorTS),                         ! Vel of fish in the CFD x-direction as seen by CFD ref frame
     .                  VelFishCFD(2,FN,NPorTS),                         ! Vel of fish in the CFD y-direction as seen by CFD ref frame
     .                  VelFishCFD(3,FN,NPorTS),                         ! Vel of fish in the CFD z-direction as seen by CFD ref frame
     .                  AgtUtil(1,FN,NPorTS),                            ! NPorTS = TimeStep when this is output
     .                  FN,                                              ! Fish ID#
     .                  FhSpdRes(FN,NPorTS),                             ! Fish speed (m/s)
     .                  FhAttrb(1,FN,NPorTS),
     .                  fishSensoryVelocity(1,1,FN),                     ! u_flow
     .                  fishSensoryVelocity(2,1,FN),                     ! v_flow
     .                  fishSensoryVelocity(3,1,FN)                      ! w_flow
 9973   FORMAT(3(E21.12),I6,E18.8,I6,3(E18.6),1(F6.1),1(I6),1(E18.6),
     .        1(I6),3(E18.6))

      ! Account for different array numbering - to be set at begin
      FN = FN-1
      TimeStep = TimeStep-1

      RETURN
      END SUBROUTINE OutputFishData_ZonesAreTime


!***********************************************************************


      SUBROUTINE OutputFishPassageAndDecisions
     .(
     I  NFISH,
     I  TimeStep,
     I  DT,
     I  NUMEXITROUTES,
     I  numBoundaries,
     I  EXITTALLY
     .)

      IMPLICIT NONE

      integer ::    NFISH,
     .              TimeStep,
     .              NUMEXITROUTES,
     .              numBoundaries,
     .              EXTS,
     .              EXITTALLY(2,numBoundaries)
!
      real*8 ::     DT

      ! Account for different array indexing - to be reversed at end
      TimeStep = TimeStep+1

! Write-out the number of virtual fish exiting each passage route:

 9950 format (A, F9.1, A, I8, A, F9.1, A)

            IF (TimeStep > 2) WRITE(73,*) ' '

            WRITE(73,9950) 'Zone T="Time',REAL(TimeStep)*DT,
     .                  '", I=',NUMEXITROUTES,
     .                  ', SOLUTIONTIME=',REAL(TimeStep)*DT,
     .                  ', STRANDID=1, F=POINT'                          ! I = Number of bar charts to be displayed

          DO EXTS=1, numBoundaries
            if (EXITTALLY(2,EXTS) /= 0) then
              WRITE(73,9999) EXTS, !POS1,
     .                       REAL(EXITTALLY(1,EXTS))/REAL(NFISH)*100.0
            end if

          END DO

 9999     FORMAT (I6, '     ', F5.1)

      ! Account for different array numbering - to be set at begin
      TimeStep = TimeStep-1

      RETURN
      END SUBROUTINE OutputFishPassageAndDecisions


!***********************************************************************


      SUBROUTINE ENDOFRUNOUTPUT
     .(
     I      NFISH,
     I      TimeStep,
     I      nTimeSteps,
     I      DT,
     I      NumDecisions,
     I      OUTINTV,
     I      FishLocation,
     I      VelFishCFD,
     I      FhSpdRes,
     I      FhAttrb,
     I      SVAOCFDXYZ,
     I      longResults,
     I      statesStatistics,
     I      NumDecisionsSum
     .)

      IMPLICIT NONE

      integer ::    TimeStep,
     .              nTimeSteps,
     .              NFISH,
     .              NumDecisions(NFISH),
     .              FN6,
     .              ITBOUT,
     .              OUTINTV,
     .              NPorTS,
     .              FhAttrb(1,nFish,nTimeSteps),
     .              longResults,
     .              NumDecisionsSum,
     .              row


      real*8 ::     DT,
     .              FishLocation(3,NFISH,nTimeSteps+1),
     .              VelFishCFD(3,NFISH,nTimeSteps),
     .              FhSpdRes(nFish,nTimeSteps),
     .              SVAOCFDXYZ(2,NFISH,nTimeSteps),
     .              statesStatistics(8,NFISH),
     .              outputMatrix(4,NumDecisionsSum)


      if (longResults /= 1) then
        goto 7301
      endif

!Write ASCII-Formatted Fish Data (only if activated).

      ! Opening a file to output simulated fish tracks for Tecplot format
      OPEN(UNIT=40000,FILE=
     .      "output/v_TecTrack_ZonesAreIndividualFish.dat",
     .      STATUS='REPLACE')
      WRITE(40000,9977)                                                  ! Creating an output file header for Tecplot
 9977 FORMAT('TITLE = "Virtual Fish Tracks - Organized by Individual
     .Fish"')
      WRITE(40000,9998)
 9998 FORMAT('VARIABLES = "x_fish", "y_fish", "z_fish", "WSC", "Time",
     . "Node BC", "uFishCFD", "vFishCFD", "wFishCFD",
     . "FID", "FhSpdRes",
     . "SVaoCFD_XY", "SVaoCFD_Z",
     . "FhAttrb1"')

 9980 format (A, I7, A, I10, A)

      ! Write
      DO FN6=1,NFISH
        IF (FN6 > 1) WRITE(40000,*) ' '

        WRITE(40000,9980) 'Zone T="Fish ID#',FN6,'", I=',
     .                  NumDecisions(FN6),', F=POINT'
        ITBOUT = 0                                                       ! ITBOUT = Number of iterations since last output

        DO TimeStep=2,NumDecisions(FN6)+1
          ITBOUT = ITBOUT + 1
          IF (TimeStep == 2) ITBOUT = OUTINTV                            ! Ensures fish data is output for the first time step for which there's movement
          IF (ITBOUT /= OUTINTV) GOTO 9953
          NPorTS = TimeStep
          WRITE(40000,9996) FishLocation(1,FN6,NPorTS),                  ! X-position (m) of virtual fish at time step
     .                      FishLocation(2,FN6,NPorTS),                  ! Y-position (m) of virtual fish at time step
     .                      FishLocation(3,FN6,NPorTS),                  ! Z-position (m) of virtual fish at time step
     .                      0,                                           ! Must be zero so virtual fish is not "value-blanked" in Tecplot
     .                      REAL(TimeStep)*DT,                           ! Time (sec)
     .                      0,                                           ! Must be zero so virtual fish is not "value-blanked" in Tecplot
     .                      VelFishCFD(1,FN6,NPorTS),                    ! Vel of fish in the CFD x-direction as seen by CFD ref frame
     .                      VelFishCFD(2,FN6,NPorTS),                    ! Vel of fish in the CFD y-direction as seen by CFD ref frame
     .                      VelFishCFD(3,FN6,NPorTS),                    ! Vel of fish in the CFD z-direction as seen by CFD ref frame
     .                      FN6,                                         ! Fish ID#
     .                      FhSpdRes(FN6,NPorTS),                        ! Fish speed (m/sec)
     .                      SVAOCFDXYZ(1,FN6,NPorTS),                    ! (S)wim (V)ector (A)ngle (O)ff CFD X-Y axis orientation
     .                      SVAOCFDXYZ(2,FN6,NPorTS),                    ! (S)wim (V)ector (A)ngle (O)ff CFD xy-plane
     .                      FhAttrb(1,FN6,NPorTS)

 9953     CONTINUE
          IF (ITBOUT == OUTINTV) ITBOUT = 0                              ! Set counter ITBOUT = 0
        END DO

      END DO

 9996 FORMAT(3(E21.12),I12,E25.12,I12,3(E18.6),1(I6),3(E18.6),1(I6))

      CLOSE(40000)

 7301 continue


!-----------------------------------------------------------------------
! The following portions are
! Copyright (c) 2021 David Gisen, Bundesanstalt fuer Wasserbau
!   See license file for information on usage
!-----------------------------------------------------------------------


!Write Binary-Formatted Fish Data (always).

      ! Gather data
      row = 1
      do FN6=1,nFish
        do TimeStep=2,NumDecisions(FN6)+1
            outputMatrix(1,row) = real(TimeStep)*dt
            outputMatrix(2:4,row) = FishLocation(1:3,FN6,TimeStep)
            row = row+1
        end do
      end do

      ! Opening a file to output simulated fish tracks for Matlab format (binary)
      open(unit=40020,file=
     .      "output/time_xyz_byFish_binary.dat", action="write",
     .      status='replace', form="unformatted", access="stream")

      write(40020) outputMatrix

      close(unit=40020)



!Write internal states

      open(unit=120,file="output/internalStates.txt",
     .  status='replace')
      write(120,1201)
 1201 format('"FishNumber" "Track duration" "Min. Motivation"
     .  "Max. Motivation" "Mean Motivation" "Min. Fatigue"
     .  "Max. Fatigue" "Mean Fatigue"')

      do FN6=1,nFish
        write(120,1210) int(statesStatistics(1,FN6)),
     .  statesStatistics(2,FN6),
     .  statesStatistics(3,FN6),
     .  statesStatistics(4,FN6),
     .  statesStatistics(5,FN6),
     .  statesStatistics(6,FN6),
     .  statesStatistics(7,FN6),
     .  statesStatistics(8,FN6)
      end do
 1210 format (I8, F18.4, 6(F12.6))
      close(120)


      RETURN
      END SUBROUTINE ENDOFRUNOUTPUT
