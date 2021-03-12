!/*--------------------------------------------------------------------*\
!Application
!   VectorRelation
!
!Purpose
!   Transform vectors relative to CFD model coordinates
!
! This file is provided courtesy of the United States Army Corps of
!   Engineers, Engineer Research and Development Center, Environmental
!   Laboratory. See license file for information on usage.
!
!\*--------------------------------------------------------------------*/

                                                                                                  ! Equivalence in BehaviorRule code:
      Subroutine VectorRelation(
     I                          VectorType,                              ! Local scalar             => PositionVector=0 ; VelocityVector=1
     I                          MVComponentsCFD,                         ! Local array              MVComponentsCFD(1:3)
     I                          VldRVOrientXYZ,                          ! Local array              VldRVOrientXYZ(1:2)
     O                          VldMVOrientXYZ,                          ! Local array              VldMVOrientXYZ(1:2)
     I                          RVaoCFDXYZ,                              ! Local array              RVaoCFDXYZ(1:2)
     I                          RVvoCFDXYZ,                              ! Local array              RVvoCFDXYZ(1:3)
     O                          MVaoRVXY,                                ! Local scalar             => (M)aster (V)ector (a)ngle (o)ff (R)eference (V)ector's xy-axis
     O                          MVaoCFDXYZ,                              ! Local array              MVaoCFDXYZ(1:2)
     O                          MVvoRVXYZ                                ! Local array              MVvoRVXYZ(1:3)       => (M)aster (V)ector (v)elocity (o)ff (R)eference (V)ector's xy-axis
     .                          )
      IMPLICIT NONE


      INTEGER VldMVOrientXYZ(2),                                         ! Arrays
     .        VldRVOrientXYZ(2),
     .        VectorType,XYZ                                             ! Scalars

      real*8  MVComponentsCFD(3),                                        ! Arrays
     .        RVaoCFDXYZ(2),
     .        MVaoCFDXYZ(2),
     .        RVvoCFDXYZ(3),
     .        MVvoRVCFDXYZ(3),
     .        MVvoRVXYZ(3),
     .        AbsMVaoCFDXY,                                              ! Scalars
     .        MVvroCFDXY,AbsMVaoCFDZ,MVaoRVXY,MVvroRVXY,MVaoRVXYRad


! Determine angle of the Master Vector relative to the
! CFD xy-axis:

      IF ((MVComponentsCFD(1) .EQ. 0.0) .AND.
     .    (MVComponentsCFD(2) .EQ. 0.0)) THEN
        VldMVOrientXYZ(1) = 0                                            ! No valid xy-plane orientation of master vector
        MVaoCFDXYZ(1)     = 0.0
      ELSE
        VldMVOrientXYZ(1) = 1                                            ! Valid xy-plane orientation of master vector can be determined
        AbsMVaoCFDXY = ABS(ATAN2(ABS(MVComponentsCFD(2)),
     .                           ABS(MVComponentsCFD(1))))

        IF      ((MVComponentsCFD(1) >= 0.0) .AND.
     .           (MVComponentsCFD(2) >= 0.0)) THEN
          MVaoCFDXYZ(1) =         AbsMVaoCFDXY * (180.0 / 3.14159)       ! (M)aster (V)ector (a)ngle (o)ff CFD xy-axis orientation
        ELSE IF ((MVComponentsCFD(1) <= 0.0) .AND.
     .           (MVComponentsCFD(2) >= 0.0)) THEN
          MVaoCFDXYZ(1) = 180.0 - AbsMVaoCFDXY * (180.0 / 3.14159)       ! (M)aster (V)ector (a)ngle (o)ff CFD xy-axis orientation
        ELSE IF ((MVComponentsCFD(1) <= 0.0) .AND.
     .           (MVComponentsCFD(2) <= 0.0)) THEN
          MVaoCFDXYZ(1) = 180.0 + AbsMVaoCFDXY * (180.0 / 3.14159)       ! (M)aster (V)ector (a)ngle (o)ff CFD xy-axis orientation
        ELSE IF ((MVComponentsCFD(1) >= 0.0) .AND.
     .           (MVComponentsCFD(2) <= 0.0)) THEN
          MVaoCFDXYZ(1) = 360.0 - AbsMVaoCFDXY * (180.0 / 3.14159)       ! (M)aster (V)ector (a)ngle (o)ff CFD xy-axis orientation
        END IF
      END IF

! Determine angle of the Master Vector vertical component relative to the
! CFD xy-plane:

      MVvroCFDXY = SQRT( MVComponentsCFD(1) * MVComponentsCFD(1) +
     .                   MVComponentsCFD(2) * MVComponentsCFD(2) )
      IF ((   MVvroCFDXY      .EQ. 0.0) .AND.
     .    (MVComponentsCFD(3) .EQ. 0.0)) THEN
        VldMVOrientXYZ(2) = 0                                            ! No valid vertical orientation of master vector
        MVaoCFDXYZ(2)     = 0.0
      ELSE
        VldMVOrientXYZ(2) = 1                                            ! Valid vertical orientation of master vector can be determined
        AbsMVaoCFDZ = ABS(ATAN2(ABS(MVComponentsCFD(3)),
     .                          ABS(MVvroCFDXY)))
        IF (MVComponentsCFD(3) >= 0.0) THEN
          MVaoCFDXYZ(2) =  AbsMVaoCFDZ * (180.0 / 3.14159)               ! (M)aster (V)ector (a)ngle (o)ff CFD xy-plane
        ELSE
          MVaoCFDXYZ(2) = -AbsMVaoCFDZ * (180.0 / 3.14159)               ! (M)aster (V)ector (a)ngle (o)ff CFD xy-plane
        END IF
      END IF

! Determine angle of the Master Vector relative to the
! Reference Vector's xy-axis:

      IF ((VldRVOrientXYZ(1) .EQ. 1) .AND.
     .    (VldMVOrientXYZ(1) .EQ. 1)) THEN
        MVaoRVXY = MVaoCFDXYZ(1) - RVaoCFDXYZ(1)                         ! (M)aster (V)ector (a)ngle (o)ff (R)eference (V)ector's xy-axis

        IF      (MVaoRVXY < 0.0) THEN
          MVaoRVXY = MVaoRVXY + 360.0
        ELSE IF (MVaoRVXY > 360.0) THEN
          MVaoRVXY = MVaoRVXY - 360.0
        END IF
        IF (MVaoRVXY > 180.0) MVaoRVXY = MVaoRVXY - 360.0                ! Make angles go from 0.0 to 180.0 and 0.0 to -180.0
      ELSE
        MVaoRVXY = 0.0                                                   ! (M)aster (V)ector (a)ngle (o)ff (R)eference (V)ector's xy-axis
      END IF

! Determine velocity of the Master Vector relative to the
! Reference Vector's xyz-axis velocity:

      IF (VectorType .EQ. 0) THEN                                        ! Skip velocity calcs if vector is a 'position vector' (and not a 'velocity vector')
        MVvoRVXYZ(1:3) = 0.0
        GOTO 2332
      END IF
      DO XYZ=1,3
        MVvoRVCFDXYZ(XYZ) = MVComponentsCFD(XYZ) - RVvoCFDXYZ(XYZ)       ! RVvoCFDXYZ(XYZ) units are m/s, so calcs here only valid if MVComponentsCFD(XYZ) units are also m/s
      END DO
      MVvroRVXY = SQRT( MVvoRVCFDXYZ(1) * MVvoRVCFDXYZ(1) +              ! (M)aster (V)ector (v)elocity (r)esultant (o)ff (R)eference (V)ector's xy-axis
     .                  MVvoRVCFDXYZ(2) * MVvoRVCFDXYZ(2) )
      IF      ((MVaoRVXY >=  0.0) .AND.
     .         (MVaoRVXY <= 90.0)) THEN
        MVaoRVXYRad  =          MVaoRVXY  * (3.14159 / 180.0)
        MVvoRVXYZ(1) =  MVvroRVXY * COS(MVaoRVXYRad)                     ! (M)aster (V)ector (v)elocity (o)ff (R)eference (V)ector's xy-axis
        MVvoRVXYZ(2) =  MVvroRVXY * SIN(MVaoRVXYRad)                     ! (M)aster (V)ector (v)elocity (o)ff (R)eference (V)ector's xy-axis
      ELSE IF ((MVaoRVXY >=  90.0) .AND.
     .         (MVaoRVXY <= 180.0)) THEN
        MVaoRVXYRad  = (180.0 - MVaoRVXY) * (3.14159 / 180.0)
        MVvoRVXYZ(1) = -MVvroRVXY * COS(MVaoRVXYRad)
        MVvoRVXYZ(2) =  MVvroRVXY * SIN(MVaoRVXYRad)
      ELSE IF ((MVaoRVXY <=   0.0) .AND.
     .         (MVaoRVXY >= -90.0)) THEN
        MVaoRVXYRad  =         -MVaoRVXY  * (3.14159 / 180.0)
        MVvoRVXYZ(1) =  MVvroRVXY * COS(MVaoRVXYRad)
        MVvoRVXYZ(2) = -MVvroRVXY * SIN(MVaoRVXYRad)
      ELSE IF ((MVaoRVXY <=  -90.0) .AND.
     .         (MVaoRVXY >= -180.0)) THEN
        MVaoRVXYRad  = (180.0 + MVaoRVXY) * (3.14159 / 180.0)
        MVvoRVXYZ(1) = -MVvroRVXY * COS(MVaoRVXYRad)
        MVvoRVXYZ(2) = -MVvroRVXY * SIN(MVaoRVXYRad)
      ELSE
        WRITE(*,*) 'vect.Rel. error: Invalid MVaoRVXY value', MVaoRVXY
        STOP
      END IF
      MVvoRVXYZ(3) = MVvoRVCFDXYZ(3)
 2332 CONTINUE

      RETURN
      END
