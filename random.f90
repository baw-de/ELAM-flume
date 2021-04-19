!/*--------------------------------------------------------------------*\
! Application
!   random
!
! Purpose
!   Return pseudo-random numbers
!
! This file is provided courtesy of the United States Army Corps of
!   Engineers, Engineer Research and Development Center, Environmental
!   Laboratory. See license file for information on usage.
!
!\*--------------------------------------------------------------------*/

      Subroutine RandomFromSeed(RNDX,SEED)

!Algorithm from:
! "Structured Fortran 77 for Engineers and Scientists, 5th Edition", Author: Delores M. Etter
! Random # is between 0.0 and 1.0

      IMPLICIT NONE

      INTEGER :: SEED

      real*8 :: RNDX


      SEED = 2045*SEED + 1
      SEED = SEED - (SEED/1048576)*1048576
      RNDX = REAL(SEED+1)/1048577.0


!    ! switch between 0 and 1 for tests
!      if (SEED == 1) then
!        SEED = 0
!        RNDX = 0
!      else ! SEED == 0
!        SEED = 1
!        RNDX = 1
!      end if


      RETURN
      END
