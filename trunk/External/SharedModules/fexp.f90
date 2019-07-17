FUNCTION fexp(x)
!-----------------------------------------------------------------------
!
!    File:         fexp
!    Module:       fexp
!    Type:         Function
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/28/02
!
!    Purpose:
!      To exponentiate a given value without over or underflow
!
!----------------------------------------------------------------------c

USE wlKindModule, ONLY: dp

IMPLICIT NONE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

REAL(dp), INTENT(in)  :: x               ! function argument

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(dp)              :: fexp            ! declare function

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

REAL(dp), PARAMETER   :: expmax = 300.d0 ! maximum absolute value of function argument
REAL(dp), PARAMETER   :: expmin = -300.d0 ! maximum absolute value of function argument

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

fexp            = DEXP( DMIN1( expmax, DMAX1( expmin, x ) ) )

END FUNCTION fexp
