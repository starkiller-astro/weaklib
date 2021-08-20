FUNCTION ff( x ,y )
!-----------------------------------------------------------------------
!
!    File:         ff
!    Module:       ff
!    Type:         Function
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/28/02
!
!    Purpose:
!      To multiply two fermi functions
!
!----------------------------------------------------------------------c

USE wlKindModule, ONLY: dp
USE numerical_module, ONLY: one

IMPLICIT NONE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

REAL(dp), INTENT(IN)  :: x               ! function argument
REAL(dp), INTENT(IN)  :: y               ! function argument

!-----------------------------------------------------------------------
!        Output variables
!-----------------------------------------------------------------------

REAL(dp)              :: ff              ! declare function

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

REAL(dp)              :: fexp            ! exponential

EXTERNAL fexp

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

ff              = one/( ( one + fexp(x) ) * ( one + fexp(y) ) )

END FUNCTION ff
