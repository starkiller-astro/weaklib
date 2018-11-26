MODULE vector_functions_module

PUBLIC f10_vec, abflog_vec, fexp_vec

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

CONTAINS

FUNCTION abflog_vec(x,d1,d2)

!----------------------------------------------------------
!
!    File:         abflog_vec
!    Module:       vector_fucntions_module
!    Type:         Function
!    Author:       E. J. Lentz, U. Tennessee, Knoxville
!
!    Date:         5/20/10
!
!    Purpose:
!      To take the base 10 log of the abs of a given array without underflow
!
!----------------------------------------------------------------------

USE kind_module
USE numerical_module, ONLY: epsilon

IMPLICIT NONE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in) :: d1, d2  ! dimensions
REAL(KIND=double), DIMENSION(d1,d2), INTENT(in)  :: x               ! function argument

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), DIMENSION(d1,d2)              :: abflog_vec      ! declare function

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER :: n1, n2  ! loop counters
REAL(KIND=double), DIMENSION(d1,d2) :: fabs

!-----------------------------------------------------------------------

DO n2 = 1,d2
  DO n1 = 1,d1
    fabs(n1,n2)          = DABS( x(n1,n2) )
  END DO !n1 loop
END DO !n2 loop

DO n2 = 1,d2
  DO n1 = 1,d1
    abflog_vec(n1,n2)    = DLOG10( fabs(n1,n2) + epsilon )
  END DO !n1 loop
END DO !n2 loop

END FUNCTION abflog_vec

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

FUNCTION f10_vec(x,d1,d2)
!-----------------------------------------------------------------------
!
!    File:         f10_vec
!    Module:       vector_fucntions_module
!    Type:         Function
!    Author:       E. J. Lentz, U. Tennessee, Knoxville
!
!    Date:         3/29/10
!
!    Purpose:
!      To raise a given array to the power of 10 without over or underflow
!
!----------------------------------------------------------------------

USE kind_module

IMPLICIT NONE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in) :: d1, d2  ! dimensions
REAL(KIND=double), DIMENSION(d1,d2), INTENT(in)  :: x               ! function argument

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), DIMENSION(d1,d2)              :: f10_vec         ! declare function

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

REAL(KIND=double), PARAMETER   :: expmax  =  300.d0 ! maximum absolute value of function argument
REAL(KIND=double), PARAMETER   :: nexpmax = -300.d0 ! (-) maximum absolute value of function argument

INTEGER :: n1, n2  ! loop counters
REAL(KIND=double), DIMENSION(d1,d2) :: fmax, fmin

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

DO n2 = 1,d2
  DO n1 = 1,d1
    fmax(n1,n2)     = DMAX1( nexpmax, x(n1,n2) )
  END DO !n1 loop
END DO !n2 loop

DO n2 = 1,d2
  DO n1 = 1,d1
    fmin(n1,n2)     = DMIN1( expmax, fmax(n1,n2) )
  END DO !n1 loop
END DO !n2 loop

DO n2 = 1,d2
  DO n1 = 1,d1
    f10_vec(n1,n2)  = 10.d0 ** fmin(n1,n2)
  END DO !n1 loop
END DO !n2 loop


END FUNCTION f10_vec

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

FUNCTION fexp_vec(x,d1,d2)
!----------------------------------------------------------
!
!    Function:     fexp_vec
!    Module:       vector_fucntions_module
!    Type:         Function
!    Author:       E. J. Lentz, U. Tennessee, Knoxville
!
!    Date:         5/20/10
!
!    Purpose:
!      To exponentiate a given array without over or underflow
!
!----------------------------------------------------------------------

USE kind_module

IMPLICIT NONE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in) :: d1, d2  ! dimensions
REAL(KIND=double), DIMENSION(d1,d2), INTENT(in)  :: x               ! function argument

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), DIMENSION(d1,d2)              :: fexp_vec        ! declare function

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

REAL(KIND=double), PARAMETER   :: expmax  =  300.d0 ! maximum absolute value of function argument
REAL(KIND=double), PARAMETER   :: nexpmax = -300.d0 ! (-) maximum absolute value of function argument

INTEGER :: n1, n2  ! loop counters
REAL(KIND=double), DIMENSION(d1,d2) :: fmax, fmin

!-----------------------------------------------------------------------

DO n2 = 1,d2
  DO n1 = 1,d1
    fmax(n1,n2)     = DMAX1( nexpmax, x(n1,n2) )
  END DO !n1 loop
END DO !n2 loop

DO n2 = 1,d2
  DO n1 = 1,d1
    fmin(n1,n2)     = DMIN1( expmax, fmax(n1,n2) )
  END DO !n1 loop
END DO !n2 loop

DO n2 = 1,d2
  DO n1 = 1,d1
    fexp_vec(n1,n2) = DEXP( fmin(n1,n2) )
  END DO !n1 loop
END DO !n2 loop


END FUNCTION fexp_vec


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

END MODULE vector_functions_module
