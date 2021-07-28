!-----------------------------------------------------------------------
!    Module:       numerical_module
!    Author:       S. W. Bruenn
!    Date:         8/16/02
!
!    Common numerical constants
!-----------------------------------------------------------------------

module numerical_module

USE kind_module

REAL(double), PARAMETER  :: zero      = 0.0d0
REAL(double), PARAMETER  :: half      = 0.5d0
REAL(double), PARAMETER  :: one       = 1.0d0
REAL(double), PARAMETER  :: epsilon   = 1.0d-100

REAL(double), PARAMETER  :: pi        = 3.1415926535897932385d0 ! pi

!-----------------------------------------------------------------------
!  Derived constants
!-----------------------------------------------------------------------

REAL(double), PARAMETER  :: pi2       = pi * pi
REAL(double), PARAMETER  :: twpi      = 2.d0 * pi
REAL(double), PARAMETER  :: frpi      = 4.d0 * pi
REAL(double), PARAMETER  :: third     = 1.d0/3.d0
REAL(double), PARAMETER  :: twothd    = 2.d0 * third
REAL(double), PARAMETER  :: fourthd   = 4.d0 * third
REAL(double), PARAMETER  :: frpith    = third * frpi
REAL(double), PARAMETER  :: sxtnpi    = 16.d0 * pi

!-----------------------------------------------------------------------
!  Equation of state coefficients
!
!   kfm : [ # nucleons g^{-1} ][ cm^{3}fm^{-3} ]
!-----------------------------------------------------------------------

REAL(double)             ::  log_e
REAL(double)             ::  ln_10

END module numerical_module
