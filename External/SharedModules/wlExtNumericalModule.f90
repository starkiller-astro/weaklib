module wlExtNumericalModule

USE wlKindModule, ONLY: dp

REAL(dp), PARAMETER  :: zero      = 0.0d0
REAL(dp), PARAMETER  :: half      = 0.5d0
REAL(dp), PARAMETER  :: one       = 1.0d0
REAL(dp), PARAMETER  :: epsilon   = 1.0d-100

REAL(dp), PARAMETER  :: pi        = 3.1415926535897932385d0 ! pi

!-----------------------------------------------------------------------
!  Derived constants
!-----------------------------------------------------------------------

REAL(dp), PARAMETER  :: pi2       = pi * pi
REAL(dp), PARAMETER  :: twpi      = 2.d0 * pi
REAL(dp), PARAMETER  :: frpi      = 4.d0 * pi
REAL(dp), PARAMETER  :: third     = 1.d0/3.d0
REAL(dp), PARAMETER  :: twothd    = 2.d0 * third
REAL(dp), PARAMETER  :: fourthd   = 4.d0 * third
REAL(dp), PARAMETER  :: frpith    = third * frpi
REAL(dp), PARAMETER  :: sxtnpi    = 16.d0 * pi

!-----------------------------------------------------------------------
!  Equation of state coefficients
!
!   kfm : [ # nucleons g^{-1} ][ cm^{3}fm^{-3} ]
!-----------------------------------------------------------------------

REAL(dp)             ::  log_e
REAL(dp)             ::  ln_10

END module wlExtNumericalModule
