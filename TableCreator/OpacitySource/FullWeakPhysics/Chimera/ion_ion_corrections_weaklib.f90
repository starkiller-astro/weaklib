SUBROUTINE ion_ion_corrections_weaklib( rho, t, e_nu, xion, aion, zion, ciicr )
!-----------------------------------------------------------------------
!
!    File:         scatiicr_H
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         6/22/00
!
!    Purpose:
!      To compute correction factors for the isoenergetic
!         scattering rates for ion-ion correlation effects
!         using a prescription given by C. Horowitz.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!
!   rho       : matter density [g cm^{-3}]
!   t         : matter temperature [K]
!   e_nu      : neutrino energy
!   xion      : mass fraction of heavy ions
!   aion      : mass number of heavy ions
!   zion      : charge number of heavy ions
!
!    Output arguments:
!
!   ccicr     : correction factor for ion-ion correlation effects
!
!-----------------------------------------------------------------------

USE wlKindModule, ONLY: dp
USE numerical_module, ONLY: zero, one, third, pi
USE physcnst_module, ONLY: kmev, hbarc, rmu

USE prb_cntl_module

IMPLICIT none

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

REAL(dp), INTENT(in)     :: rho           ! density [g cm^{-3}]
REAL(dp), INTENT(in)     :: t             ! temperature [K]
REAL(dp), INTENT(in)     :: e_nu          ! neutrino energy
REAL(dp), INTENT(in)     :: xion          ! heavy nucleus mass fraction
REAL(dp), INTENT(in)     :: aion          ! heavy nucleus mass number
REAL(dp), INTENT(in)     :: zion          ! heavy nucleus charge number

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(dp), INTENT(OUT)    :: ciicr         ! correlation correction

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                  :: i             ! do index

REAL(dp), PARAMETER      :: cgamma = 1.96d-5 ! e**2/4*pi*rmu (e**2 = 4*pi*alpha (natural units) 
                                                ! = 4*pi*alpha*hbar*c (cgs units) = 4*pi*(4.8e-10)**2
REAL(dp), PARAMETER      :: cx = 9.44d+9
REAL(dp), PARAMETER      :: amin = 1.d-10

REAL(dp)                 :: tmev          ! temperature [MeV]
REAL(dp)                 :: gamma         ! zion**2*e**2/4*pi*a*kt = cgamma*zion**2/tmev*(rho*xion/aion)**1/3
REAL(dp)                 :: gamma12       ! gamma^{1/2}
REAL(dp)                 :: nion          ! number density of ions
REAL(dp)                 :: a             ! (3/4*pi*n)**1/3 = (3*rmu*aion/4*pi*rho*xion)**1/3
REAL(dp)                 :: hbarc_cm      ! hbar (Mev cm)
REAL(dp)                 :: aE            ! a/hbarc_cm
REAL(dp)                 :: ebar          ! e_nu * aE - Wigner cell size/neutrino wavelength
REAL(dp)                 :: estar         ! 3.d0 + 4.d0/gamma12
REAL(dp)                 :: arg           ! argument of exponential

REAL(dp)                 :: beta0
REAL(dp), DIMENSION(6)   :: beta
REAL(dp), DIMENSION(4,4) :: beta3

REAL(dp)                 :: fexp                 ! exponential

EXTERNAL fexp
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Return if xion < amin  or  aion < amin  or  zion < amin
!-----------------------------------------------------------------------

IF ( DMIN1( xion, aion, zion ) < amin ) THEN
  ciicr = one
  RETURN
END IF

!-----------------------------------------------------------------------
!  Fitting coeficients from Horowitz 2002
!-----------------------------------------------------------------------

beta3(1,1) = -7.362056d+0
beta3(1,2) =  0.5371365d0
beta3(1,3) = -0.1078845d0
beta3(1,4) =  4.189612d-3
               
beta3(2,1) =  3.4489581d0
beta3(2,2) = -0.40251656d0
beta3(2,3) =  9.0877878d-2
beta3(2,4) = -3.4353581d-3
               
beta3(3,1) = -0.74128645d0
beta3(3,2) =  0.11019855d0
beta3(3,3) = -2.5359361d-2
beta3(3,4) =  9.0487744d-4
               
beta3(4,1) =  5.9573285d-2
beta3(4,2) = -1.0186552d-2
beta3(4,3) =  2.2791369d-3
beta3(4,4) = -7.4614597d-5

!-----------------------------------------------------------------------
!  Compute gamma, a, and ebar
!-----------------------------------------------------------------------

tmev               = kmev * t
gamma              = cgamma * ( zion**2/tmev ) * ( rho * xion/aion )**third
gamma              = DMAX1( gamma, one )
gamma              = DMIN1( gamma, 150.d0 )
gamma12            = DSQRT(gamma)
nion               = rho * xion/( rmu * aion )
a                  = ( 3.d0/( 4.d0 * pi * nion ) )**third
hbarc_cm           = 1.d-13 * hbarc
aE                 = a/hbarc_cm
ebar               = e_nu * aE

!-----------------------------------------------------------------------
!  ciicr = 1 if ebar > estar
!-----------------------------------------------------------------------

estar              = 3.d0 + 4.d0/gamma12
IF ( ebar >= estar ) THEN
  ciicr = 1.0d0
  RETURN
END IF

!-----------------------------------------------------------------------
!  Compute ciicr if ebar > estar
!-----------------------------------------------------------------------

beta0              = DLOG(0.3d0/( 0.3d0 + 3.d0 * gamma ) )
beta(1)            = zero
beta(2)            = 6.667d0
DO i = 3,6
  beta(i)          = beta3(i-2,1) + beta3(i-2,2) * gamma12 + beta3(i-2,3) * gamma + beta3(i-2,4) &
&                  * gamma * gamma12
END DO

arg                = beta0
DO i = 1,6
  arg              = arg + beta(i) * ebar**(DBLE(i))
END DO

ciicr              = one/( one + fexp( -arg ) )

RETURN
END SUBROUTINE ion_ion_corrections_weaklib
