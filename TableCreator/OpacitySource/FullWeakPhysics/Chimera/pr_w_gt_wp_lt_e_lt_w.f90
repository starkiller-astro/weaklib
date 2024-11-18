SUBROUTINE pr_w_gt_wp_lt_e_lt_w( xl, xu, eta, beta, j0i, j0ii, j1i, j1ii )
!-----------------------------------------------------------------------
!
!    File:         pr_w_gt_wp_lt_e_lt_w
!    Module:       pr_w_gt_wp_lt_e_lt_w
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/05/03
!
!    Purpose:
!      To integrate for the case w > w', w' < e < w, the quantities
!
!          Fe(Ee)*Febar(Enu + Enubar - Ee)*Phi(Enu,Enubar,Ee)
!
!      and
!
!          [1-Fe(Ee)]*[1-Febar(Enu + Enubar - Ee)]*Phi(Enu,Enubar,Ee)
!
!  e       : (electron energy)/kt    (integration variable)
!  w       : (in beam neutrino energy)/kt
!  wp      : (out beam neutrino energy)/kt
!  eta     : (electron chemical potential - mc2)/kt
!
!    Subprograms called:
!      none
!
!    Input arguments:
!  xl      : lower limit of integration
!  eta     : (electron chemical potential - mc2)/kt
!
!    Output arguments:
!  j0i     : zero moment of the "i" in neutrino scattering function
!  j0ii    : zero moment of the "ii" in neutrino scattering function
!  j1i     : first moment of the "i" in neutrino scattering function
!  j1ii    : first moment of the "ii" in neutrino scattering function
!
!    Variables that must be passed through common:
!      none
!
!    Modules:
!  wlKindModule, numerical_module
!  nes_module
!
!-----------------------------------------------------------------------

USE wlKindModule, ONLY: dp
USE numerical_module, ONLY: zero, half, one

USE nes_module
USE pair_module, ONLY: nleg, xe, wte

IMPLICIT none

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

REAL(dp), INTENT(in)      :: xl            ! lower limit of integration
REAL(dp), INTENT(in)      :: xu            ! upper limit of integration
REAL(dp), INTENT(in)      :: eta           ! (electron chemical potential - mc2)/kt
REAL(dp), INTENT(IN)      :: beta          ! 1/kt

!-----------------------------------------------------------------------
!        Input-Output variables.
!-----------------------------------------------------------------------

REAL(dp), INTENT(inout)   :: j0i           ! partial integral contributed by this subroutine
REAL(dp), INTENT(inout)   :: j0ii          ! partial integral contributed by this subroutine
REAL(dp), INTENT(inout)   :: j1i           ! partial integral contributed by this subroutine
REAL(dp), INTENT(inout)   :: j1ii          ! partial integral contributed by this subroutine

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                            :: i             ! summation index

REAL(dp)                  :: h0i           ! integrand of zero moment of pair function i
REAL(dp)                  :: h0ii          ! integrand of zero moment of pair function ii

REAL(dp)                  :: su0i          ! partial integral
REAL(dp)                  :: su0ii         ! partial integral
REAL(dp)                  :: su1i          ! partial integral
REAL(dp)                  :: su1ii         ! partial integral

REAL(dp)                  :: e_mid         ! midpoint energy
REAL(dp)                  :: e_del         ! half the energy width
REAL(dp)                  :: e_var         ! scaled integration point

REAL(dp)                  :: arg1          ! argument of an exponential
REAL(dp)                  :: arg2          ! argument of an exponential
REAL(dp)                  :: exp1          ! value of an exponential
REAL(dp)                  :: exp2          ! value of an exponential
REAL(dp)                  :: eblock        ! blocking factor for neutrino-antineutrino annihilation

REAL(dp)                  :: fexp          ! exponential

EXTERNAL fexp

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Return if xl = xu
!-----------------------------------------------------------------------

IF ( xl == xu ) RETURN

!-----------------------------------------------------------------------
!  Initialize for integration
!-----------------------------------------------------------------------

su0i               = zero
su0ii              = zero
su1i               = zero
su1ii              = zero

e_mid              = half * ( xu + xl )
e_del              = half * ( xu - xl )

CALL gquad( nleg, xe, wte, nleg )

!-----------------------------------------------------------------------
!  Integrate
!-----------------------------------------------------------------------

DO i = 1,nleg

  e_var            = xe(i) * e_del
  e                = e_mid + e_var

  arg1             = beta * e - eta
  exp1             = fexp(arg1)
  arg2             = beta * ( w_p_wp - e ) + eta
  exp2             = fexp(arg2)
  eblock           = ( exp1/( exp1 + one ) ) * ( exp2/( exp2 + one ) )

  h0i              = r4_15 * wp5 - r4_3 * wp4 * e + r8_3 * wp3 * e**2
  h0ii             =   wp3 * ( r8_3 * w2 + 4.d0 * w * wp + 1.6d0 * wp2 ) - e * wp * ( r16_3 * wp2 * w + 4.d0 * wp3 ) &
&                  + r8_3 * e**2 * wp3

  su0i             = su0i  + eblock * h0i  * wte(i)
  su0ii            = su0ii + eblock * h0ii * wte(i)

END DO

su0i               = e_del * su0i
su0ii              = e_del * su0ii
su1i               = e_del * su1i
su1ii              = e_del * su1ii


j0i                = j0i  + su0i
j0ii               = j0ii + su0ii
j1i                = j1i  + su1i
j1ii               = j1ii + su1ii

RETURN
END SUBROUTINE pr_w_gt_wp_lt_e_lt_w
