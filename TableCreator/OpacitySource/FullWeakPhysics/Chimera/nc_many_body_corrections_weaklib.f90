SUBROUTINE nc_manybody_corrections_weaklib( brydns, tmev, Yp, S_tot )
!-----------------------------------------------------------------------
!
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         4/06/21
!
!    Purpose:
!      To compute the modification of the neutral current scattering
!       rates due to many-body effects.
!
!    Subprograms called:
!       none
!
!    Input arguments:
!  brydns         : baryon density [fm^{-3}]
!  tmev           : temperature [MeV]
!  Yp             : proton fraction
!
!    Output arguments:
!  S_tot          : many-body modification of the axial vector response
!
!-----------------------------------------------------------------------

USE wlKindModule, ONLY: dp

!USE math_functions_module, ONLY: fexp

IMPLICIT none

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

REAL(dp), INTENT(in)        :: brydns ! baryon density [fm^{-3}]
REAL(dp), INTENT(in)        :: tmev   ! temperature [MeV]
REAL(dp), INTENT(in)        :: Yp     ! proton fraction

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(dp), INTENT(out)       :: S_tot  ! many-body modification of the axial vector response

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

REAL(dp), PARAMETER         :: A_0 = 920.d0
REAL(dp), PARAMETER         :: B_0 = 3.05d0
REAL(dp), PARAMETER         :: C_0 = 6140.d0
REAL(dp), PARAMETER         :: D_0 = 1.5d+13
REAL(dp), PARAMETER         :: ga = 1.26d0   ! axial charge of the nucleon
REAL(dp), PARAMETER         :: ga52 = 5.d0 * ga**2
REAL(dp)                    :: A             ! fitting formula parameter
REAL(dp)                    :: B             ! fitting formula parameter
REAL(dp)                    :: C             ! fitting formula parameter
REAL(dp)                    :: S_A           ! fitting formula parameter
REAL(dp)                    :: S_V           ! fitting formula parameter

REAL(dp)                    :: fexp                 ! exponential

EXTERNAL fexp

!-----------------------------------------------------------------------
!
!         \\\\\ GENERAL NEUTRINO-NUCLEON WEAK MAGNETISM /////
!         \\\\\       AND RECOIL CORRECTION TEMRS       /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Dimensionless energies
!-----------------------------------------------------------------------

A       = A_0 * ( brydns * ( 1.d0 - Yp + Yp**2 ) )/( tmev**1.22d0 )
B       = b_0/( tmev**0.75d0 )
C       = C_0 * ( brydns * Yp * ( 1.d0 - Yp ) )/ ( tmev**0.5d0 ) + D_0 * brydns**4/( tmev**6 )
S_A     = 1.d0/( 1.d0 + A * ( 1.d0 + B * fexp(-C) ) )
S_V     = 1.d0

S_tot   = ( ga52 * S_A + ( 1.d0 - Yp ) * S_V )/( ga52 + ( 1.d0 - Yp ) )

RETURN
END SUBROUTINE nc_manybody_corrections_weaklib
