SUBROUTINE abemfrnp( n, e_in, rho, t, xneut, xprot, cmpe, absornp, emitnp )
!-----------------------------------------------------------------------
!
!    File:         abemfrnp
!    Module:       abemfrnp
!    Type:         subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         8/15/93
!
!    Purpose:
!      To compute the inverse mean free paths for the absorption
!      and emission of n-type neutrinos on free neutrons and
!      protons. The inverse mean free paths are given in Bruenn,
!      1985, Ap. J. Suupl., 58, 771. Nucleon final state
!      occupancies are not included.
!
!    Subprograms called:
!       none
!
!    Input arguments:
!  n           : neutrino type (1, e-neutrino; 2, e-antineutrino; 3, t-neutrino)
!  e_in        : neutrino energy [MeV]!       
!  rho         : matter density [g cm^{-3}]!  t           : matter temperature [K]
!  xneut       : free neutron mass fraction
!  xprot       : free proton mass fraction
!  cmpe        : electron chemical potential (including rest mass) [MeV]
!
!    Output arguments:
!  absor       : absorption inverse mean free path (/cm)
!  emitnp      : emission inverse mean free path (/cm)
!
!    Modules used:
!  wlKindModule
!  numerical_module
!  physcnst_module
! 
!  prb_cntl_module
!-----------------------------------------------------------------------

USE wlKindModule, ONLY: dp
USE numerical_module, ONLY: zero, one, pi
USE physcnst_module, ONLY: cvel, hbar, gv, ga, mp, me, Gw, kmev, dmnp, rmu

USE prb_cntl_module, ONLY: iaefnp, rhoaefnp

IMPLICIT none

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(IN)         :: n             ! neutrino flavor index

REAL(dp), INTENT(in)    :: e_in          ! zone centered incoming neutrino energy [MeV]
REAL(dp), INTENT(in)    :: rho           ! density (g/cm^3)
REAL(dp), INTENT(in)    :: t             ! temperature [K]
REAL(dp), INTENT(in)    :: xneut         ! free neutron mass fraction
REAL(dp), INTENT(in)    :: xprot         ! free proton mass fraction
REAL(dp), INTENT(in)    :: cmpe          ! electron chemical porential

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(dp), INTENT(out)   :: absornp       ! inverse mean free path for absorption on free nucleons
REAL(dp), INTENT(out)   :: emitnp        ! inverse mean free path for emission from free nucleons

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

REAL(dp), PARAMETER     :: g2     = ( Gw/mp**2 )**2 * ( hbar * cvel )**2
REAL(dp), PARAMETER     :: coef   = g2 * ( gv**2 + 3.d0 * ga**2 )/pi
REAL(dp), PARAMETER     :: me2    = me * me
REAL(dp)                :: tmev          !  temperature [MeV]
REAL(dp)                :: ron           !  real or effective neutron number
REAL(dp)                :: rop           !  real or effective proton number
REAL(dp)                :: enupdm        !  e_in + dmnp
REAL(dp)                :: enupdm2       !  enupdm^2
REAL(dp)                :: pemitnp       !  ( emission mean free path )/( rop*F_e )
REAL(dp)                :: femitnp       !  ( emission mean free path )/( F_e )
REAL(dp)                :: femitnpl      !  ln(femitnp)
REAL(dp)                :: F_em          !  1 - F_e, 1/( 1 + exp( ( enu + dmnp - cmpe )/t ) )
REAL(dp)                :: F_el          !  ln(F_e)
REAL(dp)                :: eta           !  ( e_in + dmnp - cmpe )/tmev
REAL(dp)                :: ex            !  exp(eta)
REAL(dp)                :: cmppos        !  positron chemical potential
REAL(dp)                :: F_p           !  1/( 1 + exp( ( e_in - dmnp - cmppos )/t ) )
REAL(dp)                :: enumdm        !  e_in - dmnp
REAL(dp)                :: enumdm2       !  enumdm^2

REAL(dp)                :: emitnpl       !  ln(emitnp)

REAL(dp), PARAMETER     :: expmax = 3.d02

REAL(dp), EXTERNAL      :: fexp

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  .Return if iaefnp = 0 or rho  >  rhoaefnp or n > 2
!-----------------------------------------------------------------------

IF ( iaefnp == 0  .or.  rho  >  rhoaefnp  .or.  n > 2 ) THEN
  emitnp             = zero
  absornp            = zero
  RETURN
END IF

!-----------------------------------------------------------------------
!  Set emitnp and absornp to zero and return if both xneut and
!   xprot are zero.
!-----------------------------------------------------------------------

IF ( xneut == zero  .and.  xprot == zero ) THEN
  emitnp             = zero
  absornp            = zero
  RETURN
END IF

!-----------------------------------------------------------------------
!  ron : number of free neutrons per unit volume (/cm**3).
!  rop : number of free protons per unit volume (/cm**3).
!-----------------------------------------------------------------------

tmev                 = kmev * t
ron                  = xneut * rho/rmu
rop                  = xprot * rho/rmu

!-----------------------------------------------------------------------
!
!    \\\\\ ABSORPTION AND EMISSION INV MFP'S FOR E-NEUTRINOS /////
!
!-----------------------------------------------------------------------

ex                   = zero

!-----------------------------------------------------------------------
!  pemitnp: ( emission mean free path )/( rop*F_e ).
!  femitnp: ( emission mean free path )/( F_e ).
!-----------------------------------------------------------------------

IF ( n == 1 ) THEN

  enupdm             = e_in + dmnp
  enupdm2            = enupdm * enupdm
  pemitnp            = coef * enupdm2 * dsqrt( one - me2/enupdm2 )
  femitnp            = rop * pemitnp

!-----------------------------------------------------------------------
!  Case that the free proton fraction is zero
!-----------------------------------------------------------------------

  IF ( femitnp == zero ) THEN
    emitnp           = zero
    eta              = ( e_in + dmnp - cmpe )/tmev
    ex               = fexp(eta)
    F_em             = ex/( one + ex )
    absornp          = ron * pemitnp * F_em
    RETURN
  END IF

!-----------------------------------------------------------------------
!  General case compute logs of inv mfp's
!-----------------------------------------------------------------------

  femitnpl           = dlog(femitnp)
  eta                = ( e_in + dmnp - cmpe )/tmev
  IF ( eta > expmax ) THEN
    F_el             = -eta
    F_em             = one
  ELSE
    ex               = fexp(eta) 
    F_el             = -DLOG( one + ex )
    F_em             = ex/( one + ex )
  END IF

!-----------------------------------------------------------------------
!  Log of emission inverse mean free path; emission inverse mean free
!   path.
!-----------------------------------------------------------------------

  emitnpl            = femitnpl + F_el
  emitnp             = fexp(emitnpl)

!-----------------------------------------------------------------------
!  Compute the asorption inverse mean free path
!-----------------------------------------------------------------------

  absornp            = ron * pemitnp * F_em

  RETURN
  
END IF

!-----------------------------------------------------------------------
!
!  \\\\\ ABSORPTION AND EMISSION INV MFP'S FOR E-ANTINEUTRINOS /////
!
!-----------------------------------------------------------------------

IF ( n == 2 ) THEN

  IF ( e_in - me <= DABS(dmnp) ) THEN
    emitnp           = zero
    absornp          = zero
    RETURN
  END IF

!-----------------------------------------------------------------------
!  Absorption and emission inv mfp's for e-antineutrinos
!
!  pemitnp: ( emission mean free path )/( ron*F_p )
!-----------------------------------------------------------------------

  enumdm             = e_in - dmnp
  enumdm2            = enumdm * enumdm
  pemitnp            = coef * enumdm2 * dsqrt( one - me2/enumdm2 )

!-----------------------------------------------------------------------
!  Positron occupation number for positrons of energy e_in - dmnp
!-----------------------------------------------------------------------

  cmppos             = -cmpe
  eta                = ( e_in - dmnp - cmppos )/tmev
  ex                 = fexp(eta)
  F_p                = one/( one + ex )

!-----------------------------------------------------------------------
!  inv mfp's
!-----------------------------------------------------------------------

  emitnp             = ron * pemitnp * F_p
  absornp            = rop * pemitnp * ( one - F_p )

  RETURN
  
END IF

END SUBROUTINE abemfrnp
