SUBROUTINE abemfrnpetanp( n, e_in, rho, t, xneut, xprot, cmpn, cmpp, cmpe, &
& absornp, emitnp)
!-----------------------------------------------------------------------
!
!    File:         abemfrnpetanp
!    Module:       abemfrnpetanp
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
!      occupancies are included.
!
!    Subprograms called:
!  terminate   : terminates the run in the event of an error
!
!    Input arguments:
!  n           : neutrino type (1, e-neutrino; 2, e-antineutrino; 3, t-neutrino)
!  e_in        : neutrino energy [MeV]!       
!  rho         : matter density [g cm^{-3}]!  t           : matter temperature [K]
!  xneut       : free neutron mass fraction
!  xprot       : free proton mass fraction
!  cmpn        : free neutron chemical potential (excluding rest mass) [MeV]
!  cmpp        : free proton chemical potential (excluding rest mass) [MeV]
!  cmpe        : electron chemical potential (including rest mass) [MeV]
!
!    Output arguments:
!  absor       : absorption inverse mean free path (/cm)
!  emitnp      : emission inverse mean free path (/cm)!  
!
!    Modules used:
!  wlKindModule
!  numerical_module
!  physcnst_module
! 
!  prb_cntl_module
!
!-----------------------------------------------------------------------

USE wlKindModule, ONLY: dp
USE numerical_module, ONLY: zero, one, pi
USE physcnst_module, ONLY: cvel, hbar, gv, ga, mp, me, Gw, kmev, dmnp, rmu

USE prb_cntl_module, ONLY: iaefnp, rhoaefnp

IMPLICIT none

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)         :: n             ! neutrino flavor index

REAL(dp), INTENT(in)    :: e_in          ! zone centered incoming neutrino energy [MeV]
REAL(dp), INTENT(in)    :: rho           ! density (g/cm^3)
REAL(dp), INTENT(in)    :: t             ! temperature [K]
REAL(dp), INTENT(in)    :: xneut         ! free neutron mass fraction
REAL(dp), INTENT(in)    :: xprot         ! free proton mass fraction
REAL(dp), INTENT(in)    :: cmpn          ! neutron chemical porential
REAL(dp), INTENT(in)    :: cmpp          ! proton chemical porential
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
REAL(dp)                :: etam          !  ( e_in + dmnp + cmpn - cmpp - cmpe )/tmev 
REAL(dp)                :: cmppos        !  positron chemical potential
REAL(dp)                :: F_p           !  1/( 1 + exp( ( e_in - dmnp - cmppos )/t ) )
REAL(dp)                :: enumdm        !  e_in - dmnp
REAL(dp)                :: enumdm2       !  enumdm^2

REAL(dp)                :: absornpl      !  ln(absornp)
REAL(dp)                :: emitnpl       !  ln(emitnp)

REAL(dp), PARAMETER     :: expmax = 3.d02
REAL(dp), PARAMETER     :: e_min  = 1.d-100

REAL(dp), EXTERNAL      :: fexp

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Return if iaefnp = 0 or rho  >  rhoaefnp or n > 2
!-----------------------------------------------------------------------

IF ( iaefnp == 0  .or.  rho  >  rhoaefnp  .or.  n > 2 ) THEN
  emitnp             = zero
  absornp            = zero
  RETURN
END IF

!-----------------------------------------------------------------------
!  Set emitnp and absornp to zero and return if both xneut and xprot are
!   zero.
!-----------------------------------------------------------------------

IF ( xneut == zero  .and.  xprot == zero ) THEN
  emitnp             = zero
  absornp            = zero
  RETURN
END IF

!-----------------------------------------------------------------------
!  ron : number of free neutrons per unit volume (/cm**3) corrected for
!         final state free proton occupancy.
!  rop : number of free protons per unit volume (/cm**3) corrected for
!         final state free neutron occupancy.
!-----------------------------------------------------------------------

tmev                 = kmev * t

IF ( DABS( cmpn - cmpp ) < dmnp ) then

!-----------------------------------------------------------------------
!  Compute neutron and proton number
!-----------------------------------------------------------------------

  ron                = ( xneut/rmu ) * rho
  rop                = ( xprot/rmu ) * rho

ELSE

!-----------------------------------------------------------------------
!  Compute corrected neutron and proton number
!-----------------------------------------------------------------------

  ron                = DABS( ( rho/rmu ) * ( xprot - xneut )/( fexp( ( cmpp - cmpn )/tmev ) - one ) )
  rop                = DABS( ( rho/rmu ) * ( xneut - xprot )/( fexp( ( cmpn - cmpp )/tmev ) - one ) )

END IF

!-----------------------------------------------------------------------
!
!    \\\\\ ABSORPTION AND EMISSION INV MFP'S FOR E-NEUTRINOS /////
!
!-----------------------------------------------------------------------

IF ( n == 1 ) THEN

  enupdm             = e_in + dmnp
  enupdm2            = enupdm * enupdm
  pemitnp            = coef * enupdm2 * dsqrt( one - me2/enupdm2 )
  femitnp            = rop * pemitnp

!-----------------------------------------------------------------------
!  Case that the free proton fraction is zero
!-----------------------------------------------------------------------

  IF ( femitnp < e_min ) THEN
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
  ELSE
    ex               = fexp(eta)
    F_el             = -DLOG( one + ex )
  END IF

!-----------------------------------------------------------------------
!  Log of emission inverse mean free path; emission inverse mean free
!   path.
!-----------------------------------------------------------------------

  emitnpl            = femitnpl + F_el
  emitnp             = fexp(emitnpl)

!-----------------------------------------------------------------------
!  Log of asorption inverse mean free path (detailed balance);
!   asorption inverse mean free path.
!
!  absornp = dexp( etam )*emitnp
!  etam  = ( e_in + dmnp + cmpn - cmpp - cmpe )/tmev
!-----------------------------------------------------------------------

  etam               = ( e_in + dmnp + cmpn - cmpp - cmpe )/tmev
  absornpl           = emitnpl + etam
  absornp            = fexp(absornpl)
  RETURN
  
END IF

!-----------------------------------------------------------------------
!
!  \\\\\ ABSORPTION AND EMISSION INV MFP'S FOR E-ANTINEUTRINOS /////
!
!-----------------------------------------------------------------------

IF ( n == 2 ) THEN

!-----------------------------------------------------------------------
!  Return if the neutrino energy is below the threshold for absorption
!   on fre protons.
!-----------------------------------------------------------------------

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

STOP
END SUBROUTINE abemfrnpetanp
