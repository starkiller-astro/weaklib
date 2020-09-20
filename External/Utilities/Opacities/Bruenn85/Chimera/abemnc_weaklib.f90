SUBROUTINE abemnc_weaklib( n, nez, e_in, rho, t, xh, ah, zh, cmpn, cmpp, cmpe,  &
& absrnc, emitnc )
!-----------------------------------------------------------------------
!
!    Purpose:
!      To computes the inverse mean free path for the absorption of e-type neutrinos
!       on nuclei and the inverse process (emission of e-type neutrinos
!       by electron capture on nuclei). The Fuller, Fowler, Neuman 1982, Ap. J. 252, 715
!       approximation for the heavy nucleus matrix element as given in Bruenn 1985,
!       Ap. J. Suupl., 58, 771 is used.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!
!  n           : neutrino type (1, e-neutrino; 2, e-antineutrino; 3, t-neutrino)
!  nez         : number of energy groups
!  e_in        : neutrino energy [MeV]
!  rho         : matter density [g cm^{-3}]!  t           : matter temperature [K]
!  xh          : heavy nuclei mass fraction
!  ah          : heavy nuclei mass number
!  zh          : heavy nuclei charge number
!  cmpn        : free neutron chemical potential (excluding rest mass) [MeV]
!  cmpp        : free proton chemical potential (excluding rest mass) [MeV]
!  cmpe        : electron chemical potential (including rest mass) [MeV]
!                 between the exitation energy of daughter and parent nucleus [MeV]
!
!    Output arguments:
!  absrnc      : absorption inverse mean free path on nuclei (/cm)
!  emitnc      : emission inverse mean free path on nuclei (/cm)
!
!    Modules used:
!  kind_module
!  numerical_module
!  physcnst_module
!
!  prb_cntl_module
!
!-----------------------------------------------------------------------

USE kind_module
USE numerical_module, ONLY: zero, one, epsilon, pi
USE physcnst_module, ONLY: Gw, mp, hbar, cvel, ga, kmev, me, dmnp, rmu

USE prb_cntl_module, ONLY: iaence, edmpe

IMPLICIT none

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(IN)         :: n             ! neutrino flavor index
INTEGER, INTENT(IN)         :: nez           ! number of energy groups

REAL(double), INTENT(in), DIMENSION(nez) :: e_in ! zone centered incoming neutrino energy [MeV]
REAL(double), INTENT(in)    :: rho           ! density (g/cm^3)
REAL(double), INTENT(in)    :: t             ! temperature [K]
REAL(double), INTENT(in)    :: xh            ! heavy nuclei mass fraction
REAL(double), INTENT(in)    :: ah            ! heavy nuclei mass number
REAL(double), INTENT(in)    :: zh            ! heavy nuclei charge number
REAL(double), INTENT(in)    :: cmpn          ! neutron chemical porential
REAL(double), INTENT(in)    :: cmpp          ! proton chemical porential
REAL(double), INTENT(in)    :: cmpe          ! electron chemical porential

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(double), INTENT(out), DIMENSION(nez) :: absrnc ! inverse mean free path for absorption on free nucleons
REAL(double), INTENT(out), DIMENSION(nez) :: emitnc ! inverse mean free path for emission from free nucleons

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                     :: k             ! energy index counter
REAL(double), PARAMETER     :: g2    = ( Gw/mp**2 )**2 * ( hbar * cvel )**2   ! square of the Fermi constant in units of cm^{2} MeV^{-2}
REAL(double), PARAMETER     :: coef  = ( 2.d0/7.d0 ) * ( g2/pi ) * ( ga**2 )  ! coefficient of the emis and abs inverse mean free paths (cm^{-1})
REAL(double), PARAMETER     :: me2   = me * me    ! square of the electron rest mass
REAL(double)                :: nh            ! number of neutrons
REAL(double)                :: coefn         ! number of neutron holes
REAL(double)                :: zm20          ! number of available protons
REAL(double)                :: eelec         ! electron energy [MeV]
REAL(double)                :: eelec2        ! electron energy squared [MeV]
REAL(double)                :: ronc          ! number density of nuclei (cm^{-3})
REAL(double)                :: ceelec        ! coefficient for computing the emission inverse mean free path
REAL(double)                :: pemit         ! emission inverse mean free path per electron
REAL(double)                :: eta           ! (eelec - cmpe)/tmev
REAL(double)                :: F_e           ! electron occupation number
REAL(double)                :: tmev          ! temperature [MeV]
REAL(double)                :: etam          ! ( e_in + dmnp + cmpn - cmpp - cmpe )/tmev
REAL(double), EXTERNAL      :: fexp          ! exponential function
REAL(double), PARAMETER     :: maxeta = 4.d1 ! limit where 1/(1+e(x)) = e^-x (ln(1e16))

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Set rates to zero and return if iaence = 0 or if n ne 1
!-----------------------------------------------------------------------

IF ( iaence == 0  .or.  n /= 1 ) THEN
  emitnc           = zero
  absrnc           = zero
  RETURN
END IF ! iaence == 0  .or.  n /= 1 .or. nse == 0

!-----------------------------------------------------------------------
!  Compute coefn, the number of neutron holes
!-----------------------------------------------------------------------

nh                 = ah - zh
coefn              = DMIN1( 40.d+00 - nh, 6.d+00 )
IF ( coefn < zero ) THEN
  emitnc           = zero
  absrnc           = zero
  RETURN
END IF ! oefn < zero

!-----------------------------------------------------------------------
!  The Bowers and Wilson prescription is as follows:
!
!  if (coefn .lt. 6./16.) coefn = 6./16.
!  coefn        = coefn/5.
!
!  Compute the number of available protons.
!-----------------------------------------------------------------------

zm20               = DMIN1( zh - 20.d+00, 8.d+00 )
IF ( zm20 <= zero ) THEN
  emitnc           = zero
  absrnc           = zero
  RETURN
END IF ! zm20 <= zero

tmev               = kmev * t
ronc               = xh * rho/( ah * rmu )

DO k = 1, nez
!-----------------------------------------------------------------------
!  Compute the electron energy necessary to produce an e-neutrino
!   of energy e_in.
!-----------------------------------------------------------------------

  eelec            = e_in(k) + cmpn - cmpp + edmpe + dmnp
  eelec2           = eelec * eelec
  IF ( eelec < zero  .or.  one - me2/eelec2 < zero ) THEN
    emitnc(k)        = zero
    absrnc(k)        = zero
    CYCLE
  END IF ! eelec < zero

!-----------------------------------------------------------------------
!  Compute the emission rate, pemit, per electron
!-----------------------------------------------------------------------

  ceelec           = DSQRT( DMAX1( one - me2/eelec2, epsilon ) )
  pemit            = coef * ronc * zm20 * eelec2 * ceelec

!-----------------------------------------------------------------------
!  Compute the electron occupation number, F_e
!-----------------------------------------------------------------------

  eta              = ( eelec - cmpe )/tmev
  F_e              = one/( one + fexp(eta) )

!-----------------------------------------------------------------------
!  Compute the emission inverse mean free path
!-----------------------------------------------------------------------

  emitnc(k)        = pemit * F_e * coefn

!-----------------------------------------------------------------------
!  Compute the absorption inverse mean free path using detailed
!   balance.
!-----------------------------------------------------------------------

  etam             = ( e_in(k) + dmnp + cmpn - cmpp - cmpe )/tmev
  IF ( eta < maxeta ) THEN
    absrnc(k)      = emitnc(k) * fexp(etam)
  ELSE
    absrnc(k)      = pemit * coefn * fexp( etam - eta )
  END IF ! eta < maxeta

END DO ! k = 1,nez

RETURN
END SUBROUTINE abemnc_weaklib
