!***********************************************************************
!
! Physical constants module
!
!***********************************************************************

MODULE wlExtPhysicalConstantsModule

USE wlKindModule, ONLY: dp
USE wlExtNumericalModule, ONLY: pi

REAL(dp), PARAMETER  :: Gw        = 1.027d-5      ! weak interaction constant (proton mass).

REAL(dp), PARAMETER  :: Gw_MeV    = 1.16636d-11   ! weak interaction constant [MeV].

REAL(dp), PARAMETER  :: sin2W     = 0.23116d+00   ! weak mixing angle.

REAL(dp), PARAMETER  :: gv        = 1.0d+00       ! Fermi beta decay constant. This is unity as
!                                                       the coupling for Fermi decay is absorbed 
!                                                       by the weak interaction constant.

REAL(dp), PARAMETER  :: ga        = 1.21d+00      ! Gamow-Teller beta decay constant. This is 
!                                                       not unity as the oupling for Gamow-Teller
!                                                       decay is renormalized. 

REAL(dp), PARAMETER  :: cv        = 0.96d+00      ! weak interaction constant, which in the 
!                                                       standard model is given by
!                                                               1       2
!                                                          cv = _ + 2sin (W)
!                                                               2
!                                                       where W is the 'Weinberg angle'.  

REAL(dp), PARAMETER  :: ca        = 0.5d+00       ! weak interaction constant, which in the 
!                                                       standard model is given by
!                                                          cv = 1/2 .

REAL(dp), PARAMETER  :: mu_p      = 1.79285d0     ! anomalous proton magnetic moment

REAL(dp), PARAMETER  :: mu_n      = -1.91304d0    ! anomalous neutron magnetic moment

REAL(dp), PARAMETER  :: cvel      = 2.99792d+10   ! velocity of light [cm s^{-1}]

REAL(dp), PARAMETER  :: cvel_inv  = 1.d0/cvel     ! 1/( velocity of light [cm s^{-1}])

REAL(dp), PARAMETER  :: h         = 4.13567d-21   ! Planck's constant [MeV s]

REAL(dp), PARAMETER  :: hbar      = h/( 2.d0*pi ) ! Planck's constant divided by 2*pi [MeV s]

REAL(dp), PARAMETER  :: hbarc     = 197.327d+00   ! hbar * c (GeV fm)

REAL(dp), PARAMETER  :: g         = 6.67384d-08   ! Gravitational constant [dynes cm^{2} g^{-2}]

REAL(dp), PARAMETER  :: ergmev    = 1.602177d-6   ! ergs per MeV

REAL(dp), PARAMETER  :: ergfoe    = 1.d-51        ! ergs/s per foe

REAL(dp), PARAMETER  :: kmev      = 8.61733d-11   ! Boltzmann's constant [MeV K^{-1}]

REAL(dp), PARAMETER  :: kmev_inv  = 1.d0/kmev     ! Boltzmann's constant inverse [K MeV^{-1}]

REAL(dp), PARAMETER  :: cm3fm3    = 1.d-39        ! cm**3/fm**3

REAL(dp), PARAMETER  :: avn       = 6.022141d+23  ! avogadro's number

REAL(dp), PARAMETER  :: rmu       = 1.0d0/avn     ! atomic mass unit

REAL(dp), PARAMETER  :: me        = 0.510998d+00  ! electron mass [MeV]

REAL(dp), PARAMETER  :: mp        = 938.272d+00   ! proton mass [MeV]

REAL(dp), PARAMETER  :: mn        = 939.565d+00   ! neutron mass [MeV]

REAL(dp), PARAMETER  :: mb        = 931.494d+00   ! baryon mass [MeV]

REAL(dp), PARAMETER  :: mp_ex     = 7.2889705d+00 ! proton mass excess [MeV]

REAL(dp), PARAMETER  :: mn_ex     = 8.0713171d+00 ! neutron mass excess [MeV]

REAL(dp), PARAMETER  :: dmnp      = 1.29333d+00   ! mn - mp

REAL(dp), PARAMETER  :: msolar    = 1.98855d+33   ! mass of the sun [g]

REAL(dp), PARAMETER  :: gamhv     = 2.5d+00       ! mass of the sun [g]

REAL(dp), PARAMETER  :: wnm       = - 16.0d+00    ! mass of the sun [g]

REAL(dp), PARAMETER  :: ws        = 31.5d+00      ! mass of the sun [g]

REAL(dp), PARAMETER  :: xk0       = 180.0d+00     ! mass of the sun [g]

REAL(dp), PARAMETER  :: xkzafac   = 2.0d+00       ! mass of the sun [g]

REAL(dp), PARAMETER  :: asig      = 8.56345d-8    ! radiation constant [MeV fm^{-3} MeV^{-4}]

REAL(dp), PARAMETER  :: kfm       = cm3fm3/rmu    ! # nucleons/gram )( cm3/fm3 ) [g^{-1} cm^{3} fm^{-3}]

REAL(dp), PARAMETER  :: kp        = ergmev/cm3fm3 ! ( erg/cm3 ) / ( mev/fm3 ) [erg cm^{-3} MeV^{-1} fm^{3}]

REAL(dp), PARAMETER  :: ku        = ergmev/rmu    ! ( # nucleons/gram )( erg/mev ) [g^{-1} erg MeV^{-1}]

REAL(dp), PARAMETER  :: cv_p      = 0.5d0 - 2.d0 * sin2W    ! vector current neutrino-proton coupling constant

REAL(dp), PARAMETER  :: cv_n      = - 0.5d0       ! vector current neutrino-neutron coupling constant

REAL(dp), PARAMETER  :: ca_p      = ga/2.d0       ! axial vector current neutrino-proton coupling constant

REAL(dp), PARAMETER  :: ca_n      = -ga/2.d0      ! axial vector current neutrino-neutron coupling constant


!-----------------------------------------------------------------------
! Constant combinations for weak magnetism corrections
!-----------------------------------------------------------------------

REAL(dp), PARAMETER  :: CV_p2     = CV_p * CV_p

REAL(dp), PARAMETER  :: CV_n2     = CV_n * CV_n

REAL(dp), PARAMETER  :: CA_p2     = CA_p * CA_p

REAL(dp), PARAMETER  :: CA_n2     = CA_n * CA_n

REAL(dp), PARAMETER  :: F2_p      = 0.5d0 * ( mu_p - mu_n ) - 2.d0 * sin2W * mu_p    ! proton form factor

REAL(dp), PARAMETER  :: F2_n      = -0.5d0 * ( mu_p - mu_n ) - 2.d0 * sin2W * mu_n    ! neutron form factor

REAL(dp), PARAMETER  :: CV_F2CA_p = (CV_p + F2_p) * CA_p

REAL(dp), PARAMETER  :: CV_F2CA_n = (CV_n + F2_n) * CA_n

REAL(dp), PARAMETER  :: CVF2_p    = CV_p * F2_p

REAL(dp), PARAMETER  :: CVF2_n    = CV_n * F2_n

REAL(dp), PARAMETER  :: F22_p     = F2_p * F2_p

REAL(dp), PARAMETER  :: F22_n     = F2_n * F2_n

REAL(dp), PARAMETER  :: CV5CA_p   = ( 2.d0/3.d0 ) * ( CV_p2 + 5.d0 * CA_p2 )

REAL(dp), PARAMETER  :: CV5CA_n   = ( 2.d0/3.d0 ) * ( CV_n2 + 5.d0 * CA_n2 )

!-----------------------------------------------------------------------
!  Coefficients for neutrino number and energy densities, and fluxes
!-----------------------------------------------------------------------

REAL(dp), PARAMETER  ::  ncoef    = 4.d+00 * pi/( h * cvel )**3
REAL(dp), PARAMETER  ::  ecoef    = 4.d+00 * pi * ergmev/( h * cvel )**3
REAL(dp), PARAMETER  ::  csqinv   = 1.d0/cvel**2


!-----------------------------------------------------------------------
! Constant combinations for BoltzTran
!-----------------------------------------------------------------------
REAL(dp), PARAMETER  :: gf      = 8.957d-44   ! Fermi's constant [MeV cm3].
REAL(dp), PARAMETER  :: hbarc_MeV   = 1.97326946d-11 ! [MeV cm]
REAL(dp), PARAMETER  :: mbg         = mb*ergmev*cvel_inv*cvel_inv ! baryon mass [g]
REAL(dp), PARAMETER  :: stherm  = 1.0d0/(pi*hbarc_MeV**4)
REAL(dp), PARAMETER  :: therm1  = gf*gf*(gv*gv+3.d0*ga*ga)*stherm
REAL(dp), PARAMETER  :: therm2  = 2.d0*gf*gf*ga*ga*stherm/7.0d0
REAL(dp), PARAMETER  :: therm3  = gf*gf*stherm

END MODULE wlExtPhysicalConstantsModule

