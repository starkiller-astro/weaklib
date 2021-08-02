!***********************************************************************
!
! physcnst_module
!
!***********************************************************************

MODULE physcnst_module

USE kind_module, ONLY: double
USE numerical_module, ONLY: pi

REAL(double), PARAMETER  :: Gw        = 1.027d-5      ! weak interaction constant (proton mass).

REAL(double), PARAMETER  :: Gw_MeV    = 1.16636d-11   ! weak interaction constant [MeV].

REAL(double), PARAMETER  :: sin2W     = 0.23116d+00   ! weak mixing angle.

REAL(double), PARAMETER  :: gv        = 1.0d+00       ! Fermi beta decay constant. This is unity as
!                                                       the coupling for Fermi decay is absorbed 
!                                                       by the weak interaction constant.

REAL(double), PARAMETER  :: ga        = 1.21d+00      ! Gamow-Teller beta decay constant. This is 
!                                                       not unity as the oupling for Gamow-Teller
!                                                       decay is renormalized. 

REAL(double), PARAMETER  :: cv        = 0.96d+00      ! weak interaction constant, which in the 
!                                                       standard model is given by
!                                                               1       2
!                                                          cv = _ + 2sin (W)
!                                                               2
!                                                       where W is the 'Weinberg angle'.  

REAL(double), PARAMETER  :: ca        = 0.5d+00       ! weak interaction constant, which in the 
!                                                       standard model is given by
!                                                          cv = 1/2 .

REAL(double), PARAMETER  :: mu_p      = 1.79285d0     ! anomalous proton magnetic moment

REAL(double), PARAMETER  :: mu_n      = -1.91304d0    ! anomalous neutron magnetic moment

REAL(double), PARAMETER  :: cvel      = 2.99792d+10   ! velocity of light [cm s^{-1}]

REAL(double), PARAMETER  :: cvel_inv  = 1.d0/cvel     ! 1/( velocity of light [cm s^{-1}])

REAL(double), PARAMETER  :: h         = 4.13567d-21   ! Planck's constant [MeV s]

REAL(double), PARAMETER  :: hbar      = h/( 2.d0*pi ) ! Planck's constant divided by 2*pi [MeV s]

REAL(double), PARAMETER  :: hbarc     = 197.327d+00   ! hbar * c (GeV fm)

REAL(double), PARAMETER  :: g         = 6.67384d-08   ! Gravitational constant [dynes cm^{2} g^{-2}]

REAL(double), PARAMETER  :: ergmev    = 1.602177d-6   ! ergs per MeV

REAL(double), PARAMETER  :: ergfoe    = 1.d-51        ! foe per erg

REAL(double), PARAMETER  :: kmev      = 8.61733d-11   ! Boltzmann's constant [MeV K^{-1}]

REAL(double), PARAMETER  :: kmev_inv  = 1.d0/kmev     ! Boltzmann's constant inverse [K MeV^{-1}]

REAL(double), PARAMETER  :: cm3fm3    = 1.d-39        ! cm**3/fm**3

REAL(double), PARAMETER  :: avn       = 6.022141d+23  ! avogadro's number

REAL(double), PARAMETER  :: rmu       = 1.0d0/avn     ! mean baryon mass [g]

REAL(double), PARAMETER  :: me        = 0.510998d+00  ! electron mass [MeV]

REAL(double), PARAMETER  :: mp        = 938.272d+00   ! proton mass [MeV]

REAL(double), PARAMETER  :: mn        = 939.565d+00   ! neutron mass [MeV]

REAL(double), PARAMETER  :: mb        = 931.494d+00   ! baryon mass [MeV]

REAL(double), PARAMETER  :: mp_ex     = 7.2889705d+00 ! proton mass excess [MeV]

REAL(double), PARAMETER  :: mn_ex     = 8.0713171d+00 ! neutron mass excess [MeV]

REAL(double), PARAMETER  :: dmnp      = 1.29333d+00   ! mn - mp

REAL(double), PARAMETER  :: msolar    = 1.98855d+33   ! mass of the sun [g]

REAL(double), PARAMETER  :: gamhv     = 2.5d+00       ! mass of the sun [g]

REAL(double), PARAMETER  :: wnm       = - 16.0d+00    ! mass of the sun [g]

REAL(double), PARAMETER  :: ws        = 31.5d+00      ! mass of the sun [g]

REAL(double), PARAMETER  :: xk0       = 180.0d+00     ! mass of the sun [g]

REAL(double), PARAMETER  :: xkzafac   = 2.0d+00       ! mass of the sun [g]

REAL(double), PARAMETER  :: asig      = 8.56345d-8    ! radiation constant [MeV fm^{-3} MeV^{-4}]

REAL(double), PARAMETER  :: kfm       = cm3fm3/rmu    ! # nucleons/gram )( cm3/fm3 ) [g^{-1} cm^{3} fm^{-3}]

REAL(double), PARAMETER  :: kp        = ergmev/cm3fm3 ! ( erg/cm3 ) / ( mev/fm3 ) [erg cm^{-3} MeV^{-1} fm^{3}]

REAL(double), PARAMETER  :: ku        = ergmev/rmu    ! ( # nucleons/gram )( erg/mev ) [g^{-1} erg MeV^{-1}]

REAL(double), PARAMETER  :: cv_p      = 0.5d0 - 2.d0 * sin2W    ! vector current neutrino-proton coupling constant

REAL(double), PARAMETER  :: cv_n      = - 0.5d0       ! vector current neutrino-neutron coupling constant

REAL(double), PARAMETER  :: ca_p      = ga/2.d0       ! axial vector current neutrino-proton coupling constant

REAL(double), PARAMETER  :: ca_n      = -ga/2.d0      ! axial vector current neutrino-neutron coupling constant


!-----------------------------------------------------------------------
! Constant combinations for weak magnetism corrections
!-----------------------------------------------------------------------

REAL(double), PARAMETER  :: CV_p2     = CV_p * CV_p

REAL(double), PARAMETER  :: CV_n2     = CV_n * CV_n

REAL(double), PARAMETER  :: CA_p2     = CA_p * CA_p

REAL(double), PARAMETER  :: CA_n2     = CA_n * CA_n

REAL(double), PARAMETER  :: F2_p      = 0.5d0 * ( mu_p - mu_n ) - 2.d0 * sin2W * mu_p    ! proton form factor

REAL(double), PARAMETER  :: F2_n      = -0.5d0 * ( mu_p - mu_n ) - 2.d0 * sin2W * mu_n    ! neutron form factor

REAL(double), PARAMETER  :: CV_F2CA_p = (CV_p + F2_p) * CA_p

REAL(double), PARAMETER  :: CV_F2CA_n = (CV_n + F2_n) * CA_n

REAL(double), PARAMETER  :: CVF2_p    = CV_p * F2_p

REAL(double), PARAMETER  :: CVF2_n    = CV_n * F2_n

REAL(double), PARAMETER  :: F22_p     = F2_p * F2_p

REAL(double), PARAMETER  :: F22_n     = F2_n * F2_n

REAL(double), PARAMETER  :: CV5CA_p   = ( 2.d0/3.d0 ) * ( CV_p2 + 5.d0 * CA_p2 )

REAL(double), PARAMETER  :: CV5CA_n   = ( 2.d0/3.d0 ) * ( CV_n2 + 5.d0 * CA_n2 )

!-----------------------------------------------------------------------
!  Coefficients for neutrino number and energy densities, and fluxes
!-----------------------------------------------------------------------

REAL(double), PARAMETER  ::  ncoef    = 4.d+00 * pi/( h * cvel )**3
REAL(double), PARAMETER  ::  ecoef    = 4.d+00 * pi * ergmev/( h * cvel )**3
REAL(double), PARAMETER  ::  csqinv   = 1.d0/cvel**2

END MODULE physcnst_module

