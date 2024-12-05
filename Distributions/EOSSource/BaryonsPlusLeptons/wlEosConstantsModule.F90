
MODULE wlEosConstantsModule

USE wlKindModule, ONLY: dp

REAL(dp), PARAMETER  :: pi        = 3.1415926535897932385d0 ! pi

REAL(dp), PARAMETER  :: h_cgs     = 6.62607535359d-27 ! Planck's constant [g cm^2 s^{-1}]

REAL(dp), PARAMETER  :: hbar_cgs  = h_cgs/( 2.d0*pi ) ! Planck's constant divided by 2*pi [MeV s]

REAL(dp), PARAMETER  :: cvel      = 2.99792d+10   ! velocity of light [cm s^{-1}]

REAL(dp), PARAMETER  :: cvel_inv  = 1.d0/cvel     ! 1/( velocity of light [cm s^{-1}])

REAL(dp), PARAMETER  :: g         = 6.67384d-08   ! Gravitational constant [dynes cm^{2} g^{-2}]

REAL(dp), PARAMETER  :: ergmev    = 1.602177d-6   ! ergs per MeV

REAL(dp), PARAMETER  :: ergfoe    = 1.d-51        ! ergs/s per foe

REAL(dp), PARAMETER  :: kmev      = 8.61733d-11   ! Boltzmann's constant [MeV K^{-1}]

REAL(dp), PARAMETER  :: kev       = kmev*1.0d6   ! Boltzmann's constant [eV K^{-1}]

REAL(dp), PARAMETER  :: kmev_inv  = 1.d0/kmev     ! Boltzmann's constant inverse [K MeV^{-1}]

REAL(dp), PARAMETER  :: kerg      = 1.380648792741d-16

REAL(dp), PARAMETER  :: cm3fm3    = 1.d-39        ! cm**3/fm**3

REAL(dp), PARAMETER  :: avn       = 6.022141d+23  ! avogadro's number

REAL(dp), PARAMETER  :: rmu       = 1.0d0/avn     ! atomic mass unit

REAL(dp), PARAMETER  :: me        = 0.510998d+00  ! electron mass [MeV]

REAL(dp), PARAMETER  :: mmu       = 105.65837d+00 ! muon mass [MeV]

REAL(dp), PARAMETER  :: mp        = 938.27208816d+00   ! proton mass [MeV]

REAL(dp), PARAMETER  :: mn        = 939.56542052d+00   ! neutron mass [MeV]

REAL(dp), PARAMETER  :: mb        = 931.494d+00   ! baryon mass [MeV]

REAL(dp), PARAMETER  :: qe        = 4.8032042712d-10 ! electron charge [statcoulomb]

REAL(dp), PARAMETER  :: sigma_sb  = 5.670374419d-5 ! Stefan Bolzmann constant in csg [erg cm^{−2} s^{−1} K^{−4}]

REAL(dp), PARAMETER  :: asol      = 4.0d0 * sigma_sb / cvel

END MODULE wlEosConstantsModule