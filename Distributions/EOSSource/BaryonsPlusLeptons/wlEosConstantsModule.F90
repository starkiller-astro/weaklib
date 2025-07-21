
MODULE wlEosConstantsModule

USE wlKindModule, ONLY: dp

REAL(dp), PARAMETER  :: pi        = 3.1415926535897932385d0 ! pi

REAL(dp), PARAMETER  :: ergmev    = 1.602177d-6   ! ergs per MeV

REAL(dp), PARAMETER  :: h_cgs     = 6.62607535359d-27 ! Planck's constant [g cm^2 s^{-1}]

REAL(dp), PARAMETER  :: hbar_cgs  = h_cgs/( 2.d0*pi ) ! Planck's constant divided by 2*pi [MeV s]

REAL(dp), PARAMETER  :: h         = h_cgs / ergmev   ! Planck's constant [MeV s]

REAL(dp), PARAMETER  :: hbar      = h/( 2.d0*pi ) ! Planck's constant divided by 2*pi [MeV s]

REAL(dp), PARAMETER  :: cvel      = 2.99792d+10   ! velocity of light [cm s^{-1}]

REAL(dp), PARAMETER  :: hbarc     = hbar*cvel     ! hbar * c (MeV cm)

REAL(dp), PARAMETER  :: cvel_inv  = 1.d0/cvel     ! 1/( velocity of light [cm s^{-1}])

REAL(dp), PARAMETER  :: g         = 6.67384d-08   ! Gravitational constant [dynes cm^{2} g^{-2}]

REAL(dp), PARAMETER  :: ergfoe    = 1.d-51        ! ergs/s per foe

REAL(dp), PARAMETER  :: kmev      = 8.61733d-11   ! Boltzmann's constant [MeV K^{-1}]

REAL(dp), PARAMETER  :: kev       = kmev*1.0d6    ! Boltzmann's constant [eV K^{-1}]

REAL(dp), PARAMETER  :: kmev_inv  = 1.d0/kmev     ! Boltzmann's constant inverse [K MeV^{-1}]

REAL(dp), PARAMETER  :: kerg      = 1.380648792741d-16

REAL(dp), PARAMETER  :: cm3fm3    = 1.d-39        ! cm**3/fm**3

REAL(dp), PARAMETER  :: avn       = 6.022141d+23  ! avogadro's number

REAL(dp), PARAMETER  :: rmu       = 1.0d0/avn     ! atomic mass unit

REAL(dp), PARAMETER  :: me        = 0.510998d+00  ! electron mass [MeV]

REAL(dp), PARAMETER  :: mmu       = 105.65837d+00 ! muon mass [MeV]

REAL(dp), PARAMETER  :: mp        = 938.27208816d+00   ! proton mass [MeV]

REAL(dp), PARAMETER  :: mn        = 939.56542052d+00   ! neutron mass [MeV]

REAL(dp), PARAMETER  :: mpi       = 139.57d0+00   ! charged pion mass [MeV]

REAL(dp), PARAMETER  :: mb        = 931.494d+00   ! baryon mass [MeV]

REAL(dp), PARAMETER  :: qe        = 4.8032042712d-10 ! electron charge [statcoulomb]

REAL(dp), PARAMETER  :: sigma_sb  = 5.670374419d-5 ! Stefan Bolzmann constant in csg [erg cm^{−2} s^{−1} K^{−4}]

REAL(dp), PARAMETER  :: asol      = 4.0d0 * sigma_sb / cvel

! Weak Interactions Parameters
REAL(dp), PARAMETER  :: Gw        = 1.027d-5      ! weak interaction constant (proton mass).

REAL(dp), PARAMETER  :: Gw_MeV    = 1.16636d-11   ! weak interaction constant [MeV].

REAL(dp), PARAMETER  :: sin2W     = 0.23147d+00   ! weak mixing angle.
!                                                   https://link.springer.com/article/10.1007/JHEP12(2024)026

REAL(dp), PARAMETER  :: gamma_p   =  2.79284734d0 ! proton magnetic moment from PDG 2024

REAL(dp), PARAMETER  :: gamma_n   = -1.9130427d0  ! neutron magnetic moment from PDG 2024

REAL(dp), PARAMETER  :: massA     = 1014.0d0      ! MeV. dipole axial mass. Poorly constrained
!                                                   https://www.nature.com/articles/s41586-022-05478-3

REAL(dp), PARAMETER  :: massV     = 840.0d0       ! MeV. dipole vector mass. Very robustly estimated
!                                                   https://link.springer.com/article/10.1007/JHEP08(2024)187

REAL(dp), PARAMETER  :: Vud       = 0.97367d0     ! ud element of the CKM matrix from PDG 2024

REAL(dp), PARAMETER  :: gv        = 1.0d+00       ! Fermi beta decay constant. This is unity as
!                                                       the coupling for Fermi decay is absorbed 
!                                                       by the weak interaction constant.

REAL(dp), PARAMETER  :: ga        = 1.2723d+00    ! Gamow-Teller beta decay constant. This is 
!                                                   not unity as the oupling for Gamow-Teller
!                                                   decay is renormalized. 

REAL(dp), PARAMETER  :: cv        = 0.96d+00      ! weak interaction constant, which in the 
!                                                       standard model is given by
!                                                               1       2
!                                                          cv = _ + 2sin (W)
!                                                               2
!                                                       where W is the 'Weinberg angle'.  

REAL(dp), PARAMETER  :: ca        = 0.5d+00       ! weak interaction constant, which in the 
!                                                       standard model is given by
!                                                          cv = 1/2 .

END MODULE wlEosConstantsModule