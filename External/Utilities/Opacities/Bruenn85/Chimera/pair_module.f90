!-----------------------------------------------------------------------
!    Module:       pair_module
!    Author:       S. W. Bruenn
!    Date:         8/16/02
!-----------------------------------------------------------------------

MODULE pair_module

USE kind_module
USE numerical_module, ONLY : zero
USE physcnst_module, ONLY: Gw, mp, hbar, cvel, pi

USE edit_module, ONLY : nlog
USE parallel_module, ONLY : myid

!-----------------------------------------------------------------------
!  Cube indices
!-----------------------------------------------------------------------
!  idrpp(j), itrpp(j), and iyrpp(j) are integers defining the location of
!   log(rho), log(t), and ye for radial zone j on the grid for the
!   electron-positron pair annihilation process, i.e.,
!
!     idrpp(j,ij_ray,ik_ray)/dgrid < log(rho(j)) < ( idrpp(j,ij_ray,ik_ray) + 1 )/dgrid 
!
!     itrpp(j,ij_ray,ik_ray)/tgrid <  log(t(j))  < ( itrpp(j,ij_ray,ik_ray) + 1 )/tgrid
!
!     0.5 - iyrpp(j,ij_ray,ik_ray)/ygrid < ye < 0.5 - ( iyrpp(j,ij_ray,ik_ray) + 1 )/ygrid
!
!  The eight grid points surrounding log(rho), log(t), and ye for radial
!   zone j are referred to as the unit cube j. Electron-positron pair
!   annihilation rates for radial zone j are stored at the corners of
!   unit cube j. Rates for zone j are interpolated from the rates stored
!   at the corners.
!-----------------------------------------------------------------------

INTEGER, ALLOCATABLE, DIMENSION(:,:,:)                     :: idrpp
INTEGER, ALLOCATABLE, DIMENSION(:,:,:)                     :: itrpp
INTEGER, ALLOCATABLE, DIMENSION(:,:,:)                     :: iyrpp

!-----------------------------------------------------------------------
!  Electron-positron pair annihilation functions at cube corners
!-----------------------------------------------------------------------
!  paira0i(kp,k,id,it,iy,j,ij_ray,ik_ray),
!  paira0ii(kp,k,id,it,iy,j,ij_ray,ik_ray) : arrays of the zero
!   moments of the  electron-positron pair annihilation functions for
!   neutrinos of energy k, antineutrinos of energy kp, in radial zone
!   j at the unit cube corners id, it, and iy (id, it, iy = 1,2). 
!   Interpolations of the electron-positron pair annihilation functions
!   are computed from this table.
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:,:,:,:) :: paira0i
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:,:,:,:) :: paira0ii

!-----------------------------------------------------------------------
!  Interpolated electron-positron pair annihilation functions
!-----------------------------------------------------------------------
!  paf(i,j,k,n)        : electron-positron pair annihilation function i
!
!     i = 1: a0w, i = 2: b0w, i = 3: c0w,
!     i = 4: a1w, i = 5: b1w, i = 6: c1w.
!
!  pafd(i,j,k,n)     : d(paf(i,j,k,n))/d(density)
!  paft(i,j,k,n)     : d(paf(i,j,k,n))/d(temperature)
!  pafy(i,j,k,n)     : d(paf(i,j,k,n))/d(electron fraction)
!  pafp0(i,kp,j,k,n) : d(paf(i,j,k,n))/d(psi0(j,kp,n))
!  pafp1(i,kp,j,k,n) : d(paf(i,j,k,n))/d(psi1(j,kp,n))
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)         :: paf
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)         :: pafd
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)         :: paft
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)         :: pafy
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:,:)       :: pafp0
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:,:)       :: pafp1

!-----------------------------------------------------------------------
!  EOS change flag array for table regeneration tests
!-----------------------------------------------------------------------
!  c_pair_prev(j,i_ray) : array for storing c_eos_arry for comparison
!   comparison with c_eos_arry in the next cycle.
!-----------------------------------------------------------------------

CHARACTER (LEN=1), ALLOCATABLE, DIMENSION(:,:,:)           :: c_pair_prev

!-----------------------------------------------------------------------
!  Edit information for editng
!-----------------------------------------------------------------------
! rtapepr(i,k,n,j)     : coefficients for computing mean pair rates
!                      i=(1): artpe, (2): brtpe (3): crtpe
!                        (4): artae, (5): brtae (6): crtae
! rmdntspr(i,k,n,j)  : moments of m.f.p.
!                        i=(1): rmdnts, (2): rmdnts0, (3): rmdnts1
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)         :: rtapepr
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)         :: rmdntspr

!-----------------------------------------------------------------------
!  Coefficients
!-----------------------------------------------------------------------
! nnud     : upper limit to the neutrino flavor extent
! cpair1   : neutrino flavor coeeficient for computing pair rates
! cpair2   : neutrino flavor coeeficient for computing pair rates
!-----------------------------------------------------------------------

INTEGER, PARAMETER                                         :: nnud=4

REAL(KIND=double), DIMENSION(nnud)                         :: cpair1
REAL(KIND=double), DIMENSION(nnud)                         :: cpair2

!-----------------------------------------------------------------------
!   Gauss-Legendre parameters for Pairs
!-----------------------------------------------------------------------
! nleg     : number of points of Gauss-Lagendre quadrature
! xe       : points of Gauss-Lagendre quadrature
! wte      : weights of Gauss-Lagendre quadrature
!-----------------------------------------------------------------------

INTEGER, PARAMETER                                         :: nleg = 24

REAL(KIND=double), DIMENSION(nleg)                         :: xe
REAL(KIND=double), DIMENSION(nleg)                         :: wte

!-----------------------------------------------------------------------
!   Constants for paircal
!-----------------------------------------------------------------------
! g2       : combination of physical constants
! coef     : combination of physical constants
!-----------------------------------------------------------------------

REAL(KIND=double)    :: g2 = ( Gw/mp**2 )**2 * hbar**2 * cvel**3 ! cm3/MeV2 s
REAL(KIND=double)    :: coef = 2.d0 * pi * ( 1.d0/cvel ) * ( 1.d0/( 2.d0 * pi * hbar * cvel ) )**3 &
                               * ( Gw/mp**2 )**2 * hbar**2 * cvel**3 / pi


END module pair_module
