!-----------------------------------------------------------------------
!    Author:       Ran Chu
!    Date:         11/28/2019
!-----------------------------------------------------------------------

MODULE pair_module

USE kind_module
USE numerical_module, ONLY : zero
USE physcnst_module, ONLY: Gw, mp, hbar, cvel, pi

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
