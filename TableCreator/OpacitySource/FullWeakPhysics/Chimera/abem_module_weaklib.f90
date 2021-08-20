MODULE abem_module_weaklib

USE wlKindModule, ONLY: dp


!-----------------------------------------------------------------------
!  Energy and angle quadrature points and weights
!-----------------------------------------------------------------------

INTEGER, PARAMETER           :: nleg_a = 64   ! number of points of ouitgoing lepton angular Gauss-Lagendre quadrature
REAL(dp), DIMENSION (nleg_a) :: x_a           ! scaled points of angular quadrature
REAL(dp), DIMENSION (nleg_a) :: wt_a          ! scaled weights of angular quadrature

INTEGER, PARAMETER           :: nleg_e = 64   ! number of points of ouitgoing lepton energy Gauss-Lagendre quadrature
REAL(dp), DIMENSION (nleg_e) :: x_e           ! scaled points of energy quadrature
REAL(dp), DIMENSION (nleg_e) :: wt_e          ! scaled weights of energy quadrature

END module abem_module_weaklib

