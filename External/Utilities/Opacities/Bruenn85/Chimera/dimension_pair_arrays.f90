SUBROUTINE dimension_pair_arrays( nx, nez, nnu, ij_ray_dim, ik_ray_dim )
!-----------------------------------------------------------------------
!
!    File:         dimension_pair_arrays
!    Module:       dimension_pair_arrays
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/10/04
!
!    Purpose:
!      To allocate dimensions to the pair production arrays.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  nx         : x (radial) array dimension
!  nez        : neutrino energy array dimension
!  nnu        : neutrino flavor array dimension
!  ij_ray_dim : number of y-zones on a processor before swapping
!  ik_ray_dim : number of z-zones on a processor before swapping
!
!    Output arguments:
!        none
!
!    Include files:
!  numerical_module
!  physcnst_module
!
!  edit_module
!  pair_module
!  parallel_module
!
!-----------------------------------------------------------------------

USE numerical_module, ONLY : zero
USE physcnst_module, ONLY: Gw, mp, hbar, cvel, pi

USE edit_module, ONLY : nlog
USE pair_module
USE parallel_module, ONLY : myid

IMPLICIT none

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: nx            ! radial array dimension
INTEGER, INTENT(in)              :: nez           ! neutrino energy array dimension
INTEGER, INTENT(in)              :: nnu           ! neutrino flavor array dimension
INTEGER, INTENT(in)              :: ij_ray_dim    ! number of y-zones on a processor before swapping
INTEGER, INTENT(in)              :: ik_ray_dim    ! number of z-zones on a processor before swapping

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=11)               :: var_name

INTEGER                          :: istat         ! allocation status

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 FORMAT (' Electron-positron pair annihilation arrays have been dimensioned')
 1001 FORMAT (' Allocation problem for array ',a10,' in dimension_pair_arrays')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!               \\\\\ ALLOCATE PAIR_MODULE ARRAYS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Cube indices
!-----------------------------------------------------------------------

ALLOCATE (idrpp(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'idrpp      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (itrpp(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'itrpp      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (iyrpp(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'iyrpp      '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Electron-positron pair annihilation functions at cube corners
!-----------------------------------------------------------------------

ALLOCATE (paira0i(nez,nez,2,2,2,nx,ij_ray_dim*ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'paira0i    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (paira0ii(nez,nez,2,2,2,nx,ij_ray_dim*ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'paira0ii   '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Interpolated electron-positron pair annihilation functions
!-----------------------------------------------------------------------

ALLOCATE (paf(6,nez,nnu,nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'paf        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (pafd(6,nez,nnu,nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'pafd       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (paft(6,nez,nnu,nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'paft       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (pafy(6,nez,nnu,nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'pafy       '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (pafp0(6,nez,nez,nnu,nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'pafp0      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (pafp1(6,nez,nez,nnu,nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'pafp1      '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  EOS change flags for table regeneration tests
!-----------------------------------------------------------------------

ALLOCATE (c_pair_prev(nx,ij_ray_dim,ik_ray_dim), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'c_pair_prev'; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!   Arrays for electron-positron edits
!-----------------------------------------------------------------------

ALLOCATE (rtapepr(6,nez,nnu,nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rtapepr    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (rmdntspr(3,nez,nnu,nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rmdntspr   '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!             \\\\\ INITIALIZE PAIR_MODULE ARRAYS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Cube indices
!-----------------------------------------------------------------------

idrpp                     = 0
itrpp                     = 0
iyrpp                     = 0

!-----------------------------------------------------------------------
!  Electron-positron pair annihilation functions at cube corners
!-----------------------------------------------------------------------

paira0i                   = -100.d0
paira0ii                  = -100.d0

!-----------------------------------------------------------------------
!  Interpolated electron-positron pair annihilation functions
!-----------------------------------------------------------------------

paf                       = zero
pafd                      = zero
paft                      = zero
pafy                      = zero
pafp0                     = zero
pafp1                     = zero

!-----------------------------------------------------------------------
!  EOS change flags for table regeneration tests
!-----------------------------------------------------------------------

c_pair_prev               = ' '

!-----------------------------------------------------------------------
!  Get quadrature points and weights and initialize constants
!-----------------------------------------------------------------------

CALL gquad( nleg, xe, wte, nleg )

!-----------------------------------------------------------------------
!        Initialize.
!
!        coef = (1/(2*pi)**2)*(1/!)*(1/(hb*!)**3)*g**2/pi
!
!        g2 has units of cm3/MeV2 s.
!        coef has units of 1/[ MeV5 cm ].
!-----------------------------------------------------------------------

g2               = ( Gw/mp**2 )**2 * hbar**2 * cvel**3 ! cm3/MeV2 s
coef             = 2.d0 * pi * ( 1.d0/cvel ) * ( 1.d0/( 2.d0 * pi * hbar * cvel ) )**3 * g2 / pi

IF ( myid == 0 ) WRITE (nlog,101)

RETURN
END SUBROUTINE dimension_pair_arrays
