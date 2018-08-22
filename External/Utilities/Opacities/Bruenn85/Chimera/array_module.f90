!-----------------------------------------------------------------------
!    Module:       array_module
!    Author:       S. W. Bruenn
!    Date:         8/16/02
!
!    Specification of the array dimensions
!-----------------------------------------------------------------------

MODULE array_module

!-----------------------------------------------------------------------
!  nx, ny, nz, nez, nnu, nnc
!-----------------------------------------------------------------------
!  nx    : is the maximum number of radial zones
!  ny    : is the number of y-zones
!  nz    : is the number of z-zones
!  nez   : is the number of energy zones
!  nezext: is the number of energy zones added for interpolation
!  nnu   : is the number of neutrino species
!  nnc   : the maximum number of nuclei for matter not in nuclear statistical
!   equilibrium
!  n_eos : number of equations of state
!  ngrid : number of grids (1 for "unigrid", 2 for "baseball" )
!-----------------------------------------------------------------------

INTEGER                    :: nx
INTEGER                    :: ny
INTEGER                    :: nz
INTEGER                    :: nez
INTEGER                    :: nezext
INTEGER                    :: nnu
INTEGER                    :: nnc
INTEGER                    :: n_eos
INTEGER                    :: ngrid

!-----------------------------------------------------------------------
!  Tracer particles
!-----------------------------------------------------------------------
!  n_part : tracer particle array extent
!-----------------------------------------------------------------------

INTEGER                    :: n_part

!-----------------------------------------------------------------------
!  n_proc
!-----------------------------------------------------------------------
!  n_proc     : the number of processors assigned to the run
!  n_proc_y   : the number of processors assigned to the y-zones
!  n_proc_z   : the number of processors assigned to the z-zones
!  ij_ray_dim : the number of y-zones on a processor before swapping with y
!  ik_ray_dim : the number of z-zones on a processor before swapping with z
!  j_ray_dim  : the number of radial zones on a processor after swapping
!                with y
!  k_ray_dim  : the number of radial zones on a processor after swapping
!                with z
!-----------------------------------------------------------------------

INTEGER                    :: n_proc
INTEGER                    :: n_proc_y
INTEGER                    :: n_proc_z
INTEGER                    :: ij_ray_dim
INTEGER                    :: ik_ray_dim
INTEGER                    :: j_ray_dim
INTEGER                    :: k_ray_dim

!-----------------------------------------------------------------------
!  pool_size  : dimension of opacity pools
!-----------------------------------------------------------------------

INTEGER                    :: pool_size

END module array_module
