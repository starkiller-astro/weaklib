!-----------------------------------------------------------------------
!    Module:       parallel_module
!    Author:       S. W. Bruenn
!    Date:         2/25/03
!
!    Declares all variables needed for parallelization
!-----------------------------------------------------------------------

MODULE wlParallelModule

IMPLICIT none

!-----------------------------------------------------------------------
!  MPI Variables
!-----------------------------------------------------------------------

INTEGER                          :: ierr          ! initialization variable for MPI
INTEGER                          :: myid          ! rank of each processor (MPI)          
INTEGER                          :: myid_g        ! rank of each processor when split into grids
INTEGER                          :: myid_y        ! rank of each processor when split wrt k_block         
INTEGER                          :: myid_z        ! rank of each processor when split wrt j_block 
INTEGER                          :: my_grid       ! number of grid for each processor
INTEGER                          :: MPI_COMM_GRID ! communicator handle for local grid
INTEGER                          :: MPI_COMM_ROW  ! communicator handle for local xy-slab
INTEGER                          :: MPI_COMM_COL  ! communicator handle for local xz-slab
INTEGER                          :: j_block_col   ! j block (y coordinate in XZ transpose) on a given processor
INTEGER                          :: k_block_row   ! k block (z coordinate in XY transpose) on a given processor

!-----------------------------------------------------------------------
!  The Parallel Decomposition
!  For a Parallel job, the total number of processors is factorized into 
!  n_proc_y * n_proc_z.  The establishment of the MPI_COMM_ROW and MPI_COMM_COL
!  Communicators establishes a processor grid, to which bundles of rays are assigned.
!  For each sweep orientation, the dimensions and locations in this grid are given 
!  by different variables.
!  Orientation | X dimen / location |  Y dimen / location |  Z dimen / location
!     X               1  /  1         n_proc_y / myid_y     n_proc_z / myid_z
!     Y         n_proc_y / myid_y           1  /  1         n_proc_z / k_block_row
!     Z         n_proc_z / myid_z     n_proc_y / j_block_col      1  /  1            
!-----------------------------------------------------------------------

END MODULE wlParallelModule
