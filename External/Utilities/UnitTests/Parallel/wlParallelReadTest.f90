PROGRAM wlParallelReadTest

!  The goal of this test is to read an instance of the Equation of State table
!  and then to broadcast it to all other parallel processes using MPI. 
!  
!--------------------------------------------------------------------------------

  USE wlKindModule, ONLY: dp 
  USE HDF5
  USE wlExtEOSWrapperModule, ONLY: wlGetFullEOS
  USE wlEquationOfStateTableModule
  USE wlInterpolationModule
  USE wlExtTableFinishingModule
  USE wlParallelModule, ONLY: myid, ierr
  USE wlExtNumericalModule, ONLY: epsilon, zero
  USE wlIOModuleHDF, ONLY: InitializeHDF, FinalizeHDF,          & 
                           ReadEquationOfStateTableHDF,         & 
                           WriteEquationOfStateTableHDF,        &
                           ReadCHIMERAHDF,                      &
                           ReadEquationOfStateTableParallelHDF 

  USE MPI

  implicit none

  INTEGER  :: i, num_procs
  TYPE(EquationOfStateTableType) :: EOSTable

!-----------------------------------------------------------------------
!  Initialize MPI
!-----------------------------------------------------------------------

  CALL MPI_INIT( ierr )
  CALL MPI_COMM_RANK( MPI_COMM_WORLD, myid , ierr )
  CALL MPI_COMM_SIZE( MPI_COMM_WORLD, num_procs, ierr )

  IF ( myid == 0 ) CALL InitializeHDF( )

  CALL ReadEquationOfStateTableParallelHDF( EOSTable,                                    &
                                            "StandardResEquationOfStateTable7-31-15.h5", &
                                            0, MPI_COMM_WORLD )

  CALL DeAllocateEquationOfStateTable( EOSTable )

  IF ( myid == 0 ) CALL FinalizeHDF( )

END PROGRAM wlParallelReadTest
