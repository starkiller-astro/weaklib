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
  USE wlIOModuleHDF, ONLY: InitializeHDF, FinalizeHDF
  USE wlEOSIOModuleHDF, ONLY: ReadEquationOfStateTableHDF,         & 
                              WriteEquationOfStateTableHDF,        &
                              MatchTableStructure,                 &
                              BroadcastEquationOfStateTableParallel

  USE MPI

  implicit none

  INTEGER  :: i, num_procs
  TYPE(EquationOfStateTableType) :: EOSTable
  TYPE(EquationOfStateTableType) :: OldEOSTable
  TYPE(DVIDType)                 :: NewDVID
  INTEGER                        :: NewnVariables

!-----------------------------------------------------------------------
!  Initialize MPI
!-----------------------------------------------------------------------

  CALL MPI_INIT( ierr )
  CALL MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
  CALL MPI_COMM_SIZE( MPI_COMM_WORLD, num_procs, ierr )

IF ( myid == 0 ) THEN
  CALL InitializeHDF( )

  CALL ReadEquationOfStateTableHDF( OldEOSTable, 'EquationOfStateTable.h5' )

  NewnVariables = 14

  NewDVID % iPressure = 1
  NewDVID % iEntropyPerBaryon = 3
  NewDVID % iInternalEnergyDensity = 2
  NewDVID % iElectronChemicalPotential = 6
  NewDVID % iProtonChemicalPotential = 5
  NewDVID % iNeutronChemicalPotential = 4
  NewDVID % iProtonMassFraction = 8
  NewDVID % iNeutronMassFraction = 7
  NewDVID % iAlphaMassFraction = 0
  NewDVID % iHeavyMassFraction = 9
  NewDVID % iHeavyChargeNumber = 11
  NewDVID % iHeavyMassNumber = 10
  NewDVID % iHeavyBindingEnergy = 13
  NewDVID % iThermalEnergy = 14
  NewDVID % iGamma1 = 12
  
  CALL MatchTableStructure( OldEOSTable, EOSTable, NewDVID, NewnVariables )
  
  WRITE(*,*) "Table Restructured"
  
END IF

  CALL BroadcastEquationOfStateTableParallel( EOSTable, 0, myid, ierr, MPI_COMM_WORLD )

IF ( myid == 0 ) THEN
  CALL DeAllocateEquationOfStateTable( OldEOSTable )
  CALL FinalizeHDF( )
END IF

END PROGRAM wlParallelReadTest
