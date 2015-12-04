MODULE wlEOSIOModuleHDF

  USE wlKindModule, ONLY: dp
  USE HDF5 
  USE wlThermoStateModule
  USE wlDependentVariablesModule
  USE wlEquationOfStateTableModule
  USE wlIOModuleHDF

  implicit none
  PRIVATE
  INTEGER                                     :: hdferr

  PUBLIC WriteEquationOfStateTableHDF
  PUBLIC ReadEquationOfStateTableHDF
  PUBLIC DescribeEquationOfStateTable
  PUBLIC ReadEquationOfStateTableParallelHDF
  PUBLIC MatchTableStructure
  PUBLIC TransferDependentVariables

CONTAINS

  SUBROUTINE WriteEquationOfStateTableHDF( EOSTable )

    TYPE(EquationOfStateTableType), INTENT(inout) :: EOSTable
    
    INTEGER(HID_T)                                :: file_id
    INTEGER(HID_T)                                :: group_id
    
    CALL OpenFileHDF( "EquationOfStateTable.h5", .true., file_id )

    CALL OpenGroupHDF( "ThermoState", .true., file_id, group_id )
    CALL WriteThermoStateHDF( EOSTable % TS, group_id )
    CALL CloseGroupHDF( group_id )

    CALL OpenGroupHDF( "DependentVariables", .true., file_id, group_id )
    CALL WriteDependentVariablesHDF( EOSTable % DV, group_id )
    CALL CloseGroupHDF( group_id )

    CALL CloseFileHDF( file_id )

  END SUBROUTINE WriteEquationOfStateTableHDF

  SUBROUTINE DescribeEquationOfStateTable( EOSTable )

    TYPE(EquationOfStateTableType), INTENT(inout) :: EOSTable
    INTEGER :: i

    DO i = 1, 3
      WRITE (*,*) "Independent Variable", i, "=", &
                  EOSTable % TS % Names(i), "Units:", EOSTable % TS % Units(i)
      WRITE (*,*) "Independent Variable", i, " Minimum = " , EOSTable % TS % minValues(i)
      WRITE (*,*) "Independent Variable", i, "Maximum =" , EOSTable % TS % maxValues(i)
      WRITE (*,*) "Number of Independent Variable", i, "Points =" , EOSTable % TS % nPoints(i) 
      IF ( EOSTable % TS % LogInterp(i) == 1 ) THEN
         WRITE (*,*) "Independent Variable ", i, "Grid Logarithmically Spaced" 
         ELSE
         WRITE (*,*) "Independent Variable", i, "Grid Linearly Spaced" 
      END IF  
    END DO

    DO i = 1, EOSTable % nVariables
      WRITE (*,*) "Dependent Variable", i, "=", &
                  EOSTable % DV % Names(i), "Units:", EOSTable % DV % Units(i)
    END DO

  END SUBROUTINE DescribeEquationOfStateTable
  
  SUBROUTINE ReadEquationOfStateTableHDF( EOSTable, FileName )

    CHARACTER(len=*), INTENT(in)                  :: FileName
    INTEGER, DIMENSION(3)                         :: nPoints
    INTEGER                                       :: nVariables
    TYPE(EquationOfStateTableType), INTENT(inout) :: EOSTable
    INTEGER(HID_T)                                :: file_id
    INTEGER(HID_T)                                :: group_id

    CALL OpenFileHDF( FileName, .false., file_id )

    CALL OpenGroupHDF( "DependentVariables", .false., file_id, group_id )

    CALL ReadDimensionsHDF( nPoints, group_id )
    CALL ReadNumberVariablesHDF( nVariables, group_id )
    CALL CloseGroupHDF( group_id )

    CALL AllocateEquationOfStateTable( EOSTable, nPoints , nVariables )

    CALL ReadThermoStateHDF( EOSTable % TS, file_id )

    CALL ReadDependentVariablesHDF( EOSTable % DV, file_id )

    CALL DescribeEquationOfStateTable( EOSTable )

    CALL CloseFileHDF( file_id )

  END SUBROUTINE ReadEquationOfStateTableHDF

  SUBROUTINE ReadEquationOfStateTableParallelHDF( EOSTable, FileName, rootproc, COMMUNICATOR )

    USE MPI
    
    implicit none

    CHARACTER(len=*), INTENT(in)                  :: FileName
    INTEGER, INTENT(in)                           :: rootproc
    INTEGER, INTENT(in)                           :: COMMUNICATOR
    TYPE(EquationOfStateTableType), INTENT(inout) :: EOSTable
    INTEGER, DIMENSION(3)                         :: nPoints
    INTEGER, DIMENSION(4)                         :: buffer
    INTEGER                                       :: nStates, nVariables, i
    INTEGER                                       :: i_count
    INTEGER                                       :: myid, ierr 
    INTEGER(HID_T)                                :: file_id
    INTEGER(HID_T)                                :: group_id

    CALL MPI_COMM_RANK( COMMUNICATOR, myid , ierr )

    IF ( myid == rootproc ) THEN

      CALL OpenFileHDF( FileName, .false., file_id )

      CALL OpenGroupHDF( "DependentVariables", .false., file_id, group_id )

      CALL ReadDimensionsHDF( nPoints, group_id )
      CALL ReadNumberVariablesHDF( nVariables, group_id )
      CALL CloseGroupHDF( group_id )
    
      buffer(1) = nPoints(1)
      buffer(2) = nPoints(2)
      buffer(3) = nPoints(3)
      buffer(4) = nVariables

  WRITE (*,*) "in: process", myid, "buffer(4)", buffer(4) 

    END IF

    i_count = SIZE(buffer) 

    CALL MPI_BCAST( buffer, i_count, MPI_INTEGER, rootproc, COMMUNICATOR, ierr )

    IF ( myid /= rootproc ) THEN
      
      nPoints(1) = buffer(1)
      nPoints(2) = buffer(2)
      nPoints(3) = buffer(3)
      nVariables = buffer(4) 

  WRITE (*,*) "out process", myid, "buffer(4)", nVariables 

    END IF

    CALL AllocateEquationOfStateTable( EOSTable, nPoints , nVariables )

    IF ( myid == rootproc ) THEN

      CALL ReadThermoStateHDF( EOSTable % TS, file_id )

      CALL ReadDependentVariablesHDF( EOSTable % DV, file_id )

      CALL CloseFileHDF( file_id )

    END IF

    i_count = PRODUCT(nPoints) 
    
    nStates = 3

    DO i= 1, nStates
      CALL MPI_BCAST(EOSTable % TS % States(i) % Values(:), nPoints(i), &
                     MPI_DOUBLE_PRECISION, rootproc, COMMUNICATOR, ierr )
    END DO

    CALL MPI_BCAST(EOSTable % TS % Names(:), nStates,     &
                     MPI_DOUBLE_PRECISION, rootproc, COMMUNICATOR, ierr )
    CALL MPI_BCAST(EOSTable % TS % Units(:), nStates,     &
                     MPI_DOUBLE_PRECISION, rootproc, COMMUNICATOR, ierr )
    CALL MPI_BCAST(EOSTable % TS % minValues(:), nStates, &
                     MPI_DOUBLE_PRECISION, rootproc, COMMUNICATOR, ierr )
    CALL MPI_BCAST(EOSTable % TS % maxValues(:), nStates, &
                     MPI_DOUBLE_PRECISION, rootproc, COMMUNICATOR, ierr )

  WRITE (*,*) "process", myid, "test TS data", EOSTable % TS % States(1) % Values(10)
  WRITE (*,*) "process", myid, "test TS data", EOSTable % TS % Names(1) 
  WRITE (*,*) "process", myid, "test TS data", EOSTable % TS % Units(1) 

    DO i= 1, nVariables  
      CALL MPI_BCAST(EOSTable % DV % Variables(i) % Values(:,:,:), i_count,  &
                     MPI_DOUBLE_PRECISION, rootproc, COMMUNICATOR, ierr )
    END DO

    CALL MPI_BCAST(EOSTable % DV % Names(:), nVariables,                   &
                     MPI_DOUBLE_PRECISION, rootproc, COMMUNICATOR, ierr )
    CALL MPI_BCAST(EOSTable % DV % Units(:), nVariables,                   &
                     MPI_DOUBLE_PRECISION, rootproc, COMMUNICATOR, ierr )
    CALL MPI_BCAST(EOSTable % DV % Offsets(:), nVariables,                 &
                     MPI_DOUBLE_PRECISION, rootproc, COMMUNICATOR, ierr )
    CALL MPI_BCAST(EOSTable % DV % Repaired(:,:,:), i_count,               &
                   MPI_INTEGER, rootproc, COMMUNICATOR, ierr )

  END SUBROUTINE ReadEquationOfStateTableParallelHDF

  SUBROUTINE TransferDependentVariables( OldDV, NewDV, NewLocation, OldLocation )

    TYPE(DependentVariablesType), INTENT(inout)  :: OldDV
    TYPE(DependentVariablesType), INTENT(inout)  :: NewDV
    INTEGER, INTENT(in)                          :: OldLocation
    INTEGER, INTENT(in)                          :: NewLocation

     IF ( NewLocation == 0 ) THEN
       WRITE (*,*) "Dependent Variable", OldDV % Names( OldLocation ), "omitted"
       RETURN
     END IF

     NewDV % Variables( NewLocation ) % Values(:,:,:) &
               = OldDV % Variables( OldLocation ) % Values(:,:,:)
     NewDV % Names( NewLocation ) = OldDV % Names( OldLocation )
     NewDV % Units( NewLocation ) = OldDV % Units( OldLocation )
     NewDV % Offsets( NewLocation ) = OldDV % OffSets( OldLocation )

  END SUBROUTINE TransferDependentVariables

  SUBROUTINE MatchTableStructure( EOSTableIn, EOSTableOut, NewDVID, NewnVariables, filename )

    TYPE(EquationOfStateTableType), INTENT(inout) :: EOSTableIn
    TYPE(EquationOfStateTableType), INTENT(out)   :: EOSTableOut
    TYPE(DVIDType), INTENT(in)                    :: NewDVID
    INTEGER, INTENT(in)                           :: NewnVariables
    CHARACTER(len=*), INTENT(in)                  :: FileName
    INTEGER(HID_T)                 :: file_id
    INTEGER(HID_T)                 :: group_id


    CALL AllocateEquationOfStateTable( EOSTableOut, EOSTableIn % nPoints, &
                                     NewnVariables )

    EOSTableOut % DV % Indices = NewDVID
    EOSTableOut % DV % Repaired(:,:,:) = EOSTableIn % DV % Repaired(:,:,:)

    EOSTableOut % TS = EOSTableIn % TS

    ASSOCIATE( &

    NewiPressure => NewDVID % iPressure, &
    NewiEntropy => NewDVID % iEntropyPerBaryon, &
    NewiIntEnergy => NewDVID % iInternalEnergyDensity, &
    NewiEChemPot => NewDVID % iElectronChemicalPotential, &
    NewiPChemPot => NewDVID % iProtonChemicalPotential, &
    NewiNChemPot => NewDVID % iNeutronChemicalPotential, &
    NewiPMassFrac => NewDVID % iProtonMassFraction, &
    NewiNMassFrac => NewDVID % iNeutronMassFraction, &
    NewiAMassFrac => NewDVID % iAlphaMassFraction, &
    NewiHMassFrac => NewDVID % iHeavyMassFraction, &
    NewiHCharNum => NewDVID % iHeavyChargeNumber, &
    NewiHMassNum => NewDVID % iHeavyMassNumber, &
    NewiHeavyBE => NewDVID % iHeavyBindingEnergy, &
    NewiThermEnergy => NewDVID % iThermalEnergy, &
    NewiGamma1 => NewDVID % iGamma1, &

    OldiPressure => EOSTableIn % DV % Indices % iPressure , &
    OldiEntropy => EOSTableIn % DV % Indices % iEntropyPerBaryon, &
    OldiIntEnergy => EOSTableIn % DV % Indices % iInternalEnergyDensity, &
    OldiEChemPot => EOSTableIn % DV % Indices % iElectronChemicalPotential, &
    OldiPChemPot => EOSTableIn % DV % Indices % iProtonChemicalPotential, &
    OldiNChemPot => EOSTableIn % DV % Indices % iNeutronChemicalPotential, &
    OldiPMassFrac => EOSTableIn % DV % Indices % iProtonMassFraction, &
    OldiNMassFrac => EOSTableIn % DV % Indices % iNeutronMassFraction, &
    OldiAMassFrac => EOSTableIn % DV % Indices % iAlphaMassFraction, &
    OldiHMassFrac => EOSTableIn % DV % Indices % iHeavyMassFraction, &
    OldiHCharNum => EOSTableIn % DV % Indices % iHeavyChargeNumber, &
    OldiHMassNum => EOSTableIn % DV % Indices % iHeavyMassNumber, &
    OldiHeavyBE => EOSTableIn % DV % Indices % iHeavyBindingEnergy, &
    OldiThermEnergy => EOSTableIn % DV % Indices % iThermalEnergy, &
    OldiGamma1 => EOSTableIn % DV % Indices % iGamma1 )
    CALL TransferDependentVariables( EOSTableIn % DV, EOSTableOut % DV, &
                                     NewiPressure, OldiPressure )
    CALL TransferDependentVariables( EOSTableIn % DV, EOSTableOut % DV, &
                                     NewiEntropy, OldiEntropy )
    CALL TransferDependentVariables( EOSTableIn % DV, EOSTableOut % DV, &
                                     NewiIntEnergy, OldiIntEnergy )
    CALL TransferDependentVariables( EOSTableIn % DV, EOSTableOut % DV, &
                                     NewiEChemPot, OldiEChemPot )
    CALL TransferDependentVariables( EOSTableIn % DV, EOSTableOut % DV, &
                                     NewiPChemPot, OldiPChemPot )
    CALL TransferDependentVariables( EOSTableIn % DV, EOSTableOut % DV, &
                                     NewiNChemPot, OldiNChemPot )
    CALL TransferDependentVariables( EOSTableIn % DV, EOSTableOut % DV, &
                                     NewiPMassFrac, OldiPMassFrac )
    CALL TransferDependentVariables( EOSTableIn % DV, EOSTableOut % DV, &
                                     NewiNMassFrac, OldiNMassFrac )
    CALL TransferDependentVariables( EOSTableIn % DV, EOSTableOut % DV, &
                                     NewiAMassFrac, OldiAMassFrac )
    CALL TransferDependentVariables( EOSTableIn % DV, EOSTableOut % DV, &
                                     NewiHMassFrac, OldiHMassFrac )
    CALL TransferDependentVariables( EOSTableIn % DV, EOSTableOut % DV, &
                                     NewiHCharNum, OldiHCharNum )
    CALL TransferDependentVariables( EOSTableIn % DV, EOSTableOut % DV, &
                                     NewiHMassNum, OldiHMassNum )
    CALL TransferDependentVariables( EOSTableIn % DV, EOSTableOut % DV, &
                                     NewiHeavyBE, OldiHeavyBE )
    CALL TransferDependentVariables( EOSTableIn % DV, EOSTableOut % DV, &
                                     NewiThermEnergy, OldiThermEnergy )
    CALL TransferDependentVariables( EOSTableIn % DV, EOSTableOut % DV, &
                                     NewiGamma1, OldiGamma1 )

    CALL OpenFileHDF( filename, .true., file_id )

    CALL OpenGroupHDF( "ThermoState", .true., file_id, group_id )
    CALL WriteThermoStateHDF( EOSTableOut % TS, group_id )
    CALL CloseGroupHDF( group_id )

    CALL OpenGroupHDF( "DependentVariables", .true., file_id, group_id )
    CALL WriteDependentVariablesHDF( EOSTableOut % DV, group_id )
    CALL CloseGroupHDF( group_id )

    CALL CloseFileHDF( file_id )

    END ASSOCIATE

  END SUBROUTINE MatchTableStructure

END MODULE wlEOSIOModuleHDF

