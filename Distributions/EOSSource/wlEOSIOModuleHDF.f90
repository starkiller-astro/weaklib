MODULE wlEOSIOModuleHDF

  USE wlKindModule, ONLY: dp
  USE wlThermoStateModule
  USE wlDependentVariablesModule
  USE wlEquationOfStateTableModule
  USE wlIOModuleHDF

  USE HDF5 

  IMPLICIT NONE
  PRIVATE

  INTEGER :: hdferr

  PUBLIC WriteEquationOfStateTableHDF
  PUBLIC ReadEquationOfStateTableHDF
  PUBLIC DescribeEquationOfStateTable
  PUBLIC BroadcastEquationOfStateTableParallel
  PUBLIC MatchTableStructure
  PUBLIC TransferDependentVariables
  PUBLIC EOSVertexQuery
  PUBLIC ReadEOSMetadataHDF
  PUBLIC WriteEOSMetadataHDF

CONTAINS


  SUBROUTINE WriteEquationOfStateTableHDF( EOSTable, EOSTableName_Option )

    TYPE(EquationOfStateTableType), INTENT(inout)        :: EOSTable
    CHARACTER(len=*),               INTENT(in), OPTIONAL :: EOSTableName_Option
    
    INTEGER(HID_T) :: file_id
    INTEGER(HID_T) :: group_id
    CHARACTER(256) :: EOSTableName
 
    IF( PRESENT( EOSTableName_Option ) )THEN
       EOSTableName = TRIM( EOSTableName_Option )
    ELSE
       EOSTableName = 'EquationOfStateTable.h5'
    END IF
 
    CALL OpenFileHDF( EOSTableName, .true., file_id )

    CALL OpenGroupHDF( "ThermoState", .true., file_id, group_id )
    CALL WriteThermoStateHDF( EOSTable % TS, group_id )
    CALL CloseGroupHDF( group_id )

    CALL OpenGroupHDF( "DependentVariables", .true., file_id, group_id )
    CALL WriteDependentVariablesHDF( EOSTable % DV, group_id )
    CALL CloseGroupHDF( group_id )

    CALL OpenGroupHDF( "Metadata", .true., file_id, group_id )
    CALL WriteEOSMetadataHDF( EOSTable % MD, group_id )
    !CALL CloseGroupHDF( group_id )

    CALL CloseFileHDF( file_id )

  END SUBROUTINE WriteEquationOfStateTableHDF


  SUBROUTINE DescribeEquationOfStateTable( EOSTable )

    TYPE(EquationOfStateTableType), INTENT(inout) :: EOSTable

    INTEGER :: i

    ASSOCIATE( TS => EOSTable % TS )

    WRITE(*,*)
    WRITE(*,'(A2,A)') ' ', 'DescribeEquationOfStateTable'
    DO i = 1, 3

      WRITE(*,*)
      WRITE(*,'(A5,A21,I1.1,A3,A32)') &
        '', 'Independent Variable ', i, ' = ', TRIM( TS % Names(i) )
      WRITE(*,'(A7,A7,A20)') &
        '', 'Units: ', TRIM( TS % Units(i) )
      WRITE(*,'(A7,A12,ES12.4E2)') &
        '', 'Min Value = ', TS % minValues(i)
      WRITE(*,'(A7,A12,ES12.4E2)') &
        '', 'Max Value = ', TS % maxValues(i)
      WRITE(*,'(A7,A12,I4.4)') &
        '', 'nPoints   = ', TS % nPoints(i)

      IF ( TS % LogInterp(i) == 1 ) THEN
        WRITE (*,'(A7,A27)') &
          '', 'Grid Logarithmically Spaced'
      ELSE
        WRITE (*,'(A7,A20)') &
          '', 'Grid Linearly Spaced'
      END IF

    END DO

    END ASSOCIATE ! TS

    ASSOCIATE( DV => EOSTable % DV )

    WRITE(*,*)
    DO i = 1, EOSTable % nVariables

      WRITE (*,'(A5,A19,I3.3,A3,A32,A9,A20)') &
        '', 'Dependent Variable ', i, ' = ', TRIM( DV % Names(i) ), &
        '  Units: ', TRIM( DV % Units(i) )
    END DO
    WRITE(*,*)

    END ASSOCIATE ! DV

  END SUBROUTINE DescribeEquationOfStateTable

  SUBROUTINE ReadEOSMetadataHDF( MD, file_id )

    TYPE(MetadataType), INTENT(inout)           :: MD
    INTEGER(HID_T), INTENT(in)                  :: file_id

    INTEGER(HID_T)                              :: group_id
    INTEGER(HSIZE_T), DIMENSION(1)              :: datasize1d

    CALL OpenGroupHDF( "Metadata", .false., file_id, group_id )

    datasize1d = 1
    CALL ReadHDF( "IDTag", MD % IDTag(:), &
                              group_id, datasize1d )

    CALL ReadHDF( "TableResolution", MD % TableResolution(:), &
                              group_id, datasize1d )

    CALL ReadHDF( "NucEOSLink", MD % NucEOSLink(:), &
                              group_id, datasize1d )

    CALL ReadHDF( "LeptonEOSLink", MD % LeptonEOSLink(:), &
                              group_id, datasize1d )

    CALL ReadHDF( "SourceLink", MD % SourceLink(:), &
                              group_id, datasize1d )

    CALL ReadHDF( "WLRevision", MD % WLRevision(:), &
                              group_id, datasize1d )

    CALL ReadHDF( "TableLink", MD % TableLink(:), &
                              group_id, datasize1d )

    CALL CloseGroupHDF( group_id )

  END SUBROUTINE ReadEOSMetadataHDF

  SUBROUTINE WriteEOSMetadataHDF( MD, group_id )

    TYPE(MetadataType), INTENT(inout)           :: MD
    INTEGER(HID_T), INTENT(in)                  :: group_id
    INTEGER(HSIZE_T), DIMENSION(1)              :: datasize1d

    datasize1d = 1
    CALL WriteHDF( "IDTag", MD % IDTag(:), &
                              group_id, datasize1d )

    CALL WriteHDF( "TableResolution", MD % TableResolution(:), &
                              group_id, datasize1d )

    CALL WriteHDF( "NucEOSLink", MD % NucEOSLink(:), &
                              group_id, datasize1d )

    CALL WriteHDF( "LeptonEOSLink", MD % LeptonEOSLink(:), &
                              group_id, datasize1d )

    CALL WriteHDF( "SourceLink", MD % SourceLink(:), &
                              group_id, datasize1d )

    CALL WriteHDF( "WLRevision", MD % WLRevision(:), &
                              group_id, datasize1d )

    CALL WriteHDF( "TableLink", MD % TableLink(:), &
                              group_id, datasize1d )

    CALL CloseGroupHDF( group_id )

  END SUBROUTINE WriteEOSMetadataHDF
  

  SUBROUTINE ReadEquationOfStateTableHDF( EOSTable, FileName )

    TYPE(EquationOfStateTableType), INTENT(inout) :: EOSTable
    CHARACTER(len=*), INTENT(in)                  :: FileName

    INTEGER, DIMENSION(3)                         :: nPoints
    INTEGER                                       :: nVariables
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

    CALL ReadEOSMetadataHDF( EOSTable % MD, file_id )

!    CALL DescribeEquationOfStateTable( EOSTable )

    CALL CloseFileHDF( file_id )

  END SUBROUTINE ReadEquationOfStateTableHDF

  SUBROUTINE BroadcastEquationOfStateTableParallel &
               ( EOSTable, rootproc, myid, ierr,  COMMUNICATOR )

    USE MPI
    
    implicit none

    INTEGER, INTENT(in)                           :: rootproc
    INTEGER, INTENT(in)                           :: COMMUNICATOR
    INTEGER, INTENT(out)                          :: ierr !initialization variable for MPI
    INTEGER, INTENT(in)                           :: myid ! rank of each processor (MPI)   
    TYPE(EquationOfStateTableType), INTENT(inout) :: EOSTable
    INTEGER, DIMENSION(3)                         :: nPoints
    INTEGER, DIMENSION(22)                        :: buffer
    INTEGER                                       :: nStates, nVariables, i
    INTEGER                                       :: i_count
    INTEGER                                       :: charlen
    CHARACTER(LEN=120), DIMENSION(7)              :: sendstring

    IF ( myid == rootproc ) THEN
    
      buffer( 1) = EOSTable % nPoints(1)
      buffer( 2) = EOSTable % nPoints(2)
      buffer( 3) = EOSTable % nPoints(3)
      buffer( 4) = EOSTable % nVariables
      buffer( 5) = EOSTable % DV % Indices % iPressure
      buffer( 6) = EOSTable % DV % Indices % iEntropyPerBaryon
      buffer( 7) = EOSTable % DV % Indices % iInternalEnergyDensity
      buffer( 8) = EOSTable % DV % Indices % iElectronChemicalPotential
      buffer( 9) = EOSTable % DV % Indices % iProtonChemicalPotential
      buffer(10) = EOSTable % DV % Indices % iNeutronChemicalPotential
      buffer(11) = EOSTable % DV % Indices % iProtonMassFraction
      buffer(12) = EOSTable % DV % Indices % iNeutronMassFraction
      buffer(13) = EOSTable % DV % Indices % iAlphaMassFraction
      buffer(14) = EOSTable % DV % Indices % iHeavyMassFraction
      buffer(15) = EOSTable % DV % Indices % iHeavyChargeNumber
      buffer(16) = EOSTable % DV % Indices % iHeavyMassNumber
      buffer(17) = EOSTable % DV % Indices % iHeavyBindingEnergy
      buffer(18) = EOSTable % DV % Indices % iThermalEnergy
      buffer(19) = EOSTable % DV % Indices % iGamma1
      buffer(20) = EOSTable % TS % Indices % iRho
      buffer(21) = EOSTable % TS % Indices % iT
      buffer(22) = EOSTable % TS % Indices % iYe

    END IF

    i_count = SIZE(buffer) 

    CALL MPI_BCAST( buffer, i_count, MPI_INTEGER, rootproc, COMMUNICATOR, ierr )
      
    IF ( myid /= rootproc ) THEN
      
      nPoints(1) = buffer(1)
      nPoints(2) = buffer(2)
      nPoints(3) = buffer(3)
      nVariables = buffer(4) 

      CALL AllocateEquationOfStateTable( EOSTable, nPoints , nVariables )

      EOSTable % DV % Indices % iPressure                  = buffer( 5)
      EOSTable % DV % Indices % iEntropyPerBaryon          = buffer( 6)
      EOSTable % DV % Indices % iInternalEnergyDensity     = buffer( 7)
      EOSTable % DV % Indices % iElectronChemicalPotential = buffer( 8)
      EOSTable % DV % Indices % iProtonChemicalPotential   = buffer( 9)
      EOSTable % DV % Indices % iNeutronChemicalPotential  = buffer(10)
      EOSTable % DV % Indices % iProtonMassFraction        = buffer(11)
      EOSTable % DV % Indices % iNeutronMassFraction       = buffer(12)
      EOSTable % DV % Indices % iAlphaMassFraction         = buffer(13)
      EOSTable % DV % Indices % iHeavyMassFraction         = buffer(14)
      EOSTable % DV % Indices % iHeavyChargeNumber         = buffer(15)
      EOSTable % DV % Indices % iHeavyMassNumber           = buffer(16)
      EOSTable % DV % Indices % iHeavyBindingEnergy        = buffer(17)
      EOSTable % DV % Indices % iThermalEnergy             = buffer(18)
      EOSTable % DV % Indices % iGamma1                    = buffer(19)
      EOSTable % TS % Indices % iRho                       = buffer(20)
      EOSTable % TS % Indices % iT                         = buffer(21)
      EOSTable % TS % Indices % iYe                        = buffer(22)

    END IF

    i_count = PRODUCT(EOSTable % nPoints) 
    
    nStates = 3

    DO i= 1, nStates
      CALL MPI_BCAST(EOSTable % TS % States(i) % Values(:), EOSTable % nPoints(i), &
                     MPI_DOUBLE_PRECISION, rootproc, COMMUNICATOR, ierr )
    END DO

    CALL MPI_BCAST(EOSTable % TS % Names(:), nStates*32,                       &
                     MPI_CHARACTER, rootproc, COMMUNICATOR, ierr )
    CALL MPI_BCAST(EOSTable % TS % Units(:), nStates*32,                       &
                     MPI_CHARACTER, rootproc, COMMUNICATOR, ierr )
    CALL MPI_BCAST(EOSTable % TS % minValues(:), nStates,                      &
                     MPI_DOUBLE_PRECISION, rootproc, COMMUNICATOR, ierr )
    CALL MPI_BCAST(EOSTable % TS % maxValues(:), nStates,                      &
                     MPI_DOUBLE_PRECISION, rootproc, COMMUNICATOR, ierr )
    CALL MPI_BCAST(EOSTable % TS % LogInterp(:), nStates,                      &
                     MPI_INTEGER, rootproc, COMMUNICATOR, ierr )


    DO i= 1, EOSTable % nVariables  
      CALL MPI_BCAST(EOSTable % DV % Variables(i) % Values(:,:,:), i_count,    &
                     MPI_DOUBLE_PRECISION, rootproc, COMMUNICATOR, ierr )
    END DO

    CALL MPI_BCAST(EOSTable % DV % Names(:), EOSTable % nVariables*32,         &
                     MPI_CHARACTER, rootproc, COMMUNICATOR, ierr )
    CALL MPI_BCAST(EOSTable % DV % Units(:), EOSTable % nVariables*32,         &
                     MPI_CHARACTER, rootproc, COMMUNICATOR, ierr )
    CALL MPI_BCAST(EOSTable % DV % Offsets(:), EOSTable % nVariables,          &
                     MPI_DOUBLE_PRECISION, rootproc, COMMUNICATOR, ierr )
    CALL MPI_BCAST(EOSTable % DV % minValues(:), EOSTable % nVariables,        &
                     MPI_DOUBLE_PRECISION, rootproc, COMMUNICATOR, ierr )
    CALL MPI_BCAST(EOSTable % DV % maxValues(:), EOSTable % nVariables,        &
                     MPI_DOUBLE_PRECISION, rootproc, COMMUNICATOR, ierr )
    CALL MPI_BCAST(EOSTable % DV % Repaired(:,:,:), i_count,                   &
                   MPI_INTEGER, rootproc, COMMUNICATOR, ierr )

    sendstring(1) = EOSTable % MD % IDTag(1)
    sendstring(2) = EOSTable % MD % TableResolution(1)
    sendstring(3) = EOSTable % MD % NucEOSLink(1)
    sendstring(4) = EOSTable % MD % LeptonEOSLink(1)
    sendstring(5) = EOSTable % MD % SourceLink(1)
    sendstring(6) = EOSTable % MD % WLRevision(1)
    sendstring(7) = EOSTable % MD % TableLink(1)
    charlen       = 7 * 120

    CALL MPI_BCAST( sendstring, charlen, &
                     MPI_CHARACTER, rootproc, COMMUNICATOR, ierr )

    EOSTable % MD % IDTag(1)           = sendstring(1)
    EOSTable % MD % TableResolution(1) = sendstring(2)
    EOSTable % MD % NucEOSLink(1)      = sendstring(3) 
    EOSTable % MD % LeptonEOSLink(1)   = sendstring(4) 
    EOSTable % MD % SourceLink(1)      = sendstring(5) 
    EOSTable % MD % WLRevision(1)      = sendstring(6) 
    EOSTable % MD % TableLink(1)       = sendstring(7) 

  END SUBROUTINE BroadcastEquationOfStateTableParallel


  SUBROUTINE TransferDependentVariables( OldDV, NewDV, NewLocation, OldLocation )

    TYPE(DependentVariablesType), INTENT(in)     :: OldDV
    TYPE(DependentVariablesType), INTENT(inout)  :: NewDV
    INTEGER, INTENT(in)                          :: OldLocation
    INTEGER, INTENT(in)                          :: NewLocation

     IF ( NewLocation == 0 ) THEN
       WRITE (*,*) "Dependent Variable  ", OldDV % Names( OldLocation ), &
                   "omitted"
       RETURN
     END IF

     NewDV % Variables( NewLocation ) % Values(:,:,:) &
               = OldDV % Variables( OldLocation ) % Values(:,:,:)
     NewDV % Names( NewLocation ) = OldDV % Names( OldLocation )
     NewDV % Units( NewLocation ) = OldDV % Units( OldLocation )
     NewDV % Offsets( NewLocation ) = OldDV % OffSets( OldLocation )

  END SUBROUTINE TransferDependentVariables


  SUBROUTINE MatchTableStructure( EOSTableIn, EOSTableOut, NewDVID, NewnVariables )

    TYPE(EquationOfStateTableType), INTENT(in)    :: EOSTableIn
    TYPE(EquationOfStateTableType), INTENT(out)   :: EOSTableOut
    TYPE(DVIDType), INTENT(in)                    :: NewDVID
    INTEGER, INTENT(in)                           :: NewnVariables

    CALL AllocateEquationOfStateTable( EOSTableOut, EOSTableIn % nPoints, &
                                     NewnVariables )

    EOSTableOut % DV % Indices = NewDVID
    EOSTableOut % DV % Repaired(:,:,:) = EOSTableIn % DV % Repaired(:,:,:)

    EOSTableOut % TS = EOSTableIn % TS
    EOSTableOut % MD = EOSTableIn % MD

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
    EOSTableOut % DV % Indices % iPressure = NewiPressure
    CALL TransferDependentVariables( EOSTableIn % DV, EOSTableOut % DV, &
                                     NewiEntropy, OldiEntropy )
    EOSTableOut % DV % Indices % iEntropyPerBaryon = NewiEntropy
    CALL TransferDependentVariables( EOSTableIn % DV, EOSTableOut % DV, &
                                     NewiIntEnergy, OldiIntEnergy )
    EOSTableOut % DV % Indices % iInternalEnergyDensity = NewiIntEnergy
    CALL TransferDependentVariables( EOSTableIn % DV, EOSTableOut % DV, &
                                     NewiEChemPot, OldiEChemPot )
    EOSTableOut % DV % Indices % iElectronChemicalPotential = NewiEChemPot
    CALL TransferDependentVariables( EOSTableIn % DV, EOSTableOut % DV, &
                                     NewiPChemPot, OldiPChemPot )
    EOSTableOut % DV % Indices % iProtonChemicalPotential = NewiPChemPot
    CALL TransferDependentVariables( EOSTableIn % DV, EOSTableOut % DV, &
                                     NewiNChemPot, OldiNChemPot )
    EOSTableOut % DV % Indices % iNeutronChemicalPotential = NewiNChemPot
    CALL TransferDependentVariables( EOSTableIn % DV, EOSTableOut % DV, &
                                     NewiPMassFrac, OldiPMassFrac )
    EOSTableOut % DV % Indices % iProtonMassFraction = NewiPMassFrac
    CALL TransferDependentVariables( EOSTableIn % DV, EOSTableOut % DV, &
                                     NewiNMassFrac, OldiNMassFrac )
    EOSTableOut % DV % Indices % iNeutronMassFraction = NewiNMassFrac
    CALL TransferDependentVariables( EOSTableIn % DV, EOSTableOut % DV, &
                                     NewiAMassFrac, OldiAMassFrac )
    EOSTableOut % DV % Indices % iAlphaMassFraction = NewiAMassFrac
    CALL TransferDependentVariables( EOSTableIn % DV, EOSTableOut % DV, &
                                     NewiHMassFrac, OldiHMassFrac )
    EOSTableOut % DV % Indices % iHeavyMassFraction = NewiHMassFrac
    CALL TransferDependentVariables( EOSTableIn % DV, EOSTableOut % DV, &
                                     NewiHCharNum, OldiHCharNum )
    EOSTableOut % DV % Indices % iHeavyChargeNumber = NewiHCharNum
    CALL TransferDependentVariables( EOSTableIn % DV, EOSTableOut % DV, &
                                     NewiHMassNum, OldiHMassNum )
    EOSTableOut % DV % Indices % iHeavyMassNumber = NewiHMassNum
    CALL TransferDependentVariables( EOSTableIn % DV, EOSTableOut % DV, &
                                     NewiHeavyBE, OldiHeavyBE )
    EOSTableOut % DV % Indices % iHeavyBindingEnergy = NewiHeavyBE
    CALL TransferDependentVariables( EOSTableIn % DV, EOSTableOut % DV, &
                                     NewiThermEnergy, OldiThermEnergy )
    EOSTableOut % DV % Indices % iThermalEnergy = NewiThermEnergy
    CALL TransferDependentVariables( EOSTableIn % DV, EOSTableOut % DV, &
                                     NewiGamma1, OldiGamma1 )
    EOSTableOut % DV % Indices % iGamma1 = NewiGamma1

    END ASSOCIATE

  END SUBROUTINE MatchTableStructure


  SUBROUTINE EOSVertexQuery( irho, iT, iYe, EOSTable, Values )

    INTEGER                                    :: i, j
    INTEGER, DIMENSION(:), INTENT(in)          :: irho
    INTEGER, DIMENSION(:), INTENT(in)          :: iT
    INTEGER, DIMENSION(:), INTENT(in)          :: iYe
    TYPE(EquationOfStateTableType), INTENT(in) :: EOSTable
    REAL(dp), DIMENSION(:,:), INTENT(out)      :: Values

    DO i = 1, SIZE(irho)
      DO j = 1, EOSTable % DV % nVariables
        Values(i,j) &
          = 10.d0**( EOSTable % DV % Variables(j) % Values&
          ( irho(i), iT(i), iYe(i) ) ) -               &
          EOSTable % DV % Offsets(j)
      END DO
    END DO

  END SUBROUTINE EOSVertexQuery


END MODULE wlEOSIOModuleHDF

