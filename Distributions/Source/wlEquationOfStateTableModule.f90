MODULE wlEquationOfStateTableModule

  USE wlKindModule, ONLY: dp
  USE HDF5
  USE wlThermoStateModule
  USE wlDependentVariablesModule

  implicit none
  PRIVATE

  TYPE, PUBLIC :: EquationOfStateTableType
    !CHARACTER(LEN=32), DIMENSION(:), ALLOCATABLE :: Name
    INTEGER                                      :: nVariables
    INTEGER, DIMENSION(3)                        :: nPoints
    TYPE(ThermoStateType)           :: TS 
    TYPE(DependentVariablesType)    :: DV
  END TYPE


  PUBLIC AllocateEquationOfStateTable
  PUBLIC DeAllocateEquationOfStateTable
  PUBLIC TableLimitFail
  PUBLIC MatchTableStructure

CONTAINS 

  SUBROUTINE AllocateEquationOfStateTable( EOSTable, nPoints, nVariables )

    TYPE(EquationOfStateTableType), INTENT(inout)    :: EOSTable
    INTEGER, INTENT(in)               :: nVariables
    INTEGER, DIMENSION(3), INTENT(in) :: nPoints
   
    EOSTable % nPoints(1:3) = nPoints(1:3)
    EOSTable % nVariables = nVariables

    CALL AllocateThermoState( EOSTable % TS, EOSTable % nPoints )
    CALL AllocateDependentVariables( EOSTable % DV, EOSTable % nPoints, &
                                     EOSTable % nVariables ) 

  END SUBROUTINE AllocateEquationOfStateTable

  SUBROUTINE DeAllocateEquationOfStateTable( EOSTable )

    TYPE(EquationOfStateTableType) :: EOSTable

    CALL DeAllocateThermoState( EOSTable % TS )
    CALL DeAllocateDependentVariables( EOSTable % DV )

  END SUBROUTINE DeAllocateEquationOfStateTable

  LOGICAL FUNCTION TableLimitFail( rho, t, ye, EOSTable )

    !LOGICAL                                    :: TableLimitFail
    REAL(dp), INTENT(in)                       :: rho, t, ye
    TYPE(EquationOfStateTableType), INTENT(in) :: EOSTable

      TableLimitFail = .false.
      IF ( rho < EOSTable % TS % minValues(1) ) TableLimitFail = .true.
      IF ( rho > EOSTable % TS % maxValues(1) ) TableLimitFail = .true.
      IF (   t < EOSTable % TS % minValues(2) ) TableLimitFail = .true.
      IF (   t > EOSTable % TS % maxValues(2) ) TableLimitFail = .true.
      IF (  ye < EOSTable % TS % minValues(3) ) TableLimitFail = .true.
      IF (  ye > EOSTable % TS % maxValues(3) ) TableLimitFail = .true.

  END FUNCTION TableLimitFail

  SUBROUTINE MatchTableStructure( EOSTable, LocalDVID )
  
    TYPE(EquationOfStateTableType), INTENT(inout)  :: EOSTable
    TYPE(DVIDType), INTENT(in)                     :: LocalDVID
    TYPE(DependentVariablesType)                   :: LocalDV

    CALL AllocateDependentVariables( LocalDV, EOSTable % nPoints, &
                                     EOSTable % nVariables )
    LocalDV % Indices = LocalDVID
    LocalDV % Repaired(:,:,:) = EOSTable % DV % Repaired(:,:,:)

    ASSOCIATE( &
    
    NewiPressure => LocalDVID % iPressure, &
    NewiEntropy => LocalDVID % iEntropyPerBaryon, &
    NewiIntEnergy => LocalDVID % iInternalEnergyDensity, &
    NewiEChemPot => LocalDVID % iElectronChemicalPotential, &
    NewiPChemPot => LocalDVID % iProtonChemicalPotential, &
    NewiNChemPot => LocalDVID % iNeutronChemicalPotential, &
    NewiPMassFrac => LocalDVID % iProtonMassFraction, &
    NewiNMassFrac => LocalDVID % iNeutronMassFraction, &
    NewiAMassFrac => LocalDVID % iAlphaMassFraction, &
    NewiHMassFrac => LocalDVID % iHeavyMassFraction, &
    NewiHCharNum => LocalDVID % iHeavyChargeNumber, &
    NewiHMassNum => LocalDVID % iHeavyMassNumber, &
    NewiHeavyBE => LocalDVID % iHeavyBindingEnergy, &
    NewiThermEnergy => LocalDVID % iThermalEnergy, &
    NewiGamma1 => LocalDVID % iGamma1, &
    OldiPressure => EOSTable % DV % Indices % iPressure , &
    OldiEntropy => EOSTable % DV % Indices % iEntropyPerBaryon, &
    OldiIntEnergy => EOSTable % DV % Indices % iInternalEnergyDensity, &
    OldiEChemPot => EOSTable % DV % Indices % iElectronChemicalPotential, &
    OldiPChemPot => EOSTable % DV % Indices % iProtonChemicalPotential, &
    OldiNChemPot => EOSTable % DV % Indices % iNeutronChemicalPotential, &
    OldiPMassFrac => EOSTable % DV % Indices % iProtonMassFraction, &
    OldiNMassFrac => EOSTable % DV % Indices % iNeutronMassFraction, &
    OldiAMassFrac => EOSTable % DV % Indices % iAlphaMassFraction, &
    OldiHMassFrac => EOSTable % DV % Indices % iHeavyMassFraction, &
    OldiHCharNum => EOSTable % DV % Indices % iHeavyChargeNumber, &
    OldiHMassNum => EOSTable % DV % Indices % iHeavyMassNumber, &
    OldiHeavyBE => EOSTable % DV % Indices % iHeavyBindingEnergy, &
    OldiThermEnergy => EOSTable % DV % Indices % iThermalEnergy, &
    OldiGamma1 => EOSTable % DV % Indices % iGamma1 )


    CALL TransferDependentVariables( EOSTable % DV, LocalDV, NewiPressure, OldiPressure )
    CALL TransferDependentVariables( EOSTable % DV, LocalDV, NewiEntropy, OldiEntropy )
    CALL TransferDependentVariables( EOSTable % DV, LocalDV, NewiIntEnergy, OldiIntEnergy )
    CALL TransferDependentVariables( EOSTable % DV, LocalDV, NewiEChemPot, OldiEChemPot )
    CALL TransferDependentVariables( EOSTable % DV, LocalDV, NewiPChemPot, OldiPChemPot )
    CALL TransferDependentVariables( EOSTable % DV, LocalDV, NewiNChemPot, OldiNChemPot )
    CALL TransferDependentVariables( EOSTable % DV, LocalDV, NewiPMassFrac, OldiPMassFrac )
    CALL TransferDependentVariables( EOSTable % DV, LocalDV, NewiNMassFrac, OldiNMassFrac )
    CALL TransferDependentVariables( EOSTable % DV, LocalDV, NewiAMassFrac, OldiAMassFrac )
    CALL TransferDependentVariables( EOSTable % DV, LocalDV, NewiHMassFrac, OldiHMassFrac )
    CALL TransferDependentVariables( EOSTable % DV, LocalDV, NewiHCharNum, OldiHCharNum )
    CALL TransferDependentVariables( EOSTable % DV, LocalDV, NewiHMassNum, OldiHMassNum )
    CALL TransferDependentVariables( EOSTable % DV, LocalDV, NewiHeavyBE, OldiHeavyBE )
    CALL TransferDependentVariables( EOSTable % DV, LocalDV, NewiThermEnergy, OldiThermEnergy )
    CALL TransferDependentVariables( EOSTable % DV, LocalDV, NewiGamma1, OldiGamma1 )

    EOSTable % DV = LocalDV

    END ASSOCIATE


  END SUBROUTINE MatchTableStructure

  SUBROUTINE SwapDependentVariables( EOSTable, TargetBuffer, IndexBuffer )

    TYPE(EquationOfStateTableType), INTENT(inout)  :: EOSTable
    CHARACTER(LEN=32)                              :: NameBuffer
    CHARACTER(LEN=32)                              :: UnitBuffer
    INTEGER, INTENT(in)                            :: IndexBuffer
    INTEGER, INTENT(in)                            :: TargetBuffer
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE        :: ValuesBuffer

    ALLOCATE( ValuesBuffer( EOSTable % DV % nPoints(1), &
               EOSTable % DV % nPoints(2), EOSTable % DV % nPoints(3) ) )

      ValuesBuffer(:,:,:) = EOSTable % DV % Variables( TargetBuffer ) % Values(:,:,:)
      NameBuffer = EOSTable % DV % Names( TargetBuffer )
      UnitBuffer = EOSTable % DV % Units( TargetBuffer )

      EOSTable % DV % Variables( TargetBuffer ) % Values(:,:,:) &
                = EOSTable % DV % Variables( IndexBuffer ) % Values(:,:,:)
      EOSTable % DV % Names( TargetBuffer ) = EOSTable % DV % Names( IndexBuffer )
      EOSTable % DV % Units( TargetBuffer ) = EOSTable % DV % Units( IndexBuffer )

      EOSTable % DV % Variables( IndexBuffer ) % Values(:,:,:) = ValuesBuffer(:,:,:)
      EOSTable % DV % Names( IndexBuffer ) = NameBuffer
      EOSTable % DV % Units( IndexBuffer ) = UnitBuffer


      CALL IndexMatch( TargetBuffer, IndexBuffer, &
                       EOSTable % DV % Indices )

  END SUBROUTINE SwapDependentVariables

  SUBROUTINE TransferDependentVariables( DV, LocalDV, NewLocation, OldLocation )

    TYPE(DependentVariablesType), INTENT(inout)  :: DV 
    TYPE(DependentVariablesType), INTENT(inout)  :: LocalDV 
    INTEGER, INTENT(in)                          :: OldLocation
    INTEGER, INTENT(in)                          :: NewLocation


      LocalDV % Variables( NewLocation ) % Values(:,:,:) &
                = DV % Variables( OldLocation ) % Values(:,:,:)
      LocalDV % Names( NewLocation ) = DV % Names( OldLocation )
      LocalDV % Units( NewLocation ) = DV % Units( OldLocation )
      LocalDV % Offsets( NewLocation ) = DV % OffSets( OldLocation )

  END SUBROUTINE TransferDependentVariables

  SUBROUTINE IndexMatch( TargetBuffer, IndexBuffer, Indices )

    INTEGER, INTENT(in)        :: TargetBuffer
    INTEGER, INTENT(in)        :: IndexBuffer
    TYPE(DVIDType), INTENT(inout) :: Indices  

      IF ( TargetBuffer == 1) THEN 
      Indices % iPressure = IndexBuffer 
      ELSE IF ( TargetBuffer == 2) THEN 
      Indices % iEntropyPerBaryon = IndexBuffer 
      ELSE IF ( TargetBuffer == 3) THEN 
      Indices % iInternalEnergyDensity = IndexBuffer 
      ELSE IF ( TargetBuffer == 4) THEN 
      Indices % iElectronChemicalPotential = IndexBuffer 
      ELSE IF ( TargetBuffer == 5) THEN 
      Indices % iProtonChemicalPotential = IndexBuffer 
      ELSE IF ( TargetBuffer == 6) THEN 
      Indices % iNeutronChemicalPotential = IndexBuffer 
      ELSE IF ( TargetBuffer == 7) THEN 
      Indices % iProtonMassFraction = IndexBuffer 
      ELSE IF ( TargetBuffer == 8) THEN 
      Indices % iNeutronMassFraction = IndexBuffer 
      ELSE IF ( TargetBuffer == 9) THEN 
      Indices % iAlphaMassFraction = IndexBuffer 
      ELSE IF (TargetBuffer == 10) THEN 
      Indices % iHeavyMassFraction = IndexBuffer 
      ELSE IF (TargetBuffer == 11) THEN 
      Indices % iHeavyChargeNumber = IndexBuffer 
      ELSE IF (TargetBuffer == 12) THEN 
      Indices % iHeavyMassNumber = IndexBuffer 
      ELSE IF (TargetBuffer == 13) THEN 
      Indices % iHeavyBindingEnergy = IndexBuffer 
      ELSE IF (TargetBuffer == 14) THEN 
      Indices % iThermalEnergy = IndexBuffer 
      ELSE IF (TargetBuffer == 15) THEN 
      Indices % iGamma1 = IndexBuffer 
      END IF

  END SUBROUTINE IndexMatch


END MODULE wlEquationOfStateTableModule
