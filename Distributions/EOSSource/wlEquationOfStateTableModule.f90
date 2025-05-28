MODULE wlEquationOfStateTableModule

  USE wlKindModule, ONLY: dp
  USE wlThermoStateModule
  USE wlDependentVariablesModule

  USE HDF5

  IMPLICIT NONE
  PRIVATE

  TYPE, PUBLIC :: MetadataType
    CHARACTER(LEN=120), DIMENSION(1) :: IDTag           ! Format: "wl-EOS-LS220-20-40-100, Energy zero, Date
    CHARACTER(LEN=120), DIMENSION(1) :: TableResolution ! "20 pts/dec rho, 40 pts/dec 
    CHARACTER(LEN=120), DIMENSION(1) :: NucEOSLink      ! Nuclear EOS ads link 
    CHARACTER(LEN=120), DIMENSION(1) :: LeptonEOSLink   ! Electron/photon EOS ads link 
    CHARACTER(LEN=120), DIMENSION(1) :: SourceLink      ! COMPOSE/SC download URL
    CHARACTER(LEN=120), DIMENSION(1) :: WLRevision      !
    CHARACTER(LEN=120), DIMENSION(1) :: TableLink       ! WeakLibTrac Link
  END TYPE

  TYPE, PUBLIC :: EquationOfStateTableType
    INTEGER                      :: nVariables
    INTEGER, DIMENSION(3)        :: nPoints
    TYPE(ThermoStateType)        :: TS
    TYPE(DependentVariablesType) :: DV
    TYPE(MetadataType)           :: MD
  END TYPE

  TYPE, PUBLIC :: EquationOfStateCompOSETableType
    INTEGER                      :: nVariables
    INTEGER, DIMENSION(3)        :: nPoints
    TYPE(ThermoStateType)        :: TS
    TYPE(DependentVariablesCompOSEType) :: DV
    TYPE(MetadataType)           :: MD
  END TYPE

  TYPE, PUBLIC :: EquationOfState4DTableType
    INTEGER                        :: nVariables
    INTEGER, DIMENSION(4)          :: nPoints
    TYPE(ThermoState4DType)        :: TS
    TYPE(DependentVariables4DType) :: DV
    TYPE(MetadataType)             :: MD
  END TYPE

  PUBLIC :: AllocateEquationOfStateTable
  PUBLIC :: DeAllocateEquationOfStateTable
  PUBLIC :: TableLimitFail
  PUBLIC :: SwapDependentVariables
  PUBLIC :: IndexMatch

  PUBLIC :: AllocateEquationOfStateCompOSETable
  PUBLIC :: DeAllocateEquationOfStateCompOSETable
  PUBLIC :: CompOSETableLimitFail
  PUBLIC :: SwapDependentVariablesCompOSE
  PUBLIC :: IndexMatchCompOSE

  PUBLIC :: AllocateEquationOfState4DTable
  PUBLIC :: DeAllocateEquationOfState4DTable

CONTAINS


  SUBROUTINE AllocateEquationOfStateTable( EOSTable, nPoints, nVariables )

    TYPE(EquationOfStateTableType), INTENT(inout) :: EOSTable
    INTEGER, INTENT(in)                           :: nVariables
    INTEGER, DIMENSION(3), INTENT(in)             :: nPoints

    EOSTable % nPoints(1:3) = nPoints(1:3)
    EOSTable % nVariables = nVariables

    CALL AllocateThermoState( EOSTable % TS, EOSTable % nPoints )
    CALL AllocateDependentVariables &
           ( EOSTable % DV, EOSTable % nPoints, EOSTable % nVariables )

  END SUBROUTINE AllocateEquationOfStateTable


  SUBROUTINE DeAllocateEquationOfStateTable( EOSTable )

    TYPE(EquationOfStateTableType) :: EOSTable

    CALL DeAllocateThermoState( EOSTable % TS )
    CALL DeAllocateDependentVariables( EOSTable % DV )

  END SUBROUTINE DeAllocateEquationOfStateTable


  LOGICAL FUNCTION TableLimitFail( rho, t, ye, EOSTable )

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

  
  ! ---------------------------------------------------------------------- !
  ! ------------------------------ CompOSE SECTION ----------------------- !
  ! ---------------------------------------------------------------------- !

  SUBROUTINE AllocateEquationOfStateCompOSETable( EOSTable, nPoints, nVariables )

    TYPE(EquationOfStateCompOSETableType), INTENT(inout) :: EOSTable
    INTEGER, INTENT(in)                           :: nVariables
    INTEGER, DIMENSION(3), INTENT(in)             :: nPoints

    EOSTable % nPoints(1:3) = nPoints(1:3)
    EOSTable % nVariables = nVariables

    CALL AllocateThermoState( EOSTable % TS, EOSTable % nPoints )
    CALL AllocateDependentVariablesCompOSE &
           ( EOSTable % DV, EOSTable % nPoints, EOSTable % nVariables )

  END SUBROUTINE AllocateEquationOfStateCompOSETable


  SUBROUTINE DeAllocateEquationOfStateCompOSETable( EOSTable )

    TYPE(EquationOfStateCompOSETableType) :: EOSTable

    CALL DeAllocateThermoState( EOSTable % TS )
    CALL DeAllocateDependentVariablesCompOSE( EOSTable % DV )

  END SUBROUTINE DeAllocateEquationOfStateCompOSETable

  LOGICAL FUNCTION CompOSETableLimitFail( rho, t, ye, EOSTable )

    REAL(dp), INTENT(in)                       :: rho, t, ye
    TYPE(EquationOfStateCompOSETableType), INTENT(in) :: EOSTable

    CompOSETableLimitFail = .false.
    IF ( rho < EOSTable % TS % minValues(1) ) CompOSETableLimitFail = .true.
    IF ( rho > EOSTable % TS % maxValues(1) ) CompOSETableLimitFail = .true.
    IF (   t < EOSTable % TS % minValues(2) ) CompOSETableLimitFail = .true.
    IF (   t > EOSTable % TS % maxValues(2) ) CompOSETableLimitFail = .true.
    IF (  ye < EOSTable % TS % minValues(3) ) CompOSETableLimitFail = .true.
    IF (  ye > EOSTable % TS % maxValues(3) ) CompOSETableLimitFail = .true.

  END FUNCTION CompOSETableLimitFail

  SUBROUTINE SwapDependentVariablesCompOSE( EOSTable, TargetBuffer, IndexBuffer )

    TYPE(EquationOfStateCompOSETableType), INTENT(inout)  :: EOSTable
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

      CALL IndexMatchCompOSE( TargetBuffer, IndexBuffer, &
                    EOSTable % DV % Indices )

  END SUBROUTINE SwapDependentVariablesCompOSE

  SUBROUTINE IndexMatchCompOSE( TargetBuffer, IndexBuffer, Indices )

    INTEGER, INTENT(in)        :: TargetBuffer
    INTEGER, INTENT(in)        :: IndexBuffer
    TYPE(DVIDCompOSEType), INTENT(inout) :: Indices  

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
      ELSE IF ( TargetBuffer == 10) THEN 
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
      ELSE IF (TargetBuffer == 16) THEN 
      Indices % iProtonEffMass = IndexBuffer 
      ELSE IF (TargetBuffer == 17) THEN 
      Indices % iNeutronEffMass = IndexBuffer 
      ELSE IF (TargetBuffer == 18) THEN 
      Indices % iProtonSelfEnergy = IndexBuffer 
      ELSE IF (TargetBuffer == 19) THEN 
      Indices % iNeutronSelfEnergy = IndexBuffer 
      END IF

  END SUBROUTINE IndexMatchCompOSE

  ! ---------------------------------------------------------------------- !
  ! --------------------------------- 4D SECTION ------------------------- !
  ! ---------------------------------------------------------------------- !

  SUBROUTINE AllocateEquationOfState4DTable( EOSTable, nPoints, nVariables )

    TYPE(EquationOfState4DTableType), INTENT(inout) :: EOSTable
    INTEGER, INTENT(in)                           :: nVariables
    INTEGER, DIMENSION(4), INTENT(in)             :: nPoints

    EOSTable % nPoints(1:4) = nPoints(1:4)
    EOSTable % nVariables = nVariables

    CALL AllocateThermoState4D( EOSTable % TS, EOSTable % nPoints )
    CALL AllocateDependentVariables4D &
           ( EOSTable % DV, EOSTable % nPoints, EOSTable % nVariables )

  END SUBROUTINE AllocateEquationOfState4DTable


  SUBROUTINE DeAllocateEquationOfState4DTable( EOSTable )

    TYPE(EquationOfState4DTableType) :: EOSTable

    CALL DeAllocateThermoState4D( EOSTable % TS )
    CALL DeAllocateDependentVariables4D( EOSTable % DV )

  END SUBROUTINE DeAllocateEquationOfState4DTable

END MODULE wlEquationOfStateTableModule
