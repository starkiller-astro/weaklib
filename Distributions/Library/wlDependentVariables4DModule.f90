MODULE wlDependentVariables4DModule

  USE wlKindModule, ONLY: dp

  implicit none
  PRIVATE


  TYPE, PUBLIC :: ValueTypeDV4D
    REAL(dp), ALLOCATABLE, DIMENSION(:,:,:,:) :: Values
  END TYPE

  TYPE, PUBLIC :: DVIDType4D
    INTEGER :: iPressure
    INTEGER :: iEntropyPerBaryon
    INTEGER :: iInternalEnergyDensity
    INTEGER :: iElectronChemicalPotential
    INTEGER :: iProtonChemicalPotential
    INTEGER :: iNeutronChemicalPotential
    INTEGER :: iProtonMassFraction
    INTEGER :: iNeutronMassFraction
    INTEGER :: iAlphaMassFraction
    INTEGER :: iHeavyMassFraction
    INTEGER :: iHeavyChargeNumber
    INTEGER :: iHeavyMassNumber
    INTEGER :: iHeavyBindingEnergy
    INTEGER :: iThermalEnergy
    INTEGER :: iGamma1
    INTEGER :: iProtonSelfEnergy
    INTEGER :: iNeutronSelfEnergy
  END TYPE

  TYPE, PUBLIC :: DependentVariables4DType
    INTEGER :: nVariables
    INTEGER, DIMENSION(4) :: nPoints
    CHARACTER(LEN=32), DIMENSION(:), ALLOCATABLE :: Names
    CHARACTER(LEN=32), DIMENSION(:), ALLOCATABLE :: Units
    REAL(dp), DIMENSION(:), ALLOCATABLE :: Offsets
    INTEGER, DIMENSION(:,:,:,:), ALLOCATABLE  :: Repaired
    REAL(dp), DIMENSION(:), ALLOCATABLE :: minValues
    REAL(dp), DIMENSION(:), ALLOCATABLE :: maxValues
    TYPE(ValueTypeDV4D), DIMENSION(:), ALLOCATABLE :: Variables
    TYPE(DVIDType4D) :: Indices
  END TYPE

  PUBLIC AllocateDependentVariables4D
  PUBLIC DeAllocateDependentVariables4D

CONTAINS


  SUBROUTINE AllocateDependentVariables4D( DV, nPoints, nVariables )
   
    TYPE(DependentVariables4DType)      :: DV 
    INTEGER,               INTENT(in) :: nVariables
    INTEGER, DIMENSION(4), INTENT(in) :: nPoints

    INTEGER :: i

    DV % nPoints    = nPoints
    DV % nVariables = nVariables

    ALLOCATE( DV % Names( nVariables ) )
    ALLOCATE( DV % Units( nVariables ) )
    ALLOCATE( DV % Offsets( nVariables ) ) 
    ALLOCATE( DV % Variables( nVariables ) ) 
    ALLOCATE( DV % minValues( nVariables ) ) 
    ALLOCATE( DV % maxValues( nVariables ) ) 

    ALLOCATE( DV % Repaired(1:nPoints(1), 1:nPoints(2), 1:nPoints(3), 1:nPoints(4)) )

    DO i = 1, nVariables
      ALLOCATE &
        ( DV % Variables(i) &
             % Values(1:nPoints(1), 1:nPoints(2), 1:nPoints(3), 1:nPoints(4)) )
    END DO

  END SUBROUTINE AllocateDependentVariables4D


  SUBROUTINE DeAllocateDependentVariables4D( DV )
  
    TYPE(DependentVariables4DType) :: DV

    INTEGER :: i

    DO i = 1, DV % nVariables
      DEALLOCATE( DV % Variables(i) % Values )
    END DO
   
    DEALLOCATE( DV % Repaired )

    DEALLOCATE( DV % Variables )
    DEALLOCATE( DV % Offsets )
    DEALLOCATE( DV % Units )
    DEALLOCATE( DV % Names )
    DEALLOCATE( DV % minValues )
    DEALLOCATE( DV % maxValues )

  END SUBROUTINE DeAllocateDependentVariables4D


END MODULE wlDependentVariables4DModule

