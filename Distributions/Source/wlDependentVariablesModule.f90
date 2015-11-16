MODULE wlDependentVariablesModule

  USE wlKindModule, ONLY: dp

  implicit none
  PRIVATE


  TYPE :: ValueType
    REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: Values
  END TYPE

  TYPE, PUBLIC :: DVIDType
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
  END TYPE


  TYPE, PUBLIC :: DependentVariablesType
    INTEGER :: nVariables
    INTEGER, DIMENSION(3) :: nPoints
    CHARACTER(LEN=32), DIMENSION(:), ALLOCATABLE :: Names
    CHARACTER(LEN=32), DIMENSION(:), ALLOCATABLE :: Units
    REAL(dp), DIMENSION(:), ALLOCATABLE :: Offsets
    INTEGER, DIMENSION(:,:,:), ALLOCATABLE   :: Repaired
    TYPE(ValueType), DIMENSION(:), ALLOCATABLE :: Variables 
    TYPE(DVIDType) :: Indices
  END TYPE

  PUBLIC AllocateDependentVariables
  PUBLIC DeAllocateDependentVariables

CONTAINS

  SUBROUTINE AllocateDependentVariables( DV, nPoints, nVariables )
   
    TYPE(DependentVariablesType) :: DV 
    INTEGER, INTENT(in) :: nVariables
    INTEGER, DIMENSION(3), INTENT(in) :: nPoints

    INTEGER :: i

    ALLOCATE( DV % Names( nVariables ) )
    ALLOCATE( DV % Units( nVariables ) )
    ALLOCATE( DV % Offsets( nVariables ) ) 
    ALLOCATE( DV % Variables( nVariables ) ) 
    
    DV % nPoints = nPoints
    DV % nVariables = nVariables

    DO i = 1, nVariables
      ALLOCATE( DV % Variables(i) &
                   % Values( 1:DV % nPoints(1), 1:DV % nPoints(2), 1:DV % nPoints(3) ) ) 
    END DO
   
    ALLOCATE( DV % Repaired( 1:DV % nPoints(1), 1:DV % nPoints(2), 1:DV % nPoints(3) ) ) 

  END SUBROUTINE AllocateDependentVariables

  SUBROUTINE DeAllocateDependentVariables( DV )
  
    TYPE(DependentVariablesType) :: DV

    INTEGER :: i

    DO i = 1, SIZE( DV % Variables )
      DEALLOCATE( DV % Variables(i) % Values )
    END DO
   
    DEALLOCATE( DV % Variables )
    DEALLOCATE( DV % Names )
    DEALLOCATE( DV % Units )
    DEALLOCATE( DV % Repaired )

  END SUBROUTINE DeAllocateDependentVariables

END MODULE wlDependentVariablesModule

