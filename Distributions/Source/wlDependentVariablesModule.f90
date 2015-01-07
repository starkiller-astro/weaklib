MODULE wlDependentVariablesModule

  USE wlKindModule, ONLY: dp

  implicit none
  PRIVATE

  TYPE :: ValueType
    REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: Values
  END TYPE

  TYPE, PUBLIC :: DependentVariablesType
    CHARACTER(LEN=32), DIMENSION(:), ALLOCATABLE :: Names
    INTEGER :: nVariables
    TYPE(ValueType), DIMENSION(:), ALLOCATABLE :: Variables 
  END TYPE

  PUBLIC AllocateDependentVariables
  PUBLIC DeAllocateDependentVariables

CONTAINS

  SUBROUTINE AllocateDependentVariables( DV, nValues, nVariables )
   
    TYPE(DependentVariablesType) :: DV 
    INTEGER, INTENT(in) :: nVariables
    INTEGER, DIMENSION(3), INTENT(in) :: nValues

    INTEGER :: i

    ALLOCATE( DV % Variables( nVariables ) ) 

    DO i = 1, nVariables
      ALLOCATE( DV % Variables(i) &
                   % Values( 1:nValues(1), 1:nValues(2), 1:nValues(3) ) ) 
    END DO

  END SUBROUTINE AllocateDependentVariables

  SUBROUTINE DeAllocateDependentVariables( DV )
  
    TYPE(DependentVariablesType) :: DV

    INTEGER :: i

    DO i = 1, SIZE( DV % Variables )
      DEALLOCATE( DV % Variables(i) % Values )
    END DO
   
    DEALLOCATE( DV % Variables )

  END SUBROUTINE DeAllocateDependentVariables

END MODULE wlDependentVariablesModule

