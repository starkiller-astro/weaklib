MODULE wlThermoStateModule
  
  USE wlKindModule, ONLY: dp

  implicit none
  PRIVATE

  TYPE :: ValueType
    REAL(dp), ALLOCATABLE, DIMENSION(:) :: Values
  END TYPE
  
  TYPE, PUBLIC :: ThermoStateType
    CHARACTER(LEN=32), DIMENSION(3) :: Names
    INTEGER, DIMENSION(3) :: nValues
    REAL(dp), DIMENSION(3) :: minValues
    REAL(dp), DIMENSION(3) :: maxValues
    TYPE(ValueType), DIMENSION(3) :: States
  END TYPE  

  PUBLIC AllocateThermoState 
  
CONTAINS 

  SUBROUTINE AllocateThermoState( ThermoState, nValues )
     
    TYPE(ThermoStateType) :: ThermoState
    INTEGER, DIMENSION(3), INTENT(in) :: nValues

    INTEGER :: i
    
    DO i = 1, 3
      ALLOCATE( ThermoState % States(i) % Values(1:nValues(i)) ) 
    END DO 
 
  END SUBROUTINE AllocateThermoState


!  SUBROUTINE DeAllocateState( State, nValues )

!    DO i = 1, 3
!      DEALLOCATE( State % States(i) % Values(1:nValues(i)) ) 
!    END DO 

!  END SUBROUTINE DeAllocateState

END MODULE wlThermoStateModule 
