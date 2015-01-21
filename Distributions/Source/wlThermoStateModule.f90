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
  PUBLIC DeAllocateThermoState 
  
CONTAINS 

  SUBROUTINE AllocateThermoState( ThermoState, nValues )
     
    TYPE(ThermoStateType) :: ThermoState
    INTEGER, DIMENSION(3), INTENT(in) :: nValues

    INTEGER :: i
   
    ThermoState % nValues = nValues    
 
    DO i = 1, 3
      ALLOCATE( ThermoState % States(i) % Values(1:ThermoState % nValues(i)) ) 
    END DO 
 
  END SUBROUTINE AllocateThermoState


  SUBROUTINE DeAllocateThermoState( ThermoState )

    TYPE(ThermoStateType) :: ThermoState

    INTEGER :: i
    
    DO i = 1, 3
      DEALLOCATE( ThermoState % States(i) % Values )
    END DO 

  END SUBROUTINE DeAllocateThermoState

END MODULE wlThermoStateModule 
