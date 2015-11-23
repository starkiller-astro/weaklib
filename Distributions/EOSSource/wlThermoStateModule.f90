MODULE wlThermoStateModule
  
  USE wlKindModule, ONLY: dp

  implicit none
  PRIVATE

  TYPE :: ValueType
    REAL(dp), ALLOCATABLE, DIMENSION(:) :: Values
  END TYPE

  TYPE :: TSIDType
    INTEGER :: iRho
    INTEGER :: iT
    INTEGER :: iYe
  END TYPE
  
  TYPE, PUBLIC :: ThermoStateType
    CHARACTER(LEN=32), DIMENSION(3) :: Names
    CHARACTER(LEN=32), DIMENSION(3) :: Units
    INTEGER, DIMENSION(3) :: nPoints
    INTEGER, DIMENSION(3) :: LogInterp
    REAL(dp), DIMENSION(3) :: minValues
    REAL(dp), DIMENSION(3) :: maxValues
    TYPE(ValueType), DIMENSION(3) :: States
    TYPE(TSIDType) :: Indices
  END TYPE  

  PUBLIC AllocateThermoState 
  PUBLIC DeAllocateThermoState 
  
CONTAINS 

  SUBROUTINE AllocateThermoState( TS, nPoints )
     
    TYPE(ThermoStateType) :: TS
    INTEGER, DIMENSION(3), INTENT(in) :: nPoints

    INTEGER :: i
   
    TS % nPoints = nPoints    
 
    DO i = 1, 3
      ALLOCATE( TS % States(i) % Values(1:TS % nPoints(i)) ) 
    END DO 
 
  END SUBROUTINE AllocateThermoState

  SUBROUTINE DeAllocateThermoState( TS )

    TYPE(ThermoStateType) :: TS

    INTEGER :: i
    
    DO i = 1, 3
      DEALLOCATE( TS % States(i) % Values )
    END DO 

  END SUBROUTINE DeAllocateThermoState

END MODULE wlThermoStateModule 
