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
  PUBLIC CopyThermoState
 
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

  SUBROUTINE CopyThermoState( TS_target, TS_source )
    
    TYPE(ThermoStateType) :: TS_target, TS_source

    INTEGER :: i

    TS_target % Names           = TS_source % Names
    TS_target % Units           = TS_source % Units
    TS_target % nPoints         = TS_source % nPoints
    TS_target % LogInterp       = TS_source % LogInterp
    TS_target % minValues       = TS_source % minValues
    TS_target % maxValues       = TS_source % maxValues
    TS_target % Indices % iRho  = TS_source % Indices % iRho 
    TS_target % Indices % iT    = TS_source % Indices % iT 
    TS_target % Indices % iYe   = TS_source % Indices % iYe 

    DO i = 1, 3

      TS_target % States(i) % Values = TS_source % States(i) % Values     

    END DO

  END SUBROUTINE CopyThermoState

END MODULE wlThermoStateModule 
