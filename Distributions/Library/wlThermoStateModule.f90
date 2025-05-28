MODULE wlThermoStateModule
  
  USE wlKindModule, ONLY: dp

  implicit none
  PRIVATE

  TYPE, PUBLIC :: ValueTypeTS
    REAL(dp), ALLOCATABLE, DIMENSION(:) :: Values
  END TYPE

  TYPE, PUBLIC :: ValueTypeTS4D
    REAL(dp), ALLOCATABLE, DIMENSION(:) :: Values
  END TYPE

  TYPE, PUBLIC :: TSIDType
    INTEGER :: iRho
    INTEGER :: iT
    INTEGER :: iYe
  END TYPE
  
  TYPE, PUBLIC :: TS4DIDType
    INTEGER :: iRho
    INTEGER :: iT
    INTEGER :: iYe
    INTEGER :: iYm
  END TYPE

  TYPE, PUBLIC :: ThermoStateType
    CHARACTER(LEN=32), DIMENSION(3) :: Names
    CHARACTER(LEN=32), DIMENSION(3) :: Units
    INTEGER, DIMENSION(3) :: nPoints
    INTEGER, DIMENSION(3) :: LogInterp
    REAL(dp), DIMENSION(3) :: minValues
    REAL(dp), DIMENSION(3) :: maxValues
    TYPE(ValueTypeTS), DIMENSION(3) :: States
    TYPE(TSIDType) :: Indices
  END TYPE  

  TYPE, PUBLIC :: ThermoState4DType
    CHARACTER(LEN=32), DIMENSION(4) :: Names
    CHARACTER(LEN=32), DIMENSION(4) :: Units
    INTEGER, DIMENSION(4) :: nPoints
    INTEGER, DIMENSION(4) :: LogInterp
    REAL(dp), DIMENSION(4) :: minValues
    REAL(dp), DIMENSION(4) :: maxValues
    TYPE(ValueTypeTS4D), DIMENSION(4) :: States
    TYPE(TS4DIDType) :: Indices
  END TYPE  

  PUBLIC AllocateThermoState 
  PUBLIC DeAllocateThermoState 
  PUBLIC CopyThermoState
 
  PUBLIC AllocateThermoState4D 
  PUBLIC DeAllocateThermoState4D 
  PUBLIC CopyThermoState4D
 
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

  SUBROUTINE AllocateThermoState4D( TS, nPoints )
     
    TYPE(ThermoState4DType) :: TS
    INTEGER, DIMENSION(4), INTENT(in) :: nPoints

    INTEGER :: i
   
    TS % nPoints = nPoints    
 
    DO i = 1, 4
      ALLOCATE( TS % States(i) % Values(1:TS % nPoints(i)) ) 
    END DO 
 
  END SUBROUTINE AllocateThermoState4D

  SUBROUTINE DeAllocateThermoState4D( TS )

    TYPE(ThermoState4DType) :: TS

    INTEGER :: i
    
    DO i = 1, 4
      DEALLOCATE( TS % States(i) % Values )
    END DO 

  END SUBROUTINE DeAllocateThermoState4D

  SUBROUTINE CopyThermoState4D( TS_target, TS_source )
    
    TYPE(ThermoState4DType) :: TS_target, TS_source

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
    TS_target % Indices % iYm   = TS_source % Indices % iYm 

    DO i = 1, 4

      TS_target % States(i) % Values = TS_source % States(i) % Values     

    END DO

  END SUBROUTINE CopyThermoState4D

END MODULE wlThermoStateModule 
