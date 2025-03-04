MODULE wlThermoState4DModule
  
  USE wlKindModule, ONLY: dp

  implicit none
  PRIVATE

  TYPE, PUBLIC :: ValueTypeTS4D
    REAL(dp), ALLOCATABLE, DIMENSION(:) :: Values
  END TYPE

  TYPE, PUBLIC :: TS4DIDType
    INTEGER :: iRho
    INTEGER :: iT
    INTEGER :: iYe
    INTEGER :: iYm
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

  PUBLIC AllocateThermoState4D 
  PUBLIC DeAllocateThermoState4D 
  PUBLIC CopyThermoState4D
 
CONTAINS 

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

END MODULE wlThermoState4DModule 
