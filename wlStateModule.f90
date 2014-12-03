MODULE wlStateModule
  
!  USE wlKindModule, ONLY: dp

  implicit none
  PRIVATE

  TYPE :: ValueType
    REAL(8), ALLOCATABLE, DIMENSION(:) :: Values
  END TYPE
  
  TYPE, PUBLIC :: StateType
    CHARACTER(LEN=32), DIMENSION(3) :: Names
    INTEGER, DIMENSION(3) :: nValues
    REAL(8), DIMENSION(3) :: minValues
    REAL(8), DIMENSION(3) :: maxValues
    TYPE(ValueType), DIMENSION(3) :: States
  END TYPE  
!  TYPE :: ValueType
!    REAL(dp), ALLOCATABLE, DIMENSION(:) :: Values
!  END TYPE

  PUBLIC AllocateState 
  
CONTAINS 

  SUBROUTINE AllocateState( State, nValues )
     
    TYPE(StateType) :: State
    INTEGER, DIMENSION(3), INTENT(in) :: nValues

    INTEGER :: i
    
    DO i = 1, 3
      ALLOCATE( State % States(i) % Values(1:nValues(i)) ) 
    END DO 
 
  END SUBROUTINE AllocateState


!  SUBROUTINE DeAllocateState( State, nValues )

!    DO i = 1, 3
!      DEALLOCATE( State % States(i) % Values(1:nValues(i)) ) 
!    END DO 

!  END SUBROUTINE DeAllocateState

END MODULE wlStateModule 
