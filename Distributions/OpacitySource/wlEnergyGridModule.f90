MODULE wlEnergyGridModule

  USE wlKindModule, ONLY: dp

  IMPLICIT NONE
  PRIVATE

  TYPE, PUBLIC :: EnergyGridType
    CHARACTER(LEN=32) :: Name
    CHARACTER(LEN=32) :: Unit
    INTEGER  :: nPoints
    INTEGER  :: LogInterp
    REAL(dp) :: minValue
    REAL(dp) :: maxValue
    REAL(dp), DIMENSION(:), ALLOCATABLE :: Values
  END TYPE

  PUBLIC :: AllocateEnergyGrid
  PUBLIC :: DeallocateEnergyGrid
  PUBLIC :: DescribeEnergyGrid

CONTAINS


  SUBROUTINE AllocateEnergyGrid( EnergyGrid, nPoints )

    TYPE(EnergyGridType), INTENT(inout) :: EnergyGrid
    INTEGER, INTENT(in)                 :: nPoints

    EnergyGrid % nPoints = nPoints

    ALLOCATE( EnergyGrid % Values(nPoints) )

  END SUBROUTINE AllocateEnergyGrid


  SUBROUTINE DeAllocateEnergyGrid( EnergyGrid )
    
    TYPE(EnergyGridType) :: EnergyGrid
 
    DEALLOCATE( EnergyGrid % Values )  

  END SUBROUTINE DeAllocateEnergyGrid


  SUBROUTINE DescribeEnergyGrid( EnergyGrid )

    TYPE(EnergyGridType), INTENT(in) :: EnergyGrid

    INTEGER :: i

    WRITE(*,*)
    WRITE(*,'(A4,A)') '', 'Energy Grid:'
    WRITE(*,*)

    WRITE(*,'(A6,A12,A)') &
      ' ', 'Name      = ', TRIM( EnergyGrid % Name )
    WRITE(*,'(A6,A12,A)') &
      ' ', 'Unit      = ', TRIM( EnergyGrid % Unit )
    WRITE(*,'(A6,A12,ES10.4E2)') &
      ' ', 'Min Value = ', EnergyGrid % minValue
    WRITE(*,'(A6,A12,ES10.4E2)') &
      ' ', 'Max Value = ', EnergyGrid % maxValue
    WRITE(*,'(A6,A12,I4.4)') &
      ' ', 'nPoints   = ', EnergyGrid % nPoints
    IF ( EnergyGrid % LogInterp == 1 ) THEN
      WRITE (*,'(A6,A27)') &
          ' ', 'Grid Logarithmically Spaced'
    ELSE
      WRITE (*,'(A6,A20)') &
          ' ', 'Grid Linearly Spaced'
    END IF
    WRITE(*,*)

    DO i = 1, EnergyGrid % nPoints
      WRITE(*,'(A8,A6,I4.4,A4,ES10.4E2)') &
        ' ','Value(', i, ') = ', EnergyGrid % Values(i)
    END DO

  END SUBROUTINE DescribeEnergyGrid


END MODULE wlEnergyGridModule
