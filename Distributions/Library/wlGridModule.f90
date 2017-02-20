MODULE wlGridModule

  USE wlKindModule, ONLY: dp

  IMPLICIT NONE
  PRIVATE

  TYPE, PUBLIC :: GridType
    CHARACTER(LEN=32) :: Name
    CHARACTER(LEN=32) :: Unit
    INTEGER  :: nPoints
    INTEGER  :: LogInterp
    REAL(dp) :: minValue
    REAL(dp) :: maxValue
    REAL(dp), DIMENSION(:), ALLOCATABLE :: Values
  END TYPE

  PUBLIC :: MakeLinearGrid
  PUBLIC :: MakeLogGrid
  PUBLIC :: AllocateGrid
  PUBLIC :: DeallocateGrid
  PUBLIC :: DescribeGrid

CONTAINS

  SUBROUTINE MakeLinearGrid( LowerBound, UpperBound, nPoints, Grid )

    INTEGER,                      INTENT(in)  :: nPoints
    REAL(dp),                     INTENT(in)  :: LowerBound
    REAL(dp),                     INTENT(in)  :: UpperBound
    REAl(dp), DIMENSION(nPoints), INTENT(out) :: Grid 

    INTEGER  :: i
    REAL(dp) :: BinWidth

    IF ( nPoints == 1 ) THEN
      Grid(1) = LowerBound
      RETURN
    END IF

    BinWidth = ( UpperBound - LowerBound ) / DBLE( nPoints - 1 )

    DO i = 1, nPoints
      Grid(i) = LowerBound + DBLE( i - 1 ) * BinWidth
    END DO
 
  END SUBROUTINE MakeLinearGrid


  SUBROUTINE MakeLogGrid( LowerBound, UpperBound, nPoints, Grid )

    Integer, INTENT(in)      :: nPoints
    REAL(dp), INTENT(in)     :: LowerBound
    Real(dp), INTENT(in)     :: UpperBound

    REAl(dp), DIMENSION(nPoints), INTENT(out)   :: Grid 

    Integer  :: m
    REAL(dp) :: ScaleFactor

    IF ( nPoints == 1 ) THEN
      Grid(1) = LowerBound
      RETURN
    END IF

    ScaleFactor =  EXP( LOG( UpperBound / LowerBound ) / DBLE( nPoints - 1 ) )

    Grid(1) = LowerBound

    DO m = 2, nPoints - 1
      Grid(m) = Grid(m-1) * ScaleFactor
    END DO

    Grid(nPoints) = UpperBound

  END SUBROUTINE MakeLogGrid

  SUBROUTINE AllocateGrid( Grid, nPoints )

    TYPE(GridType), INTENT(inout) :: Grid
    INTEGER, INTENT(in)                 :: nPoints

    Grid % nPoints = nPoints

    ALLOCATE( Grid % Values(nPoints) )

  END SUBROUTINE AllocateGrid

  SUBROUTINE DeAllocateGrid( Grid )

    TYPE(GridType), INTENT(inout) :: Grid

    DEALLOCATE( Grid % Values )

  END SUBROUTINE DeAllocateGrid

  SUBROUTINE DescribeGrid( Grid )

    TYPE(GridType), INTENT(in) :: Grid

    INTEGER :: i

    WRITE(*,*)
    WRITE(*,'(A4,A)') '', 'Grid:'
    WRITE(*,*)

    WRITE(*,'(A6,A12,A)') &
      ' ', 'Name      = ', TRIM( Grid % Name )
    WRITE(*,'(A6,A12,A)') &
      ' ', 'Unit      = ', TRIM( Grid % Unit )
    WRITE(*,'(A6,A12,ES10.3E2)') &
      ' ', 'Min Value = ', Grid % minValue
    WRITE(*,'(A6,A12,ES10.3E2)') &
      ' ', 'Max Value = ', Grid % maxValue
    WRITE(*,'(A6,A12,I4.4)') &
      ' ', 'nPoints   = ', Grid % nPoints
    IF ( Grid % LogInterp == 1 ) THEN
      WRITE (*,'(A6,A27)') &
          ' ', 'Grid Logarithmically Spaced'
    ELSE
      WRITE (*,'(A6,A20)') &
          ' ', 'Grid Linearly Spaced'
    END IF
    WRITE(*,*)

    DO i = 1, Grid % nPoints
      WRITE(*,'(A8,A6,I4.4,A4,ES10.3E2)') &
        ' ','Value(', i, ') = ', Grid % Values(i)
    END DO

  END SUBROUTINE DescribeGrid


END MODULE wlGridModule
