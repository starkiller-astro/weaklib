MODULE wlGridModule

  USE wlKindModule, ONLY: dp

  IMPLICIT NONE
  PRIVATE

  TYPE, PUBLIC :: GridType
    CHARACTER(LEN=32) :: Name
    CHARACTER(LEN=32) :: Unit
    INTEGER  :: nPoints
    INTEGER  :: LogInterp = 0
    REAL(dp) :: minValue
    REAL(dp) :: maxValue
    REAL(dp) :: Zoom = 0.0d0
    REAL(dp), DIMENSION(:), ALLOCATABLE :: Values
    REAL(dp), DIMENSION(:), ALLOCATABLE :: Edge
    REAL(dp), DIMENSION(:), ALLOCATABLE :: Width
  END TYPE

  PUBLIC :: MakeLinearGrid
  PUBLIC :: MakeLogGrid
  PUBLIC :: MakeGeometricGrid
  PUBLIC :: AllocateGrid
  PUBLIC :: DeAllocateGrid
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

    ScaleFactor &
      = EXP( LOG( UpperBound / LowerBound ) / DBLE( nPoints - 1 ) )

    Grid(1) = LowerBound

    DO m = 2, nPoints - 1
      Grid(m) = Grid(m-1) * ScaleFactor
    END DO

    Grid(nPoints) = UpperBound

  END SUBROUTINE MakeLogGrid


  SUBROUTINE MakeGeometricGrid ( LowerBound, UpperBound, Zoom, nPoints, Grid, &
                                 Width, Edge)

    INTEGER,                        INTENT(in)  :: nPoints
    REAL(dp),                       INTENT(in)  :: LowerBound
    REAL(dp),                       INTENT(in)  :: UpperBound
    REAL(dp),                       INTENT(in)  :: Zoom
    REAl(dp), DIMENSION(nPoints),   INTENT(out) :: Grid 
    REAl(dp), DIMENSION(nPoints),   INTENT(out) :: Width 
    REAl(dp), DIMENSION(nPoints+1), INTENT(out) :: Edge 

    INTEGER  :: i

    ASSOCIATE &
      ( N   =>  nPoints, &
        xL  =>  LowerBound, &
        xR  =>  UpperBound )

    Width(1) = ( xR - xL ) * ( Zoom - 1.0_DP ) / ( Zoom**N - 1.0_DP )
    Grid (1) = xL + 0.5_DP * Width(1)
    Edge (1) = xL
    Edge (2) = xL + Width(1)
    DO i = 2, N
      Width(i)   = Width(i-1) * Zoom
      Grid (i)   = xL + SUM( Width(1:i-1) ) + 0.5_DP * Width(i)
      Edge (i+1) = xL + SUM( Width(1:i) )
    END DO

    END ASSOCIATE !-- N, etc.

  END SUBROUTINE MakeGeometricGrid

  
  SUBROUTINE AllocateGrid( Grid, nPoints )

    TYPE(GridType), INTENT(inout) :: Grid
    INTEGER, INTENT(in)                 :: nPoints

    Grid % nPoints = nPoints

    ALLOCATE( Grid % Values(nPoints) )
    ALLOCATE( Grid % Width(nPoints) )
    ALLOCATE( Grid % Edge(nPoints+1) )

  END SUBROUTINE AllocateGrid

  SUBROUTINE DeAllocateGrid( Grid )

    TYPE(GridType), INTENT(inout) :: Grid

    DEALLOCATE( Grid % Edge )
    DEALLOCATE( GRID % Width )
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
    WRITE(*,'(A6,A12,ES10.3E2)') &
      ' ', 'Zoom      = ', Grid % Zoom
    WRITE(*,'(A6,A12,I4.4)') &
      ' ', 'nPoints   = ', Grid % nPoints
    IF ( Grid % LogInterp == 1 ) THEN
      WRITE (*,'(A6,A27)') &
          ' ', 'Grid Logarithmically Spaced'
    ELSE IF ( Grid % Zoom  >=  1.0d0 ) THEN
      WRITE (*,'(A6,A25)') &
          ' ', 'Grid Geometrically Spaced'
    ELSE
      WRITE (*,'(A6,A20)') &
          ' ', 'Grid Linearly Spaced'
    END IF
    WRITE(*,*)

    DO i = 1, Grid % nPoints
      WRITE(*,'(A8,A6,I4.4,A4,ES10.3E2)') &
        ' ','Value(', i, ') = ', Grid % Values(i)
    END DO
    WRITE(*,*)

    IF ( Grid % Zoom  >=  1.0d0 ) THEN
      DO i = 1, size ( Grid % Width )
        WRITE(*,'(A8,A6,I4.4,A4,ES10.3E2)') &
          ' ','Width(', i, ') = ', Grid % Width(i)
      END DO
      WRITE(*,*)
      DO i = 1, size ( Grid % Edge )
        WRITE(*,'(A8,A5,I4.4,A4,ES10.3E2)') &
          ' ','Edge(', i, ') = ', Grid % Edge(i)
      END DO
    END IF

  END SUBROUTINE DescribeGrid


END MODULE wlGridModule
