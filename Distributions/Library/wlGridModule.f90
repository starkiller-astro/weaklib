MODULE wlGridModule

  USE wlKindModule, ONLY: dp

  IMPLICIT NONE
  PRIVATE

  TYPE, PUBLIC :: GridType
    CHARACTER(LEN=32) :: Name
    CHARACTER(LEN=32) :: Unit
    INTEGER  :: nPoints
    INTEGER  :: LogInterp = 0
    REAL(dp) :: minValue = 0.0d0
    REAL(dp) :: maxValue = 0.0d0
    REAL(dp) :: minEdge  = 0.0d0 !-- only applies with Zoom
    REAL(dp) :: maxEdge  = 0.0d0 !-- only applies with Zoom
    REAL(dp) :: minWidth = 0.0d0 !-- only applies with Zoom
    REAL(dp) :: Zoom = 0.0d0     !-- non-zero value means geometric spacing
    REAL(dp), DIMENSION(:), ALLOCATABLE :: Values
    REAL(dp), DIMENSION(:), ALLOCATABLE :: Edge  !-- only applies with Zoom
    REAL(dp), DIMENSION(:), ALLOCATABLE :: Width !-- only applies with Zoom 
  END TYPE

  PUBLIC :: MakeLinearGrid
  PUBLIC :: MakeLogGrid
  PUBLIC :: MakeGeometricGrid
  PUBLIC :: AllocateGrid
  PUBLIC :: DeAllocateGrid
  PUBLIC :: DescribeGrid

  PRIVATE :: ComputeZoom
  PRIVATE :: ZeroZoom

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


  SUBROUTINE MakeGeometricGrid( LowerEdge, UpperEdge, MinWidth, nPoints, &
                                Grid, Width, Edge, Zoom, MinCenter, MaxCenter )

    INTEGER,                        INTENT(in)  :: nPoints
    REAL(dp),                       INTENT(in)  :: LowerEdge
    REAL(dp),                       INTENT(in)  :: UpperEdge
    REAL(dp),                       INTENT(in)  :: MinWidth
    REAL(dp), DIMENSION(nPoints),   INTENT(out) :: Grid 
    REAL(dp), DIMENSION(nPoints),   INTENT(out) :: Width 
    REAL(dp), DIMENSION(nPoints+1), INTENT(out) :: Edge
    REAL(dp),                       INTENT(out) :: Zoom
    REAL(dp),                       INTENT(out) :: MinCenter
    REAL(dp),                       INTENT(out) :: MaxCenter

    INTEGER  :: i

    ASSOCIATE &
      ( N   =>  nPoints, &
        xL  =>  LowerEdge, &
        xR  =>  UpperEdge )

    call ComputeZoom( LowerEdge, UpperEdge, MinWidth, nPoints, Zoom )

    Width(1) = ( xR - xL ) * ( Zoom - 1.0_DP ) / ( Zoom**N - 1.0_DP )
    Grid (1) = xL + 0.5_DP * Width(1)
    Edge (1) = xL
    Edge (2) = xL + Width(1)
    DO i = 2, N
      Width(i)   = Width(i-1) * Zoom
      Grid (i)   = xL + SUM( Width(1:i-1) ) + 0.5_DP * Width(i)
      Edge (i+1) = xL + SUM( Width(1:i) )
    END DO
    MinCenter = Grid(1)
    MaxCenter = Grid(N)

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
      ' ', 'Min Width = ', Grid % minWidth
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

  SUBROUTINE ComputeZoom( LowerEdge, UpperEdge, MinWidth, nPoints, Zoom )

    REAL(dp), INTENT(in)  :: LowerEdge
    REAL(dp), INTENT(in)  :: UpperEdge
    REAL(dp), INTENT(in)  :: MinWidth
    INTEGER,  INTENT(in)  :: nPoints
    REAL(dp), INTENT(out) :: Zoom

    INTEGER :: i
    REAL(dp) :: a, b, c
    REAL(dp) :: fa, fb, fc

    a  = 1.000001_dp
    b  = 2.0_dp
    fa = ZeroZoom ( a, LowerEdge, UpperEdge, MinWidth, nPoints )
    fb = ZeroZoom ( b, LowerEdge, UpperEdge, MinWidth, nPoints )

    !-- Bisection
    do i = 1, 100
      c  = 0.5_dp * ( a + b )
      if ( ( b - a ) / c  <  1.0e-12_dp ) then
        Zoom = c
        return
      end if
      fc = ZeroZoom ( c, LowerEdge, UpperEdge, MinWidth, nPoints )
      if ( sign ( 1.0_dp, fc )  ==  sign ( 1.0_dp, fa ) ) then
        a  = c
        fa = fc
      else
        b  = c
        fb = fc
      end if
    end do !-- i

  END SUBROUTINE ComputeZoom

  FUNCTION ZeroZoom &
             ( Zoom, LowerEdge, UpperEdge, MinWidth, nPoints ) &
             result ( ZZ )

    REAL(dp), INTENT(in) :: Zoom
    REAL(dp), INTENT(in) :: LowerEdge
    REAL(dp), INTENT(in) :: UpperEdge
    REAL(dp), INTENT(in) :: MinWidth
    INTEGER,  INTENT(in) :: nPoints
    REAL(dp)             :: ZZ

    ZZ  =  ( UpperEdge - LowerEdge ) * ( Zoom - 1.0_dp ) &
              /  ( Zoom ** nPoints  -  1.0_dp ) &
            -  MinWidth

  END FUNCTION ZeroZoom


END MODULE wlGridModule
