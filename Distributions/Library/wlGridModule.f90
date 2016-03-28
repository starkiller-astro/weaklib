MODULE wlGridModule

  USE wlKindModule, ONLY: dp

  implicit none
  PRIVATE

  PUBLIC MakeLinearGrid
  PUBLIC MakeLogGrid

  TYPE, PUBLIC :: EnergyGridType
    CHARACTER(LEN=32) :: Names='Energy'
    CHARACTER(LEN=32) :: Units='MeV'
    INTEGER  :: nPointsE
    REAL(dp) :: minValue
    REAL(dp) :: maxValue
!    REAL(dp) :: minWidth
!    REAL(dp) :: alpha  ! constant in geometric series
    REAL(dp), ALLOCATABLE, DIMENSION(:) :: Values
  END TYPE


CONTAINS 

  SUBROUTINE MakeLinearGrid( LowerBound, UpperBound, nPoints, Grid )

    INTEGER, INTENT(in)               :: nPoints  
    REAL(dp), INTENT(in)              :: LowerBound
    REAL(dp), INTENT(in)              :: UpperBound 
   
    REAl(dp), DIMENSION(nPoints), INTENT(out)   :: Grid 

    Integer  :: m
    REAL(dp) :: BinWidth

    BinWidth = ( UpperBound - LowerBound ) / DBLE( nPoints - 1 )

    Grid(1) = LowerBound
 
    DO m = 2, nPoints - 1
      Grid(m) = Grid(m-1) + BinWidth
    END DO

    Grid(nPoints) = UpperBound

  END SUBROUTINE MakeLinearGrid

  SUBROUTINE MakeLogGrid( LowerBound, UpperBound, nPoints, Grid )

    Integer, INTENT(in)      :: nPoints          
    REAL(dp), INTENT(in)     :: LowerBound
    Real(dp), INTENT(in)     :: UpperBound

    REAl(dp), DIMENSION(nPoints), INTENT(out)   :: Grid 

    Integer  :: m
    REAL(dp) :: ScaleFactor

    ScaleFactor =  EXP( LOG( UpperBound / LowerBound ) / DBLE( nPoints - 1 ) )

    Grid(1) = LowerBound

    DO m = 2, nPoints - 1
      Grid(m) = Grid(m-1) * ScaleFactor
    END DO

    Grid(nPoints) = UpperBound

  END SUBROUTINE MakeLogGrid

END MODULE wlGridModule
