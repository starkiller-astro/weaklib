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

  TYPE, PUBLIC :: EnergyGridType
    CHARACTER(LEN=32) :: Name
    CHARACTER(LEN=32) :: Unit
    INTEGER  :: nPoints
    INTEGER  :: nFaces
    INTEGER  :: LogInterp
    REAL(dp) :: minValueCenters, minValueFaces
    REAL(dp) :: maxValueCenters, maxValueFaces
    REAL(dp), DIMENSION(:), ALLOCATABLE :: Cell_centers
    REAL(dp), DIMENSION(:), ALLOCATABLE :: Cell_faces
    REAL(dp), DIMENSION(:), ALLOCATABLE :: dE
  END TYPE

  PUBLIC :: AllocateGrid
  PUBLIC :: DeAllocateGrid
  PUBLIC :: DescribeGrid

  PUBLIC :: MakeLinearGrid
  PUBLIC :: MakeLogGrid
  PUBLIC :: MakeLogGrid_Generic
  PUBLIC :: MakeLogGrid_Centers_Faces
  PUBLIC :: AllocateGenericGrid
  PUBLIC :: AllocateEnergyGrid
  PUBLIC :: DeallocateGenericGrid
  PUBLIC :: DeallocateEnergyGrid
  PUBLIC :: DescribeGenericGrid
  PUBLIC :: DescribeEnergyGrid

  INTERFACE MakeLogGrid
    MODULE PROCEDURE MakeLogGrid_Generic
    MODULE PROCEDURE MakeLogGrid_Centers_Faces
  END INTERFACE MakeLogGrid

  INTERFACE AllocateGrid
    MODULE PROCEDURE AllocateGenericGrid
    MODULE PROCEDURE AllocateEnergyGrid
  END INTERFACE AllocateGrid

  INTERFACE DeAllocateGrid
    MODULE PROCEDURE DeAllocateGenericGrid
    MODULE PROCEDURE DeAllocateEnergyGrid
  END INTERFACE DeAllocateGrid

  INTERFACE DescribeGrid
    MODULE PROCEDURE DescribeGenericGrid
    MODULE PROCEDURE DescribeEnergyGrid
  END INTERFACE DescribeGrid

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


  SUBROUTINE MakeLogGrid_Generic( LowerBound, UpperBound, nPoints, Grid )

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

  END SUBROUTINE MakeLogGrid_Generic

  SUBROUTINE MakeLogGrid_Centers_Faces( LowerBound, UpperBound, nPoints, & 
                                  Cell_centers, Cell_faces, dE )

    Integer, INTENT(in)      :: nPoints     !number of cells
    REAL(dp), INTENT(in)     :: LowerBound  !cell center of first cell
    Real(dp), INTENT(in)     :: UpperBound  !cell center of last cell

    REAl(dp), DIMENSION(nPoints),   INTENT(out) :: Cell_centers

    !lower face of first cell is at 0
    REAl(dp), DIMENSION(nPoints+1), INTENT(out) :: Cell_faces 
    REAl(dp), DIMENSION(nPoints),   INTENT(out) :: dE

    Integer  :: m
    REAL(dp) :: ScaleFactor

    IF ( nPoints == 1 ) THEN
      Cell_centers(1) = LowerBound
      Cell_faces(1)   = 0.0d0
      Cell_faces(2)   = Lowerbound + LowerBound/2.0d0
      dE(1)           = Cell_faces(2)
      RETURN
    END IF

    ScaleFactor &
      = EXP( LOG( UpperBound / LowerBound ) / DBLE( nPoints - 1 ) )

    Cell_centers(1) = LowerBound

    DO m = 2, nPoints - 1
      Cell_centers(m) = Cell_centers(m-1) * ScaleFactor
    END DO

    Cell_faces(1) = 0.0d0
    Cell_faces(2) = (Cell_centers(1) + Cell_centers(2))/2.0d0

    DO m = 3, nPoints + 1
      Cell_faces(m) = Cell_faces(m-1) * ScaleFactor
    ENDDO

    DO m = 1, nPoints
      dE(m) = Cell_faces(m+1) - Cell_faces(m)
    ENDDO 

    Cell_centers(nPoints) = UpperBound

  END SUBROUTINE MakeLogGrid_Centers_Faces

  SUBROUTINE AllocateGenericGrid( Grid, nPoints )

    TYPE(GridType), INTENT(inout) :: Grid
    INTEGER, INTENT(in)                 :: nPoints

    Grid % nPoints = nPoints

    ALLOCATE( Grid % Values(nPoints) )

  END SUBROUTINE AllocateGenericGrid

  SUBROUTINE AllocateEnergyGrid( Grid, nPoints )

    TYPE(EnergyGridType), INTENT(inout) :: Grid
    INTEGER, INTENT(in)                 :: nPoints

    Grid % nPoints = nPoints
    Grid % nFaces = nPoints + 1

    ALLOCATE( Grid % Cell_centers(nPoints  ) )
    ALLOCATE( Grid % Cell_faces  (nPoints+1) )
    ALLOCATE( Grid % dE          (nPoints  ) )

  END SUBROUTINE AllocateEnergyGrid

  SUBROUTINE DeAllocateGenericGrid( Grid )

    TYPE(GridType), INTENT(inout) :: Grid

    DEALLOCATE( Grid % Values )

  END SUBROUTINE DeAllocateGenericGrid

  SUBROUTINE DeAllocateEnergyGrid( Grid )

    TYPE(EnergyGridType), INTENT(inout) :: Grid

    DEALLOCATE( Grid % Cell_centers )
    DEALLOCATE( Grid % Cell_faces   )
    DEALLOCATE( Grid % dE           )

  END SUBROUTINE DeAllocateEnergyGrid

  SUBROUTINE DescribeGenericGrid( Grid )

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

  END SUBROUTINE DescribeGenericGrid

  SUBROUTINE DescribeEnergyGrid( Grid )

    TYPE(EnergyGridType), INTENT(in) :: Grid

    INTEGER :: i

    WRITE(*,*)
    WRITE(*,'(A4,A)') '', 'Grid:'
    WRITE(*,*)

    WRITE(*,'(A6,A20,A)') &
      ' ', 'Name              = ', TRIM( Grid % Name )
    WRITE(*,'(A6,A20,A)') &
      ' ', 'Unit              = ', TRIM( Grid % Unit )
    WRITE(*,'(A6,A20,ES10.3E2)') &
      ' ', 'Min Value centers = ', Grid % minValueCenters
    WRITE(*,'(A6,A20,ES10.3E2)') &
      ' ', 'Max Value centers = ', Grid % maxValueCenters
    WRITE(*,'(A6,A20,ES10.3E2)') &
      ' ', 'Min Value faces   = ', Grid % minValueFaces
    WRITE(*,'(A6,A20,ES10.3E2)') &
      ' ', 'Max Value faces   = ', Grid % maxValueFaces
    WRITE(*,'(A6,A20,I4.4)') &
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
      WRITE(*,'(A8,A10,I4.4,A4,ES10.3E2)') &
        ' ','Lower face(', i, ') = ', Grid % Cell_faces(i)
      WRITE(*,'(A8,A10,I4.4,A4,ES10.3E2)') &
        ' ','Cell center(', i, ') = ', Grid % Cell_centers(i)
      WRITE(*,'(A8,A10,I4.4,A4,ES10.3E2)') &
        ' ','Upper face(', i, ') = ', Grid % Cell_faces(i+1)
    END DO

  END SUBROUTINE DescribeEnergyGrid


END MODULE wlGridModule
