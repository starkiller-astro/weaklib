MODULE wlGridModule

  USE wlKindModule, ONLY: dp
  USE wlThermoStateModule

  implicit none
  PRIVATE

  PUBLIC MakeLinearGrid
  PUBLIC MakeLogGrid
  PUBLIC WriteGrid 

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

    DO m=2,nPoints - 1
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

    DO m=2,nPoints - 1
      Grid(m) = Grid(m-1) * ScaleFactor
    END DO

    Grid(nPoints) = UpperBound

  END SUBROUTINE MakeLogGrid

  SUBROUTINE WriteGrid( rank, griddims, Grid ) 
                          ! Do we want one rho t ye Grid, or 3 subgrids?
    

    INTEGER                                   :: hdferr
    INTEGER                                   :: i 
    INTEGER, INTENT(in)                       :: rank
    INTEGER, DIMENSION(rank), INTENT(in)      :: griddims          
    TYPE(ThermoStateType), INTENT(in)         :: ThermoState1    
    INTEGER(HID_T) :: file_id
    INTEGER(HID_T) :: dset_id
    INTEGER(HID_T) :: space_id
    CHARACTER(LEN=10), PARAMETER  :: filename = "wltable.h5"
    REAL(dp), ALLOCATABLE, DIMENSION(:) :: buffer
    CALL h5open_f(hdferr)
    CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file, hdferr)
    CALL h5screate_simple_f(rank, griddims, space_id, hdferr)

     DO i = 1, rank
       CALL h5dcreate_f(file, TRIM( ThermoState1 % Names(i) ), H5T_NATIVE_DOUBLE, space_id, dset_id, hdferr) 
         !DO j = 1, griddims(i)
         !  buffer(j) = Thermostate1 % States(j) % Values(j)
         !END DO
       CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, Thermostate1 % States(i) % Values(:), griddims(i), hdferr)
       CALL h5dclose_f(dset_id, hdferr)
     END DO

  END SUBROUTINE WriteGrid

END MODULE wlGridModule
