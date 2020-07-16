PROGRAM wlTestInterpolation

  USE wlKindModule, ONLY: dp
  USE wlInterpolationModule, ONLY: &
    LogInterpolateSingleVariable, &
    LogInterpolateSingleVariable_2D2D_Custom
  USE wlOpacityFieldsModule, ONLY: &
    iJi0, iJii0, iJi1, iJii1
  USE wlOpacityTableModule, ONLY: &
    OpacityTableType, &
    DeAllocateOpacityTable
  USE wlIOModuleHDF, ONLY: &
    InitializeHDF, &
    FinalizeHDF, &
    OpenFileHDF, &
    CloseFileHDF, &
    WriteHDF, &
    OpenGroupHDF, &
    CloseGroupHDF
  USE wlOpacityTableIOModuleHDF, ONLY: &
    ReadOpacityTableHDF
  USE wlGridModule, ONLY: &
    GridType, &
    AllocateGrid, &
    DescribeGrid, &
    MakeLogGrid
  USE wlExtPhysicalConstantsModule, ONLY: kMeV, ca, cv 
  USE wlExtNumericalModule, ONLY: pi, half, twpi, zero
  USE HDF5

  USE, INTRINSIC :: iso_fortran_env, only : stdin=>input_unit, &
                                            stdout=>output_unit, &
                                            stderr=>error_unit

  IMPLICIT NONE

  !--------- parameters for creating grid to test interpolation ------------
  INTEGER, PARAMETER     :: nx_lo_res = 15 !create two different base arrays to check 
  INTEGER, PARAMETER     :: nx_hi_res = 30 !for convergence of the interpolation routines
  REAL(dp)               :: x_min = 0.0d0
  REAL(dp)               :: x_max = 1.0d2

  REAL(dp), DIMENSION(nx_lo_res,nx_lo_res,nx_lo_res,nx_lo_res,nx_lo_res) :: linear_lo_res
  REAL(dp), DIMENSION(nx_hi_res,nx_hi_res,nx_hi_res,nx_hi_res,nx_hi_res) :: linear_hi_res

  REAL(dp), DIMENSION(nx_lo_res,nx_lo_res,nx_lo_res,nx_lo_res,nx_lo_res) :: general_lo_res
  REAL(dp), DIMENSION(nx_hi_res,nx_hi_res,nx_hi_res,nx_hi_res,nx_hi_res) :: general_hi_res

  REAL(dp), DIMENSION(nx_lo_res) :: axis_lo_res
  REAL(dp), DIMENSION(nx_hi_res) :: axis_hi_res

  INTEGER :: ii, jj, kk, ll, mm

  REAL(dp) :: dx_lo, dx_hi

  dx_lo = (x_max-x_min) / nx_lo_res
  dx_hi = (x_max-x_min) / nx_hi_res


  DO ii = 1, nx_lo_res

    axis_lo_res(ii) = xmin + (ii-1)*dx_lo

    write(stdout,*) 'axis_lo_res(ii)', ii, axis_lo_res(ii), dx_lo

  END DO 

  DO ii = 1, nx_hi_res

    axis_hi_res(ii) = xmin + (ii-1)*dx_hi

  END DO 

END PROGRAM wlTestInterpolation
