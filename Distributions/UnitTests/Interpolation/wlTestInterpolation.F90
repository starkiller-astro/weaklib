PROGRAM wlTestInterpolation

  USE wlKindModule, ONLY: dp
  USE wlInterpolationModule, ONLY: &
    LinearInterp_Array_Point, &
    GetIndexAndDelta

  USE, INTRINSIC :: iso_fortran_env, only : stdin=>input_unit, &
                                            stdout=>output_unit, &
                                            stderr=>error_unit

  IMPLICIT NONE

  !--------- parameters for creating grid to test interpolation ------------
  INTEGER, PARAMETER     :: nx_lo_res = 15 !create two different base arrays to check 
  INTEGER, PARAMETER     :: nx_hi_res = 30 !for convergence of the interpolation routines
  
  !--------- extents of the test array
  REAL(dp), PARAMETER    :: x1_min = 0.0d0
  REAL(dp), PARAMETER    :: x1_max = 1.0d0
  REAL(dp), PARAMETER    :: x2_min = -1.0d0
  REAL(dp), PARAMETER    :: x2_max = 2.0d0
  REAL(dp), PARAMETER    :: x3_min = 2.0d0
  REAL(dp), PARAMETER    :: x3_max = 3.0d0
  REAL(dp), PARAMETER    :: x4_min = -2.5d0
  REAL(dp), PARAMETER    :: x4_max = 2.7d0
  REAL(dp), PARAMETER    :: x5_min = 0.0d0
  REAL(dp), PARAMETER    :: x5_max = 10.0d0

  !--------- the actual points to interpolate
  REAL(dp), PARAMETER    :: x1 = 0.1d0
  REAL(dp), PARAMETER    :: x2 = -0.34d0
  REAL(dp), PARAMETER    :: x3 = 2.47d0
  REAL(dp), PARAMETER    :: x4 = -0.793d0
  REAL(dp), PARAMETER    :: x5 = 9.375d0

  REAL(dp), DIMENSION(nx_lo_res,nx_lo_res,nx_lo_res,nx_lo_res,nx_lo_res) :: linear_data

  REAL(dp), DIMENSION(nx_lo_res,nx_lo_res,nx_lo_res,nx_lo_res,nx_lo_res) :: non_linear_data_lo_res
  REAL(dp), DIMENSION(nx_hi_res,nx_hi_res,nx_hi_res,nx_hi_res,nx_hi_res) :: non_linear_data_hi_res

  REAL(dp), DIMENSION(nx_lo_res) :: axis1_lo_res
  REAL(dp), DIMENSION(nx_lo_res) :: axis2_lo_res
  REAL(dp), DIMENSION(nx_lo_res) :: axis3_lo_res
  REAL(dp), DIMENSION(nx_lo_res) :: axis4_lo_res
  REAL(dp), DIMENSION(nx_lo_res) :: axis5_lo_res

  REAL(dp), DIMENSION(nx_hi_res) :: axis1_hi_res
  REAL(dp), DIMENSION(nx_hi_res) :: axis2_hi_res
  REAL(dp), DIMENSION(nx_hi_res) :: axis3_hi_res
  REAL(dp), DIMENSION(nx_hi_res) :: axis4_hi_res
  REAL(dp), DIMENSION(nx_hi_res) :: axis5_hi_res

  INTEGER :: ii, jj, kk, ll, mm

  REAL(dp) :: dx1_lo, dx1_hi
  REAL(dp) :: dx2_lo, dx2_hi
  REAL(dp) :: dx3_lo, dx3_hi
  REAL(dp) :: dx4_lo, dx4_hi
  REAL(dp) :: dx5_lo, dx5_hi

  REAL(dp) :: interp_result_linear
  REAL(dp) :: interp_result_non_linear_lo
  REAL(dp) :: interp_result_non_linear_hi

  INTEGER  :: ix1, ix2, ix3, ix4, ix5
  REAL(dp) :: dx1, dx2, dx3, dx4, dx5

  REAL(dp) :: linear_value
  REAL(dp) :: non_linear_value

  dx1_lo = (x1_max-x1_min) / DBLE(nx_lo_res-1)
  dx2_lo = (x2_max-x2_min) / DBLE(nx_lo_res-1)
  dx3_lo = (x3_max-x3_min) / DBLE(nx_lo_res-1)
  dx4_lo = (x4_max-x4_min) / DBLE(nx_lo_res-1)
  dx5_lo = (x5_max-x5_min) / DBLE(nx_lo_res-1)

  dx1_hi = (x1_max-x1_min) / DBLE(nx_hi_res-1)
  dx2_hi = (x2_max-x2_min) / DBLE(nx_hi_res-1)
  dx3_hi = (x3_max-x3_min) / DBLE(nx_hi_res-1)
  dx4_hi = (x4_max-x4_min) / DBLE(nx_hi_res-1)
  dx5_hi = (x5_max-x5_min) / DBLE(nx_hi_res-1)

  DO ii = 1, nx_lo_res

    axis1_lo_res(ii) = x1_min + DBLE(ii-1)*dx1_lo
    axis2_lo_res(ii) = x2_min + DBLE(ii-1)*dx2_lo
    axis3_lo_res(ii) = x3_min + DBLE(ii-1)*dx3_lo
    axis4_lo_res(ii) = x4_min + DBLE(ii-1)*dx4_lo
    axis5_lo_res(ii) = x5_min + DBLE(ii-1)*dx5_lo

  END DO 

  DO ii = 1, nx_hi_res

    axis1_hi_res(ii) = x1_min + DBLE(ii-1)*dx1_hi
    axis2_hi_res(ii) = x2_min + DBLE(ii-1)*dx2_hi
    axis3_hi_res(ii) = x3_min + DBLE(ii-1)*dx3_hi
    axis4_hi_res(ii) = x4_min + DBLE(ii-1)*dx4_hi
    axis5_hi_res(ii) = x5_min + DBLE(ii-1)*dx5_hi

  END DO

  DO ii = 1, nx_lo_res
    DO jj = 1, nx_lo_res
      DO kk = 1, nx_lo_res
        DO ll = 1, nx_lo_res
          DO mm = 1, nx_lo_res

            linear_data(ii,jj,kk,ll,mm) = axis1_lo_res(ii) &
                                        + 2.0d0 * axis2_lo_res(jj) &
                                        - 1.0d0/3.0d0 * axis3_lo_res(kk) &
                                        - 7.0d0 * axis4_lo_res(ll) &
                                        + axis5_lo_res(mm)

           non_linear_data_lo_res(ii,jj,kk,ll,mm) = sin(axis1_lo_res(ii)) &
                                                  * cos(axis2_lo_res(jj)) &
                                                  * sqrt(axis3_lo_res(kk)) &
                                                  + cos(axis4_lo_res(ll) * sin(axis5_lo_res(mm)))

          END DO
        END DO
      END DO
    END DO
  END DO 

  DO ii = 1, nx_hi_res
    DO jj = 1, nx_hi_res
      DO kk = 1, nx_hi_res
        DO ll = 1, nx_hi_res
          DO mm = 1, nx_hi_res

           non_linear_data_hi_res(ii,jj,kk,ll,mm) = sin(axis1_hi_res(ii)) &
                                                  * cos(axis2_hi_res(jj)) &
                                                  * sqrt(axis3_hi_res(kk)) &
                                                  + cos(axis4_hi_res(ll) * sin(axis5_hi_res(mm)))

          END DO
        END DO
      END DO
    END DO
  END DO 

  CALL GetIndexAndDelta( x1, axis1_lo_res, ix1, dx1 )
  CALL GetIndexAndDelta( x2, axis2_lo_res, ix2, dx2 )
  CALL GetIndexAndDelta( x3, axis3_lo_res, ix3, dx3 )
  CALL GetIndexAndDelta( x4, axis4_lo_res, ix4, dx4 )
  CALL GetIndexAndDelta( x5, axis5_lo_res, ix5, dx5 )

  write(stdout,*) 'Testing interpolation with linear source data'
  write(stdout,*) 'source(x1,x2,x3,x4,x5) = x1 + 2 x2 - 1/3 x3 - 7 x4 + x5'

  !test 1D 1D linear data
  interp_result_linear = LOG10(LinearInterp_Array_Point(ix1,dx1,0.0d0,linear_data(:,1,1,1,1)))
 
  linear_value = x1 + 2.0d0 * axis2_lo_res(1) - 1.0d0/3.0d0 * axis3_lo_res(1) - 7.0d0 * axis4_lo_res(1) + axis5_lo_res(1)
  write(stdout,*) 'interp linear data 1D 1D absolute error', interp_result_linear-linear_value

  !test 2D 2D linear data
  interp_result_linear = LOG10(LinearInterp_Array_Point(ix1,ix2,dx1,dx2,0.0d0,linear_data(:,:,1,1,1)))
 
  linear_value = x1 + 2.0d0 * x2 - 1.0d0/3.0d0 * axis3_lo_res(1) - 7.0d0 * axis4_lo_res(1) + axis5_lo_res(1)
  write(stdout,*) 'interp linear data 2D 2D absolute error', interp_result_linear-linear_value

  !test 3D 3D linear data
  interp_result_linear = LOG10(LinearInterp_Array_Point(ix1,ix2,ix3, &
                                                        dx1,dx2,dx3,0.0d0,linear_data(:,:,:,1,1)))
 
  linear_value = x1 + 2.0d0 * x2 - 1.0d0/3.0d0 * x3 - 7.0d0 * axis4_lo_res(1) + axis5_lo_res(1)
  write(stdout,*) 'interp linear data 3D 3D absolute error', interp_result_linear-linear_value

  !test 4D 4D linear data
  interp_result_linear = LOG10(LinearInterp_Array_Point(ix1,ix2,ix3,ix4, &
                                                        dx1,dx2,dx3,dx4,0.0d0,linear_data(:,:,:,:,1)))

  linear_value = x1 + 2.0d0 * x2 - 1.0d0/3.0d0 * x3 - 7.0d0 * x4 + axis5_lo_res(1)
  write(stdout,*) 'interp linear data 4D 4D absolute error', interp_result_linear-linear_value

  !test 5D 5D linear data
  interp_result_linear = LOG10(LinearInterp_Array_Point( ix1, ix2, ix3, ix4, ix5, &
                                                         dx1, dx2, dx3, dx4, dx5, &
                                                         0.0d0, linear_data))

  linear_value = x1 + 2.0d0 * x2 - 1.0d0/3.0d0 * x3 - 7.0d0 * x4 + x5
  write(stdout,*) 'interp linear data 5D 5D absolute error', interp_result_linear-linear_value


  write(stdout,*) 'Testing interpolation with nonlinear source data and two different source array resolutions (factor of 2)'
  write(stdout,*) 'source(x1,x2,x3,x4,x5) = sin(x1) * cos(x2) * sqrt(x3) + cos(x4 * sin(x5))'

  !test 5D 5D nonlinear data
  interp_result_non_linear_lo = LOG10(LinearInterp_Array_Point( ix1, ix2, ix3, ix4, ix5, &
                                                                dx1, dx2, dx3, dx4, dx5, &
                                                                0.0d0, non_linear_data_lo_res))


  CALL GetIndexAndDelta( x1, axis1_hi_res, ix1, dx1 )
  CALL GetIndexAndDelta( x2, axis2_hi_res, ix2, dx2 )
  CALL GetIndexAndDelta( x3, axis3_hi_res, ix3, dx3 )
  CALL GetIndexAndDelta( x4, axis4_hi_res, ix4, dx4 )
  CALL GetIndexAndDelta( x5, axis5_hi_res, ix5, dx5 )

  interp_result_non_linear_hi = LOG10(LinearInterp_Array_Point( ix1, ix2, ix3, ix4, ix5, &
                                                                dx1, dx2, dx3, dx4, dx5, &
                                                                0.0d0, non_linear_data_hi_res))

  non_linear_value = sin(x1) * cos(x2) * sqrt(x3) + cos(x4 * sin(x5)) 

  write(stdout,*) 'interp nonlinear data 5D 5D relative error lo res', &
                  (interp_result_non_linear_lo-non_linear_value)/non_linear_value

  write(stdout,*) 'interp nonlinear data 5D 5D relative error hi res', &
                  (interp_result_non_linear_hi-non_linear_value)/non_linear_value

END PROGRAM wlTestInterpolation
