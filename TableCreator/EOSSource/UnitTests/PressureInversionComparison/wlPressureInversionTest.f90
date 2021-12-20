PROGRAM wlPressureInversionTest

  USE wlKindModule, ONLY: dp 
  USE wlGridModule
  USE HDF5
  USE wlEquationOfStateTableModule
  USE wlInterpolationModule
  USE wlIOModuleHDF
  USE wlEOSIOModuleHDF
                  

  implicit none

  INTEGER  :: i, j, k, TestUnit1, FileUnit 
  REAL(dp), DIMENSION(:), ALLOCATABLE :: r
  REAL(dp), DIMENSION(:), ALLOCATABLE :: rho
  REAL(dp), DIMENSION(:), ALLOCATABLE :: T
  REAL(dp), DIMENSION(:), ALLOCATABLE :: Ye 
  REAL(dp), DIMENSION(:), ALLOCATABLE :: p 
  REAL(dp), DIMENSION(:), ALLOCATABLE :: Temperature
  REAL(dp), DIMENSION(:), ALLOCATABLE :: L1norm
  TYPE(EquationOfStateTableType) :: EOSTable
  REAL(dp) :: Yemin, Yemax, logrho, logTmin, logTmax, logrhomin, logrhomax, epsilon
  INTEGER, PARAMETER :: NumPoints = 112 
  CHARACTER(len=1)   :: EOSFlag     ! nuclear eos selection flag

  LOGICAL            :: fail        ! did EoS fail to converge

 4018 FORMAT (i4,3(es12.5))
  epsilon = 1.d-100

!------------------------------------------------------------------------------
!
!  Note: The resolution of the table used to generate the profile must match
!        that of the table being used to test the inversion.  
!
!------------------------------------------------------------------------------
  FileUnit = 50 
  OPEN( unit = FileUnit, FILE="FidOutput100ms.d")
  OPEN( newunit = TestUnit1, FILE="t_from_p_test_high_res.d")

  CALL InitializeHDF( )

  CALL ReadEquationOfStateTableHDF( EOSTable, "EquationOfStateTableDoubleTRes.h5" )

  ALLOCATE( r( NumPoints ), rho( NumPoints ), T( NumPoints ), Ye( NumPoints ),  &
            p( NumPoints ), Temperature( NumPoints), L1norm( NumPoints ) )


  DO i = 1, NumPoints
    READ(FileUnit, '(5(es11.3,x))') r(i), rho(i), T(i), Ye(i), p(i)
  END DO

  DO i = 1, SIZE( rho )

      CALL ComputeTempFromIntEnergy( rho(i), p(i), Ye(i),              &
           EOSTable % TS % States(1) % Values,                          &
           EOSTable % TS % States(2) % Values,                          &
           EOSTable % TS % States(3) % Values,                          &
           EOSTable % TS % LogInterp,                                   &
           EOSTable % DV % Variables(1) % Values(:,:,:),                &
           EOSTable % DV % Offsets(1), Temperature(i) )

  END DO
  
    DO i = 1, SIZE(rho)
        L1norm(i) = ABS( T(i) - Temperature(i) )/ T(i)

      WRITE (TestUnit1,'( 4(es14.7,x) )' ) L1norm(i), T(i), Temperature(i), rho(i)! ,   &

    END DO


  CLOSE(TestUnit1)
  CLOSE(FileUnit)

  CALL DeAllocateEquationOfStateTable( EOSTable )

  CALL FinalizeHDF( )

END PROGRAM wlPressureInversionTest
