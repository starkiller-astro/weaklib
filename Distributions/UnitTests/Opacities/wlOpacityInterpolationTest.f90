PROGRAM wlOpacityInterpolationTest

  USE wlKindModule, ONLY: dp
  USE wlInterpolationModule
  USE wlOpacityTableModule 
  USE wlIOModuleHDF, ONLY:&
          InitializeHDF, FinalizeHDF
  USE wlOpacityTableIOModuleHDF
  USE wlEnergyGridModule
  USE wlGridModule, ONLY: MakeLogGrid
  IMPLICIT NONE

  TYPE(OpacityTableType) :: OpacityTable
  REAL(dp), DIMENSION(:), ALLOCATABLE :: Inte_rho, Inte_T, Inte_Ye
  INTEGER, DIMENSION(4)  :: LogInterp
  REAL(dp)               :: Offset = 0.0_dp  ! Will Be Modified
  REAL(dp), DIMENSION(:), ALLOCATABLE :: Interpolant
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: Derivative
!----------------------------------------
!   interpolated energy 
!----------------------------------------
  INTEGER                :: Inte_nPointE = 3
  REAL(dp)               :: Inte_Emin = 2.0d00
  REAL(dp)               :: Inte_Emax = 2.0d02
  TYPE(EnergyGridType) :: Inte_E  

  CALL AllocateEnergyGrid( Inte_E, Inte_nPointE )

  Inte_E % Unit = 'MeV                  '
  Inte_E % Name = 'Intepolated Energy   '
  Inte_E % MinValue = Inte_Emin
  Inte_E % MaxValue = Inte_Emax
  Inte_E % LogInterp = 1
  Inte_E % nPoints = Inte_nPointE
  LogInterp(1) = 1                 ! EnergyGrid is LogGrid

  CALL MakeLogGrid &
          ( Inte_E % MinValue, Inte_E % MaxValue, &
            Inte_E % nPoints, Inte_E % Values )

  CALL DescribeEnergyGrid( Inte_E ) 

!---------------------------------------
!    interpolated rho, T, Ye
!---------------------------------------
  ALLOCATE( Inte_rho(3) )
  ALLOCATE( Inte_T(3) )
  ALLOCATE( Inte_Ye(3) )
  ALLOCATE( Interpolant(3) )
  ALLOCATE( Derivative(3,4) )
  LogInterp(2:4) = (/1, 1, 0/)     ! rho and T is LogGrid, Ye is linear
  Inte_rho = (/1.019d08, 1.072d08, 1.126d08 /)
  Inte_T   = (/6.321d09, 6.265d09, 6.411d09 /)
  Inte_Ye  = (/4.959d-1, 4.959d-1, 4.959d-1 /)
 
  WRITE(*,*), 'The interpolated points is rho:', Inte_rho
  WRITE(*,*), '                            T :', Inte_T
  WRITE(*,*), '                           Ye :', Inte_Ye
 
!---------------------------------------
!    read in the reference table
!---------------------------------------
  CALL InitializeHDF( )
  CALL ReadOpacityTableHDF( OpacityTable, "OpacityTable_Log.h5" )
  CALL FinalizeHDF( )

!  CALL DescribeOpacityTable( OpacityTable )

!--------------------------------------
!   do interpolation
!--------------------------------------
  ASSOCIATE( Table => OpacityTable % thermEmAb % Absorptivity(1) % Values )

    CALL LogInterpolateDifferentiateSingleVariable & 
           (Inte_E % Values, Inte_rho, Inte_T, Inte_Ye, & 
            OpacityTable % EnergyGrid % Values, &
            OpacityTable % EOSTable % TS % States(1) % Values, &
            OpacityTable % EOSTable % TS % States(2) % Values, &
            OpacityTable % EOSTable % TS % States(3) % Values, &
            LogInterp, Offset, Table, Interpolant, Derivative, .FALSE. )

  END ASSOCIATE ! Table

  WRITE(*,*) 'The Interpolant is ', Interpolant
  WRITE(*,*) 'The Derivative is ', Derivative

END PROGRAM wlOpacityInterpolationTest
