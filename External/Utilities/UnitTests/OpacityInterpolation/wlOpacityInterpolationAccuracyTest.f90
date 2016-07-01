PROGRAM wlOpacityInterpolationAccuracyTest

  USE wlKindModule, ONLY: dp
  USE wlInterpolationModule
  USE wlOpacityTableModule 
  USE wlIOModuleHDF, ONLY:&
          InitializeHDF, FinalizeHDF
  USE wlOpacityTableIOModuleHDF
  USE wlEnergyGridModule
  USE wlGridModule, ONLY: MakeLogGrid
  USE wlExtNumericalModule, ONLY: epsilon
  USE B85, ONLY: totalECapEm
  IMPLICIT NONE

  INTEGER  :: TestUnit1, TestUnit2
!--------- parameters for creating energy grid 
  INTEGER, PARAMETER     :: Inte_nPointE = 240
  REAL(dp)               :: Inte_Emin = 2.0d00
  REAL(dp)               :: Inte_Emax = 2.0d02
  TYPE(EnergyGridType)   :: Inte_E

!-------- variables for reading opacity table
  TYPE(OpacityTableType) :: OpacityTable

!-------- variables for reading parameters data
  REAL(dp), DIMENSION(:), ALLOCATABLE :: r, Inte_rho, Inte_T, Inte_Ye, e_int, &
                                         Inte_Z, Inte_A, Inte_chem_e, &
                                         Inte_chem_n, Inte_chem_p, &
                                         Inte_xheavy, Inte_xn, Inte_xp
  REAL(dp)                            :: tempxa
  REAL(dp), DIMENSION(Inte_nPointE)   :: buffer1, buffer2, buffer3
  CHARACTER(LEN=100)                  :: Format1, Format2, Format3, Format4
  CHARACTER(LEN=30)                   :: name1, name2, name3, name4, name5, &
                                         name6, name7, name8, name9, name10,&
                                         name11, name12

  INTEGER, DIMENSION(4)               :: LogInterp
  INTEGER                             :: i, ii, datasize
  REAL(dp)                            :: Offset

!-------- output variables ------------------------
  REAL(dp), DIMENSION(:), ALLOCATABLE   :: Interpolant
  REAL(dp), DIMENSION(:), ALLOCATABLE   :: DirectCall
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: Derivative

!----------------------------------------
!   interpolated energy 
!----------------------------------------
 
  Format1 = "(14A12)"
  Format2 = "(14ES12.3)"
  Format3 = "(12A12)"
  Format4 = "(12ES12.3)"

!  OPEN(1, FILE = "FidOutput0ms.d", FORM = "formatted", ACTION = 'read')
!  datasize = 312

  OPEN(1, FILE = "FidOutput100ms.d", FORM = "formatted", ACTION = 'read')
  datasize = 212

!  OPEN(1, FILE = "Output0ms.d", FORM = "formatted", ACTION = 'read')
!  datasize = 292

! OPEN(1, FILE = "Output100ms.d", FORM = "formatted", ACTION = 'read')
! datasize = 217
  
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

  Offset = epsilon
  ALLOCATE( r( datasize ) )
  ALLOCATE( e_int( datasize ) )
  ALLOCATE( Inte_rho( datasize ) )
  ALLOCATE( Inte_T( datasize ) )
  ALLOCATE( Inte_Ye( datasize ) )
  ALLOCATE( Inte_Z( datasize ) )
  ALLOCATE( Inte_A( datasize ) )
  ALLOCATE( Inte_chem_e( datasize ) )
  ALLOCATE( Inte_chem_n( datasize ) )
  ALLOCATE( Inte_chem_p( datasize ) )
  ALLOCATE( Inte_xheavy( datasize ) )
  ALLOCATE( Inte_xn( datasize ) )
  ALLOCATE( Inte_xp( datasize ) )
  ALLOCATE( Interpolant( Inte_nPointE ) )
  ALLOCATE( DirectCall( Inte_nPointE ) )
  ALLOCATE( Derivative( Inte_nPointE, 4 ) )

  READ( 1, Format1 ) name1, name2, name3, name4, name5, name6, name7, name8, name9, name10
  WRITE( *, Format1 ) name1, name2, name3, name4, name5, name6, name7, name8, name9, name10

  DO i = 1, datasize
    READ( 1, Format2 ) r(i), Inte_rho(i), Inte_T(i), Inte_Ye(i), e_int(i), &
                      Inte_chem_n(i), Inte_chem_p(i), Inte_chem_e(i), &
                      Inte_xn(i), Inte_xp(i), tempxa, Inte_xheavy(i), &
                      Inte_Z(i), Inte_A(i)
  END DO

  CLOSE( 1, STATUS = 'keep')  

  LogInterp(2:4) = (/1, 1, 0/)     ! rho and T is LogGrid, Ye is linear
 
!---------------------------------------
!    read in the reference table
!---------------------------------------
  CALL InitializeHDF( )
  CALL ReadOpacityTableHDF( OpacityTable, "OpacityTable.h5" )
  CALL FinalizeHDF( )

!--------------------------------------
!   do interpolation
!--------------------------------------
  name5  = ('    energy  ')
  name6  = (' direct call')
  name7  = ('   Opacity  ')
  name8  = ('    deriv T ')
  name9  = ('   deriv Ye ')
  name10 = ('error prese ')
  name11 = ('   deriv E  ')
  name12 = ('   deriv rho')
!  OPEN( 10, FILE = "IntOutput0ms.d", FORM = "formatted", ACTION = 'write')
  OPEN( 10, FILE = "IntOutput100ms.d", FORM = "formatted", ACTION = 'write')
  WRITE(10, Format3) name1, name2, name3, name4, name5, name6, name7, &
                     name10, name8, name9, name11, name12

  ASSOCIATE( Table => OpacityTable % thermEmAb % Absorptivity(1) % Values, &
             Energy => Inte_E % Values )

  DO i = 1, datasize

    buffer1(:) = Inte_rho(i)
    buffer2(:) = Inte_T(i)
    buffer3(:) = Inte_Ye(i)

    CALL LogInterpolateDifferentiateSingleVariable & 
           ( Energy, buffer1, buffer2, buffer3, & 
             OpacityTable % EnergyGrid % Values, &
             OpacityTable % EOSTable % TS % States(1) % Values, &
             OpacityTable % EOSTable % TS % States(2) % Values, &
             OpacityTable % EOSTable % TS % States(3) % Values, &
             LogInterp, Offset, Table, Interpolant, Derivative, .FALSE. )
  
    DO ii = 1, Inte_nPointE

      DirectCall( ii ) = totalEcapEm( Energy(ii), Inte_rho(i), Inte_T(i), &
                               Inte_Z(i), Inte_A(i), Inte_chem_e(i), &
                               Inte_chem_n(i), Inte_chem_p(i), Inte_xheavy(i), &
                               Inte_xn(i), Inte_xp(i) )
 
      WRITE(10, Format4) r(i), buffer1(ii), buffer2(ii), buffer3(ii), &
                         Energy(ii), DirectCall(ii), Interpolant(ii),&
                      ABS( Interpolant(ii)-DirectCall(ii) ) / DirectCall(ii),&
                         Derivative(ii,3), Derivative(ii,4), &
                         Derivative(ii,1), Derivative(ii,2)
    END DO ! ii

  END DO ! i

  END ASSOCIATE ! Table

  CLOSE( 10, STATUS = 'keep')  

END PROGRAM wlOpacityInterpolationAccuracyTest
