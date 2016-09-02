PROGRAM wlGreyTableTest

  USE wlKindModule, ONLY: dp
  USE wlInterpolationModule, ONLY: &
    LogInterpolateSingleVariable
  USE wlOpacityTableModule, ONLY: &
    OpacityTableType
  USE wlIOModuleHDF, ONLY: &
    InitializeHDF, &
    FinalizeHDF
  USE wlOpacityTableIOModuleHDF, ONLY: &
    ReadOpacityTableHDF
  USE wlEnergyGridModule, ONLY: &
    EnergyGridType, &
    AllocateEnergyGrid, &
    DescribeEnergyGrid
  USE wlGridModule, ONLY: &
    MakeLogGrid

  IMPLICIT NONE

!-------- variables for reading opacity table
  TYPE(OpacityTableType) :: OpacityTable

!-------- variables for reading parameters data
  REAL(dp), DIMENSION(:), ALLOCATABLE :: r, Inte_rho, Inte_T, Inte_Ye,&
                                         database
  CHARACTER(LEN=100)                  :: Format1, Format2, Format3, Format4
  CHARACTER(LEN=30)                   :: a,b,c,d,e,f,g,h
  INTEGER, DIMENSION(3)               :: LogInterp
  INTEGER                             :: i, ii, datasize
  REAL(dp)                            :: Offset

!-------- output variables ------------------------
  REAL(dp), DIMENSION(:), ALLOCATABLE   :: GreyMN,GreyME,GreyON,GreyOE

!----------------------------------------
!   interpolated energy 
!----------------------------------------
 
  Format1 = "(5A12)"
  Format2 = "(5ES12.3)"
  Format3 = "(6A12)"
  Format4 = "(6ES12.3)"


!  OPEN(1, FILE = "CRtest.d", FORM = "formatted", ACTION = 'read')
!  datasize = 2
  OPEN(1, FILE = "Output100ms.d", FORM = "formatted", ACTION = 'read')
  datasize = 213
  
!---------------------------------------
!    interpolated rho, T, Ye
!---------------------------------------

  ALLOCATE( database( datasize*5 ) )
  ALLOCATE( r( datasize ) )
  ALLOCATE( Inte_rho( datasize ) )
  ALLOCATE( Inte_T( datasize ) )
  ALLOCATE( Inte_Ye( datasize ) )
  ALLOCATE( GreyMN( datasize ) )
  ALLOCATE( GreyME( datasize ) )
  ALLOCATE( GreyON( datasize ) )
  ALLOCATE( GreyOE( datasize ) )

  READ( 1, Format1 ) a,b,c,d,e
  READ( 1, Format2 ) database
 
  DO i = 1, datasize
    r(i) = database(i*5-4)
    Inte_rho(i) = database(i*5-3)
    Inte_T(i) = database(i*5-2)
    Inte_Ye(i) = database(i*5-1)
    WRITE(*,*) r(i), Inte_rho(i), Inte_T(i), Inte_Ye(i)
  END DO

  CLOSE( 1, STATUS = 'keep')  
 
  LogInterp(1:3) = (/1, 1, 0/)     ! rho and T is LogGrid, Ye is linear
 
!---------------------------------------
!    read in the reference table
!---------------------------------------
  CALL InitializeHDF( )
  CALL ReadOpacityTableHDF( OpacityTable, "OpTa_Grey_5quad.h5" )
  CALL FinalizeHDF( )

  Offset = OpacityTable % thermEmAb % Offset
!--------------------------------------
!   do interpolation
!--------------------------------------
  e = ('    MNO   ')
  f = ('    MEO   ')
  g = ('    GON   ')
  h = ('GreyOpa_En')  

  OPEN( 10, FILE = "GreyOutput100ms_5quad.d", FORM = "formatted", ACTION = 'write')
  WRITE(10, Format3) a,b,c,d,e,f

  ASSOCIATE( Table  => OpacityTable % thermEmAb % GreyMoment_Number_FD(1)% Values )
    CALL LogInterpolateSingleVariable & 
           ( Inte_rho, Inte_T, Inte_Ye, & 
             OpacityTable % EOSTable % TS % States(1) % Values, &
             OpacityTable % EOSTable % TS % States(2) % Values, &
             OpacityTable % EOSTable % TS % States(3) % Values, &
             LogInterp, Offset, Table, GreyMN )
  END ASSOCIATE ! Table

  ASSOCIATE( Table  => OpacityTable % thermEmAb % GreyMoment_Energy_FD(1)% Values )
    CALL LogInterpolateSingleVariable &
           ( Inte_rho, Inte_T, Inte_Ye, &
            OpacityTable % EOSTable % TS % States(1) % Values, &
            OpacityTable % EOSTable % TS % States(2) % Values, &
            OpacityTable % EOSTable % TS % States(3) % Values, &
            LogInterp, Offset, Table, GreyME )
  END ASSOCIATE ! Table

  ASSOCIATE( Table  => OpacityTable % thermEmAb % GreyOpacity_Number_FD(1)% Values )
    CALL LogInterpolateSingleVariable &
           ( Inte_rho, Inte_T, Inte_Ye, &
             OpacityTable % EOSTable % TS % States(1) % Values, &
             OpacityTable % EOSTable % TS % States(2) % Values, &
             OpacityTable % EOSTable % TS % States(3) % Values, &
             LogInterp, Offset, Table, GreyON )
  END ASSOCIATE ! Table

  ASSOCIATE( Table  => OpacityTable % thermEmAb % GreyOpacity_Energy_FD(1)% Values )
    CALL LogInterpolateSingleVariable &
           ( Inte_rho, Inte_T, Inte_Ye, &
             OpacityTable % EOSTable % TS % States(1) % Values, &
             OpacityTable % EOSTable % TS % States(2) % Values, &
             OpacityTable % EOSTable % TS % States(3) % Values, &
             LogInterp, Offset, Table, GreyOE )
  END ASSOCIATE ! Table

  DO i = 1, datasize
    WRITE(*,*)GreyMN(i), GreyME(i), GreyON(i), GreyOE(i)
    WRITE(10, Format4) r(i), Inte_rho(i), Inte_T(i), Inte_Ye(i), &
                  !     GreyMN(i), GreyME(i), GreyON(i), GreyOE(i)
                       GreyON(i)/GreyMN(i), GreyOE(i)/GreyME(i)
  END DO ! ii

  CLOSE( 10, STATUS = 'keep')  

  WRITE(*,*) 'File GreyOutput*ms.d was written/rewrtited.'

END PROGRAM wlGreyTableTest
