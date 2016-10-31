PROGRAM wlInterpolation

  USE wlKindModule, ONLY: dp
  USE wlInterpolationModule, ONLY: &
    LogInterpolateDifferentiateSingleVariable
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

!--------- parameters for creating energy grid 
  INTEGER, PARAMETER     :: Inte_nPointE = 30
  REAL(dp)               :: Inte_Emin = 2.0d00
  REAL(dp)               :: Inte_Emax = 2.0d02
  TYPE(EnergyGridType)   :: Inte_E

!-------- variables for reading opacity table
  TYPE(OpacityTableType) :: OpacityTable

!-------- variables for reading parameters data
  REAL(dp), DIMENSION(:), ALLOCATABLE :: Inte_rho, Inte_T, Inte_Ye, database
  REAL(dp), DIMENSION(Inte_nPointE)   :: buffer1, buffer2, buffer3
  CHARACTER(LEN=100)                  :: Format1, Format2, Format3, Format4
  CHARACTER(LEN=30)                   :: a,b,c,d,e,f,g,h
  INTEGER, DIMENSION(4)               :: LogInterp
  INTEGER                             :: i, ii, datasize
  REAL(dp)                            :: Offset

!-------- output variables ------------------------
  REAL(dp), DIMENSION(:), ALLOCATABLE   :: Interpolant
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: Derivative

!----------------------------------------
!   interpolated energy 
!----------------------------------------
 
  Format1 = "(3A12)"
  Format2 = "(3ES12.3)"
  Format3 = "(7A12)"
  Format4 = "(7ES12.3)"

  OPEN(1, FILE = "Inputfile_stand.d", FORM = "formatted", ACTION = 'read')
  datasize = 1

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

  ALLOCATE( database( datasize * 3) )
  ALLOCATE( Inte_rho( datasize ) )
  ALLOCATE( Inte_T( datasize ) )
  ALLOCATE( Inte_Ye( datasize ) )
  ALLOCATE( Interpolant( Inte_nPointE ) )
  ALLOCATE( Derivative( Inte_nPointE, 4 ) )

  READ( 1, Format1 ) a,b,c
  READ( 1, Format2 ) database

  CLOSE( 1, STATUS = 'keep')  

  DO i = 1, datasize  
    Inte_rho(i) = database(i*3-2)
    Inte_T(i) = database(i*3-1)
    Inte_Ye(i) = database(i*3)
  END DO 

  LogInterp(2:4) = (/1, 1, 0/)     ! rho and T is LogGrid, Ye is linear
 
!---------------------------------------
!    read in the reference table
!---------------------------------------
  CALL InitializeHDF( )
  CALL ReadOpacityTableHDF( OpacityTable, "OpacityTable.h5" )
  CALL FinalizeHDF( )

  Offset = OpacityTable % thermEmAb % Offset
!--------------------------------------
!   do interpolation
!--------------------------------------
  e = ('    energy  ')
  f = ('    Inter  ')
  g = ('   deriv T ')
  h = ('   deriv Ye')  

  OPEN( 10, FILE = "Output.d", FORM = "formatted", ACTION = 'write')
  WRITE(10, Format3) a,b,c,e,f,g,h

  ASSOCIATE( Table  => OpacityTable % thermEmAb % Absorptivity(1) % Values,&
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
      WRITE(10, Format4) buffer1(ii), buffer2(ii), buffer3(ii), &
                         Inte_E % Values(ii), Interpolant(ii), Derivative(ii,3),&
                         Derivative(ii,4)
    END DO ! ii

  END DO ! i

  END ASSOCIATE ! Table

  CLOSE( 10, STATUS = 'keep')  

  WRITE(*,*) 'File Output.d was written/rewrtited.'

END PROGRAM wlInterpolation
