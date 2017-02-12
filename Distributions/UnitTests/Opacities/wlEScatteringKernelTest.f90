PROGRAM wlEScatteringKernelTest

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

!--------- parameters for creating energy grid 
  INTEGER, PARAMETER     :: Inte_nPointE = 5
  REAL(dp)               :: Inte_Emin = 2.0d00
  REAL(dp)               :: Inte_Emax = 2.0d02
  TYPE(EnergyGridType)   :: Inte_E

!-------- variables for reading opacity table
  TYPE(OpacityTableType) :: OpacityTable

!-------- variables for reading parameters data
  REAL(dp), DIMENSION(:), ALLOCATABLE :: r, Inte_rho, Inte_T, Inte_Ye,&
                                         database
  REAL(dp), DIMENSION(Inte_nPointE)   :: buffer1, buffer2, buffer3, buffer_r
  CHARACTER(LEN=100)                  :: Format1, Format2, Format3, Format4
  CHARACTER(LEN=30)                   :: a,b,c,d,e,f,g,h,l,a1,a2,a3,a4
  INTEGER, DIMENSION(4)               :: LogInterp
  INTEGER                             :: i, ii, datasize
  REAL(dp)                            :: Offset_Em, Offset_ES
  REAL(dp)                            :: fourPi= 4.0*3.1415926
!-------- output variables ------------------------
  REAL(dp), DIMENSION(:), ALLOCATABLE :: Inte_O, Inte_R0, Inte_R1,&
                                         GONa, GONb, GOEa, GOEb

!----------------------------------------
!   interpolated energy 
!----------------------------------------
 
  Format1 = "(5A12)"
  Format2 = "(5ES12.3)"
  Format3 = "(13A16)"
  Format4 = "(13ES16.8)"

  OPEN(1, FILE = "Output100ms.d", FORM = "formatted", ACTION = 'read')
  datasize = 213

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

  ALLOCATE( database( datasize*5 ) )
  ALLOCATE( r( datasize ) )
  ALLOCATE( Inte_rho( datasize ) )
  ALLOCATE( Inte_T( datasize ) )
  ALLOCATE( Inte_Ye( datasize ) )
  ALLOCATE( Inte_O( Inte_nPointE ) )
  ALLOCATE( Inte_R0( Inte_nPointE ) )
  ALLOCATE( Inte_R1( Inte_nPointE ) )
  ALLOCATE( GONa( Inte_nPointE ) )
  ALLOCATE( GONb( Inte_nPointE ) )
  ALLOCATE( GOEa( Inte_nPointE ) )
  ALLOCATE( GOEb( Inte_nPointE ) )

  READ( 1, Format1 ) a,b,c,d,e
  READ( 1, Format2 ) database
 
  DO i = 1, datasize
    r(i) = database(i*5-4)
    Inte_rho(i) = database(i*5-3)
    Inte_T(i) = database(i*5-2)
    Inte_Ye(i) = database(i*5-1)
  END DO

  CLOSE( 1, STATUS = 'keep')  
 
  LogInterp(2:4) = (/1, 1, 0/)     ! rho and T is LogGrid, Ye is linear
 
!---------------------------------------
!    read in the reference table
!---------------------------------------
  CALL InitializeHDF( )
  CALL ReadOpacityTableHDF( OpacityTable, "OpacityTable.h5" )
  CALL FinalizeHDF( )

  Offset_Em = OpacityTable % thermEmAb % Offset
  Offset_ES = OpacityTable % scatt_Iso % Offset
!--------------------------------------
!   do interpolation
!--------------------------------------
  e = ('  energy  ')
  f = ('  Opacity ')
  g = ('    R0    ')
  h = ('    R1    ')  
  l = ('  ratio   ')
  a1= ('GON typeA ')
  a2= ('GON typeB ')
  a3= ('GOE typeA ')
  a4= ('GOE typeB ')

  OPEN( 10, FILE = "ESOutput100ms_5quad.d", FORM = "formatted", ACTION = 'write')
  WRITE(10, Format3) a,b,c,d,e,f,g,h,l,a1,a2,a3,a4

  ASSOCIATE( Table1 => OpacityTable % thermEmAb % Absorptivity(1) % Values,&
             Table4 => OpacityTable % thermEmAb % GreyOpacity_Number_FD(1)% Values,&
             Table5 => OpacityTable % thermEmAb % GreyOpacity_Energy_FD(1)% Values,&
             Table2 => OpacityTable % scatt_Iso % Kernel(1) % Values(:,:,:,:,1),&
             Table3 => OpacityTable % scatt_Iso % Kernel(1) % Values(:,:,:,:,2),&
             Table6 => OpacityTable % scatt_Iso % GreyOpacity_Number_FD(1)% Values(:,:,:,1),&
             Table7 => OpacityTable % scatt_Iso % GreyOpacity_Energy_FD(1)% Values(:,:,:,1),&
             Energy => Inte_E % Values )

  DO i = 1, datasize

    buffer1(:) = Inte_rho(i)
    buffer2(:) = Inte_T(i)
    buffer3(:) = Inte_Ye(i)

    CALL LogInterpolateSingleVariable & 
           ( Energy, buffer1, buffer2, buffer3, &
             OpacityTable % EnergyGrid % Values, &
             OpacityTable % EOSTable % TS % States(1) % Values, &
             OpacityTable % EOSTable % TS % States(2) % Values, &
             OpacityTable % EOSTable % TS % States(3) % Values, &
             LogInterp, Offset_Em, Table1, Inte_O )

    CALL LogInterpolateSingleVariable &
           ( buffer1, buffer2, buffer3, &
             OpacityTable % EOSTable % TS % States(1) % Values, &
             OpacityTable % EOSTable % TS % States(2) % Values, &
             OpacityTable % EOSTable % TS % States(3) % Values, &
             LogInterp, Offset_Em, Table4, GONa )

    CALL LogInterpolateSingleVariable &
           ( buffer1, buffer2, buffer3, &
             OpacityTable % EOSTable % TS % States(1) % Values, &
             OpacityTable % EOSTable % TS % States(2) % Values, &
             OpacityTable % EOSTable % TS % States(3) % Values, &
             LogInterp, Offset_Em, Table5, GOEa )

    CALL LogInterpolateSingleVariable &
           ( Energy, buffer1, buffer2, buffer3, &
             OpacityTable % EnergyGrid % Values, &
             OpacityTable % EOSTable % TS % States(1) % Values, &
             OpacityTable % EOSTable % TS % States(2) % Values, &
             OpacityTable % EOSTable % TS % States(3) % Values, &
             LogInterp, Offset_ES, Table2, Inte_R0 )

    CALL LogInterpolateSingleVariable &
           ( buffer1, buffer2, buffer3, &
             OpacityTable % EOSTable % TS % States(1) % Values, &
             OpacityTable % EOSTable % TS % States(2) % Values, &
             OpacityTable % EOSTable % TS % States(3) % Values, &
             LogInterp, Offset_ES, Table6, GONb )

    CALL LogInterpolateSingleVariable &
           ( buffer1, buffer2, buffer3, &
             OpacityTable % EOSTable % TS % States(1) % Values, &
             OpacityTable % EOSTable % TS % States(2) % Values, &
             OpacityTable % EOSTable % TS % States(3) % Values, &
             LogInterp, Offset_ES, Table7, GOEb )

    CALL LogInterpolateSingleVariable &
           ( Energy, buffer1, buffer2, buffer3, &
             OpacityTable % EnergyGrid % Values, &
             OpacityTable % EOSTable % TS % States(1) % Values, &
             OpacityTable % EOSTable % TS % States(2) % Values, &
             OpacityTable % EOSTable % TS % States(3) % Values, &
             LogInterp, Offset_ES, Table3, Inte_R1 )

    DO ii = 1, Inte_nPointE

      IF ( Inte_R0(ii)*fourPi == Inte_O(ii) ) THEN
        buffer_r = 1.0
      ELSE IF( Inte_R0(ii)*fourPi .gt. Inte_O(ii)*10.0**50  ) THEN
        buffer_r = 10.0**(-50)
      ELSE 
        buffer_r = Inte_R0(ii)*fourPi/Inte_O(ii)
      END IF

      WRITE(10, Format4) r(i), buffer1(ii), buffer2(ii), buffer3(ii), &
                         Inte_E % Values(ii), Inte_O(ii), Inte_R0(ii),&
                         Inte_R1(ii), buffer_r,&
                         GONa(ii), GONb(ii), GOEa(ii), GOEb(ii)
    END DO ! ii

  END DO ! i

  END ASSOCIATE 

  CLOSE( 10, STATUS = 'keep')  

  WRITE(*,*) 'File ESOutput*ms.d was written/rewrtited.'

END PROGRAM wlEScatteringKernelTest
