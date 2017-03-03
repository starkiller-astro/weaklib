PROGRAM wlNESKernelTest

  USE wlKindModule, ONLY: dp
  USE wlInterpolationModule, ONLY: &
    LogInterpolateDifferentiateSingleVariable, &
    LogInterpolateSingleVariable
  USE wlOpacityTableModule, ONLY: &
    OpacityTableType
  USE wlIOModuleHDF, ONLY: &
    InitializeHDF, &
    FinalizeHDF
  USE wlOpacityTableIOModuleHDF, ONLY: &
    ReadOpacityTableHDF
  USE wlGridModule, ONLY: &
    GridType, &
    AllocateGrid, &
    DescribeGrid, &
    MakeLogGrid

  IMPLICIT NONE

!--------- parameters for creating energy grid 
  INTEGER, PARAMETER     :: Inte_nPointE = 5
  REAL(dp)               :: Inte_Emin = 2.0d00
  REAL(dp)               :: Inte_Emax = 2.0d02
  TYPE(GridType)         :: Inte_E

!-------- variables for reading opacity table
  TYPE(OpacityTableType) :: OpacityTable

!-------- variables for reading parameters data
  REAL(dp), DIMENSION(:), ALLOCATABLE :: r, Inte_rho, Inte_Ye, Inte_T, Inte_TMeV, Inte_cmpe, Inte_eta, database
  REAL(dp), DIMENSION(Inte_nPointE)   :: buffer1, buffer2, buffer3
  CHARACTER(LEN=100)                  :: Format1, Format2, Format3, Format4
  CHARACTER(LEN=30)                   :: a,b,c,d,e,f,g,h, ic, jc
  INTEGER, DIMENSION(4)               :: LogInterp
  INTEGER                             :: i, ii, jj, datasize
  REAL(dp)                            :: Offset_NES0, Offset_NES1, Offset_cmpe, kMeV, Constant

!-------- output variables ------------------------
  REAL(dp), DIMENSION(:), ALLOCATABLE   :: Interpolant, Interpolant2

!----------------------------------------
!   interpolated energy 
!----------------------------------------
  kMev = 8.61733d-11 ![MeV K^{-1}]
  Format1 = "(5A12)"
  Format2 = "(5ES12.3)"
  Format3 = "(10A12)"
  Format4 = "(10ES12.3E3)"
  Constant = 6.28d0 ! 2*pi
 
  OPEN(1, FILE = "Output100ms.d", FORM = "formatted", ACTION = 'read')
  datasize = 213! Vary from different input files

  CALL AllocateGrid( Inte_E, Inte_nPointE )

  Inte_E % Unit = 'MeV                  '
  Inte_E % Name = 'Intepolated Energy   '
  Inte_E % MinValue = Inte_Emin
  Inte_E % MaxValue = Inte_Emax
  Inte_E % LogInterp = 1
  Inte_E % nPoints = Inte_nPointE
  LogInterp(1:4) = (/1, 1, 1, 1/) 

  CALL MakeLogGrid &
          ( Inte_E % MinValue, Inte_E % MaxValue, &
            Inte_E % nPoints, Inte_E % Values )

  CALL DescribeGrid( Inte_E ) 

!---------------------------------------
!    interpolated rho, T, Ye
!---------------------------------------

  ALLOCATE( database( datasize * 5) )
  ALLOCATE( r( datasize ) )
  ALLOCATE( Inte_rho( datasize ) )
  ALLOCATE( Inte_Ye( datasize ) )
  ALLOCATE( Inte_T( datasize ) )
  ALLOCATE( Inte_TMeV( datasize ) )
  ALLOCATE( Inte_cmpe( datasize ) )
  ALLOCATE( Inte_eta( datasize ) )
  ALLOCATE( Interpolant( Inte_nPointE ) )
  ALLOCATE( Interpolant2( Inte_nPointE ) )

  READ( 1, Format1 ) a,b,c,d,e
  READ( 1, Format2 ) database

  CLOSE( 1, STATUS = 'keep')  

  DO i = 1, datasize
    r(i) = database(i*5-4)
    Inte_rho(i) = database(i*5-3)
    Inte_T(i) = database(i*5-2)
    Inte_Ye(i) = database(i*5-1)
  END DO

!---------------------------------------
!    read in the reference table
!---------------------------------------
  CALL InitializeHDF( )
  CALL ReadOpacityTableHDF( OpacityTable, "OpacityTable.h5" )
  CALL FinalizeHDF( )

  PRINT*,'OpacityTable nPointsTS ', OpacityTable % nPointsTS
  Offset_NES0 = OpacityTable % scatt_NES % Offsets(1,1)
  Offset_NES1 = OpacityTable % scatt_NES % Offsets(1,2)
  Offset_cmpe = OpacityTable % EOSTable % DV % Offsets(4)
!--------------------------------------
!   do interpolation
!--------------------------------------
  
  c = (' T (K)      ')
  e = ('chem_e (MeV)')
  f = (' T (MeV)    ')
  g = ('  e_out(ep) ')
  h = ('  e_in(e)   ')
  ic= (' R_NESK^0   ')
  jc= (' R_NESK^1   ')

  PRINT*, "Creating Output.d file to store output data ..."

  OPEN( 10, FILE = "NESOutput100ms.d", FORM = "formatted", ACTION = 'write')
  WRITE(10, Format3) a,b,c,d,e,f,g,h,ic,jc

  ASSOCIATE( Tablecmpe => OpacityTable % EOSTable % DV % Variables(4) % Values, &
             Table1 => OpacityTable % scatt_NES % Kernel(1) % Values(:,:,:,:,1), &
             Table2 => OpacityTable % scatt_NES % Kernel(1) % Values(:,:,:,:,2), &
             Energy => Inte_E % Values )

  CALL LogInterpolateSingleVariable &
           ( Inte_rho, Inte_T, Inte_Ye, &
             OpacityTable % EOSTable % TS % States(1) % Values, &
             OpacityTable % EOSTable % TS % States(2) % Values, &
             OpacityTable % EOSTable % TS % States(3) % Values, &
             (/1,1,0/), Offset_cmpe, &
             Tablecmpe, &
             Inte_cmpe  )

  DO i = 1, datasize

    Inte_TMeV(i) = Inte_T(i)* kMeV
    buffer1(:)   = Inte_T(i)
    Inte_eta(i)  = Inte_cmpe(i) / Inte_TMeV(i)
    buffer2(:)   = Inte_eta(i) 

    DO ii = 1, Inte_nPointE
      buffer3(:) = Energy(ii) ! ep

      CALL LogInterpolateSingleVariable & 
           ( Energy, buffer3, buffer1, buffer2, & 
             OpacityTable % EnergyGrid % Values, &
             OpacityTable % EnergyGrid % Values, &
             OpacityTable % EOSTable % TS % States(2) % Values, &
             OpacityTable % EtaGrid % Values,    &
             LogInterp, Offset_NES0, Table1, Interpolant )

      DO jj = 1, Inte_nPointE
         IF ( Interpolant(jj) .lt. 0.0d0 ) THEN
         PRINT*, 'ERROR: NEGATIVE R_0 ! '
         RETURN
         END IF
      END DO   
 
      CALL LogInterpolateSingleVariable &
           ( Energy, buffer3, buffer1, buffer2, &
             OpacityTable % EnergyGrid % Values, &
             OpacityTable % EnergyGrid % Values, &
             OpacityTable % EOSTable % TS % States(2) % Values, &
             OpacityTable % EtaGrid % Values,    &
             LogInterp, Offset_NES1, Table2, Interpolant2 )

      DO jj = 1, Inte_nPointE
        WRITE(10, Format4) r(i), Inte_rho(i), Inte_T(i), Inte_Ye(i), &
                           Inte_cmpe(i), Inte_TMeV(i), &
                           buffer3(ii), Energy(jj),  &              
                           Interpolant(jj)*Constant, &
                           Interpolant2(jj)*Constant
      END DO ! jj
    END DO ! ii

  END DO ! i

  END ASSOCIATE ! Table

  CLOSE( 10, STATUS = 'keep')  

  WRITE(*,*) 'File Output.d was written/rewrtited.'

END PROGRAM wlNESKernelTest
