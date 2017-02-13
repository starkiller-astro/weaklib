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
  INTEGER, PARAMETER     :: Inte_nPointE = 30
  REAL(dp)               :: Inte_Emin = 2.0d00
  REAL(dp)               :: Inte_Emax = 2.0d02
  TYPE(GridType)         :: Inte_E

!-------- variables for reading opacity table
  TYPE(OpacityTableType) :: OpacityTable

!-------- variables for reading parameters data
  REAL(dp), DIMENSION(:), ALLOCATABLE :: Inte_T, Inte_eta, database
  REAL(dp), DIMENSION(Inte_nPointE)   :: buffer1, buffer2, buffer3
  CHARACTER(LEN=100)                  :: Format1, Format2, Format3, Format4
  CHARACTER(LEN=30)                   :: a,b,c,d,e,f,g,h
  INTEGER, DIMENSION(4)               :: LogInterp
  INTEGER                             :: i, ii, jj, datasize
  REAL(dp)                            :: Offset_NES, kMeV

!-------- output variables ------------------------
  REAL(dp), DIMENSION(:), ALLOCATABLE   :: Interpolant

!----------------------------------------
!   interpolated energy 
!----------------------------------------
  kMev = 8.61733d-11 ![MeV K^{-1}]
  Format1 = "(2A12)"
  Format2 = "(2ES12.3)"
  Format3 = "(5A12)"
  Format4 = "(5ES12.3)"

  OPEN(1, FILE = "Inputfile_stand_NESK.d", FORM = "formatted", ACTION = 'read')
  datasize = 1  ! Vary from different input files

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

  ALLOCATE( database( datasize * 2) )
  ALLOCATE( Inte_T( datasize ) )
  ALLOCATE( Inte_eta( datasize ) )
  ALLOCATE( Interpolant( Inte_nPointE ) )

  READ( 1, Format1 ) a,b
  READ( 1, Format2 ) database

  CLOSE( 1, STATUS = 'keep')  

  DO i = 1, datasize  
    Inte_T(i) = database(i*2-1)
    Inte_eta(i) = database(i*2)
  END DO 

!---------------------------------------
!    read in the reference table
!---------------------------------------
  CALL InitializeHDF( )
  CALL ReadOpacityTableHDF( OpacityTable, "OpacityTable.h5" )
  CALL FinalizeHDF( )

  Offset_NES = OpacityTable % scatt_NES % Offset

!--------------------------------------
!   do interpolation
!--------------------------------------
  e = ('  e_out(ep) ')
  f = ('  e_in(e)   ')
  g = (' R_NESK^0   ')

  PRINT*, "Creating Output.d file to store output data ..."
  OPEN( 10, FILE = "Output.d", FORM = "formatted", ACTION = 'write')
  WRITE(10, Format3) a,b,e,f,g

  ASSOCIATE( Table  => OpacityTable % scatt_NES % Kernel(1) % Values(:,:,:,:,1), &
             Energy => Inte_E % Values )

  DO i = 1, datasize

    buffer1(:) = Inte_T(i) / kMeV
    buffer2(:) = Inte_eta(i) 

    DO ii = 1, Inte_nPointE
      buffer3(:) = Energy(ii) ! ep

      CALL LogInterpolateSingleVariable & 
           ( Energy, buffer3, buffer1, buffer2, & 
             OpacityTable % EnergyGrid % Values, &
             OpacityTable % EnergyGrid % Values, &
             OpacityTable % EOSTable % TS % States(2) % Values, &
             OpacityTable % EtaGrid % Values,    &
             LogInterp, Offset_NES, Table, Interpolant )
  
      DO jj = 1, Inte_nPointE
        WRITE(10, Format4) buffer1(ii), buffer2(ii), &
                           buffer3(ii), Energy(jj),  &              
                           1.0!!Interpolant(jj)
      END DO ! jj
    END DO ! ii

  END DO ! i

  END ASSOCIATE ! Table

  CLOSE( 10, STATUS = 'keep')  

  WRITE(*,*) 'File Output.d was written/rewrtited.'

END PROGRAM wlNESKernelTest
