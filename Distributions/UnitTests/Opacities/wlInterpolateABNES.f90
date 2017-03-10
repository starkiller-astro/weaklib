PROGRAM wlInterpolateABNES

  USE wlKindModule, ONLY: dp
  USE wlInterpolationModule, ONLY: &
    LogInterpolateSingleVariable
  USE wlOpacityTableModule
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
  USE B85, ONLY: TotalNESKernel
  USE wlExtPhysicalConstantsModule, ONLY: &
    h, cvel_inv


  IMPLICIT NONE

!--------- parameters for creating energy grid 
  INTEGER, PARAMETER     :: Inte_nPointE = 40
  REAL(dp)               :: Inte_Emin = 1.0d-1
  REAL(dp)               :: Inte_Emax = 3.0d02
  TYPE(GridType)         :: Inte_E

!-------- variables for reading opacity table
  TYPE(OpacityTableType) :: OpacityTable

!-------- variables for reading parameters data
  REAL(dp), DIMENSION(:), ALLOCATABLE :: r, Inte_rho, Inte_Ye, Inte_T, &
                                        Inte_TMeV, Inte_cmpe, Inte_eta,&
                                        database
  REAL(dp), DIMENSION(Inte_nPointE)  :: buffer1, buffer2, buffer3, buffer4, &
                                        buffer_rho, buffer_T, buffer_ye,    &
                                        Inte_Em, Inte_ES 
  REAL(dp), DIMENSION(Inte_nPointE)  :: OneM
  REAL(dp), DIMENSION(Inte_nPointE, Inte_nPointE)  :: NESK
  CHARACTER(LEN=100)                  :: Format1, Format2, Format3, Format4
  CHARACTER(LEN=30)                   :: a,b,c,d,e,f,g,ih, ic, jc
  INTEGER, DIMENSION(4)               :: LogInterp
  INTEGER                             :: i, ii, jj, datasize
  REAL(dp)                            :: Offset_NES0, Offset_NES1, Offset_cmpe,&
                                         Offset_EM, Offset_ES, kMeV, FourPi,   &  
                                         buffer

!-------- output variables ------------------------
!  REAL(dp), DIMENSION(:), ALLOCATABLE   :: Interpolant, Interpolant2


  kMev = 8.61733d-11 ![MeV K^{-1}]
  Format1 = "(5A12)"
  Format2 = "(5ES12.3)"
  Format3 = "(10A12)"
  Format4 = "(10ES12.3E3)"
  FourPi  = 12.56d0 ! 4*pi
  OneM(:) = 1.0_dp

!----------------------------------------
!   interpolated energy 
!----------------------------------------
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

! ---------- data-in ----------------
  WRITE(*,*) 'Opening Output100ms.d file .... '
  OPEN(10, FILE = "Output100ms.d", FORM = "formatted", ACTION = 'read')
  datasize = 213! Vary from different input files

!---------------------------------------
!    interpolated rho, T, Ye
!---------------------------------------
  WRITE(*,*) 'Allocating variables '
  ALLOCATE( database( datasize * 5) )
  ALLOCATE( r( datasize ) )
  ALLOCATE( Inte_rho( datasize ) )
  ALLOCATE( Inte_Ye( datasize ) )
  ALLOCATE( Inte_T( datasize ) )
  ALLOCATE( Inte_TMeV( datasize ) )
  ALLOCATE( Inte_cmpe( datasize ) )
  ALLOCATE( Inte_eta( datasize ) )

! ----------------------------------

  READ(10, Format1 ) a,b,c,d,e
  READ(10, Format2 ) database

  WRITE(*,*) 'Reading in data ... '
  DO i = 1, datasize
    r(i) = database(i*5-4)
    Inte_rho(i) = database(i*5-3)
    Inte_T(i) = database(i*5-2)
    Inte_Ye(i) = database(i*5-1)
  END DO

  CLOSE(10, STATUS = 'keep')

!---------------------------------------
!    read in the reference table
!---------------------------------------

  CALL InitializeHDF( )
  CALL ReadOpacityTableHDF( OpacityTable, "OpacityTable.h5" )
  CALL FinalizeHDF( )

  PRINT*,'OpacityTable nPointsTS ', OpacityTable % nPointsTS
  Offset_Em = OpacityTable % thermEmAb % Offsets(1)
  Offset_ES = OpacityTable % scatt_Iso % Offsets(1,1)
  Offset_NES0 = OpacityTable % scatt_NES % Offsets(1,1)
  Offset_NES1 = OpacityTable % scatt_NES % Offsets(1,2)
  Offset_cmpe = OpacityTable % EOSTable % DV % Offsets(4)

!--------------------------------------
!   do interpolation
!--------------------------------------
  
  c = (' T (K)      ')
  e = ('chem_e (MeV)')
  f = (' T (MeV)    ')
  g = ('  e_in(e)   ')
  ih = (' thermEmAb  ')
  ic= (' R_NESK^0   ')
  jc= (' R_ES^0   ')

  PRINT*, "Creating Output.d file to store output data ..."

  OPEN( 99, FILE = "NESOutput100ms.d", FORM = "formatted", ACTION = 'write')
  WRITE(99, Format3) a,b,c,d,&
                     e,f, &
                     g,ih,jc, ic 

  ASSOCIATE( TableEmAb => OpacityTable % thermEmAb % Absorptivity(1) % Values, & 
             TableESR0 => OpacityTable % scatt_Iso % Kernel(1) % Values(:,:,:,:,1), &
!             TableNESR => OpacityTable % scatt_NES % Kernel(1) % Values(:,:,:,:,1), &
             Tablecmpe => OpacityTable % EOSTable % DV % Variables(4) % Values, &
             IntE    => Inte_E  % Values )

  CALL LogInterpolateSingleVariable &
           ( Inte_rho, Inte_T, Inte_Ye, &
             OpacityTable % EOSTable % TS % States(1) % Values, &
             OpacityTable % EOSTable % TS % States(2) % Values, &
             OpacityTable % EOSTable % TS % States(3) % Values, &
             (/1,1,0/), Offset_cmpe, &
             Tablecmpe, &
             Inte_cmpe  )

  DO i = 1, datasize  

    buffer_rho(:) = Inte_rho(i)
    buffer_T(:)   = Inte_T(i)
    buffer_ye(:)  = Inte_Ye(i)

    CALL LogInterpolateSingleVariable & 
           ( IntE, buffer_rho, buffer_T, buffer_ye,             &
             OpacityTable % EnergyGrid % Values, &
             OpacityTable % EOSTable % TS % States(1) % Values, &
             OpacityTable % EOSTable % TS % States(2) % Values, &
             OpacityTable % EOSTable % TS % States(3) % Values, &
             (/1,1,1,0/), Offset_Em, TableEmAb, Inte_EM )

    CALL LogInterpolateSingleVariable &
           ( IntE, buffer_rho, buffer_T, buffer_ye,             &
             OpacityTable % EnergyGrid % Values, &
             OpacityTable % EOSTable % TS % States(1) % Values, &
             OpacityTable % EOSTable % TS % States(2) % Values, &
             OpacityTable % EOSTable % TS % States(3) % Values, &
             (/1,1,1,0/), Offset_ES, TableESR0, Inte_ES )

   
    Inte_TMeV(i) = Inte_T(i)* kMeV
    Inte_eta(i)  = Inte_cmpe(i) / Inte_TMeV(i)



    CALL TotalNESKernel( IntE, Inte_TMeV(i), Inte_cmpe(i), 30, 0, NESK )

    buffer1(:) = 0.0_dp
    buffer2(:) = IntE(:) * IntE(:) ! ep*ep

    DO jj = 1, Inte_nPointE ! e loop
 
      buffer3(:) = NESK(jj,:) ! fixed e
      buffer4(:) = buffer2(:)*buffer3(:) ! ep*ep*NES(e,ep)

!------ Trapezoidal Rule for integraling ep
      buffer = 0.0_dp
      DO ii = 1, (Inte_nPointE-1) 

      buffer = buffer + &
                      ( IntE(ii+1) - IntE(ii) ) * 0.5_dp * &
                      ( buffer4(ii+1) + buffer4(ii) )

      IF ( buffer .lt. 0.0_dp) THEN
        WRITE(*,*) 'Negative buffer ', buffer
        RETURN
      ELSE 
      END IF
        
      END DO ! ii
!------ end integral

      WRITE(99, Format4) r(i), Inte_rho(i), Inte_T(i), Inte_Ye(i), &
                         Inte_cmpe(i), Inte_TMeV(i), &
                         IntE(jj), Inte_EM(jj), FourPi*Inte_ES(jj), buffer*FourPi
    END DO ! jj
    
  END DO ! i

  END ASSOCIATE ! Table

  CLOSE( 99, STATUS = 'keep')  

  WRITE(*,*) 'File Output.d was written/rewrtited.'
  
  CALL DeAllocateOpacityTable( OpacityTable )

END PROGRAM wlInterpolateABNES
