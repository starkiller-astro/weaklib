PROGRAM wlInterpolateNES

  USE wlKindModule, ONLY: dp
  USE wlInterpolationModule, ONLY: &
    LogInterpolateSingleVariable, &
    LogInterpolateSingleVariable_2D2D_Custom
  USE wlOpacityTableModule, ONLY: &
    OpacityTableType, &
    DeAllocateOpacityTable
  USE wlIOModuleHDF, ONLY: &
    InitializeHDF, &
    FinalizeHDF, &
    OpenFileHDF, &
    CloseFileHDF, &
    WriteHDF, &
    OpenGroupHDF, &
    CloseGroupHDF
  USE wlOpacityTableIOModuleHDF, ONLY: &
    ReadOpacityTableHDF
  USE wlGridModule, ONLY: &
    GridType, &
    AllocateGrid, &
    DescribeGrid, &
    MakeLogGrid
  USE wlExtPhysicalConstantsModule, ONLY: kMeV, ca, cv
  USE wlExtNumericalModule, ONLY: pi, half, twpi, zero
  USE HDF5

  IMPLICIT NONE

!--------- parameters for creating energy grid 
  INTEGER, PARAMETER     :: Inte_nPointE = 40
  REAL(dp)               :: Inte_Emin = 1.0d-1
  REAL(dp)               :: Inte_Emax = 3.0d02
  TYPE(GridType)         :: Inte_E

!-------- variables for reading opacity table
  TYPE(OpacityTableType) :: OpacityTable

!-------- variables for reading parameters data
  REAL(dp), DIMENSION(:), ALLOCATABLE :: Inte_r, Inte_rho, Inte_T, &
                                         Inte_Ye, Inte_TMeV, &
                                         Inte_cmpe, database
  CHARACTER(LEN=100)                  :: Format1, Format2, Format3, Format4
  CHARACTER(LEN=30)                   :: a
  INTEGER, DIMENSION(4)               :: LogInterp
  INTEGER                             :: i, ii, jj, datasize
  REAL(dp)                            :: Offset_cmpe
  REAL(dp), DIMENSION(2)              :: Offset_NES
   
!-------- variables for output ------------------------
  INTEGER(HID_T)                          :: file_id, group_id
  INTEGER(HSIZE_T)                        :: datasize1d(1)
  INTEGER(HSIZE_T), DIMENSION(2)          :: datasize2d
  INTEGER(HSIZE_T), DIMENSION(3)          :: datasize3d

  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE   :: InterpolantNES_nue, &
                                               InterpolantNES_nuebar

  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE   :: InterH0i, InterH0ii 
  REAL(dp), DIMENSION(:,:), ALLOCATABLE   :: SumNES_nue, SumNES_nuebar
  CHARACTER(LEN=30)                       :: outfilename = &
                                             'InterpolatedNESOutput.h5'

!-------- local variables -------------------------
  REAL(dp)                            :: sum_NES_nue, sum_NES_nuebar, &
                                         root2p, root2n
  REAL(dp), DIMENSION(Inte_nPointE)   :: buffer1, buffer2, buffer3
  REAL(dp), DIMENSION(Inte_nPointE)   :: H0i, H0ii 
  REAL(dp), DIMENSION(Inte_nPointE)   :: NES0_nue 
  REAL(dp), DIMENSION(Inte_nPointE)   :: NES0_nuebar 
  REAL(dp), DIMENSION(Inte_nPointE)   :: roots, widths
  REAL(dp)                            :: cparp = (cv+ca)**2
  REAL(dp)                            :: cparn = (cv-ca)**2

!----------------------------------------
!   interpolated energy 
!----------------------------------------

  Format1 = "(A2,I5)"
  Format2 = "(4A12)"
  Format3 = "(4ES12.3)"

  CALL AllocateGrid( Inte_E, Inte_nPointE )

  Inte_E % Unit = 'MeV                  '
  Inte_E % Name = 'Intepolated Energy   '
  Inte_E % MinValue = Inte_Emin
  Inte_E % MaxValue = Inte_Emax
  Inte_E % LogInterp = 1
  Inte_E % nPoints = Inte_nPointE

  CALL MakeLogGrid &
          ( Inte_E % MinValue, Inte_E % MaxValue, &
            Inte_E % nPoints, Inte_E % Values )

  roots = Inte_E % Values

  widths(1) = roots(1)

  DO i = 2, Inte_nPointE

    widths(i) = Inte_E % Values(i) - Inte_E % Values(i-1)

  END DO

!---------------------------------------
!    read in profile ( rho, T, Ye )
!---------------------------------------
  OPEN(1, FILE = "ProfileBruenn.d", FORM = "formatted", &
          ACTION = 'read')
  READ( 1, Format1 ) a, datasize
  READ( 1, Format2 )

  ALLOCATE( database( datasize * 4) )
  ALLOCATE( Inte_r  ( datasize ) )
  ALLOCATE( Inte_rho( datasize ) )
  ALLOCATE( Inte_T  ( datasize ) )
  ALLOCATE( Inte_Ye ( datasize ) )
  ALLOCATE( Inte_TMeV ( datasize ) )
  ALLOCATE( Inte_cmpe ( datasize ) )
  ALLOCATE( SumNES_nue( Inte_nPointE, datasize ) )
  ALLOCATE( SumNES_nuebar( Inte_nPointE, datasize ) )
  ALLOCATE( InterpolantNES_nue( Inte_nPointE, Inte_nPointE, datasize ) )
  ALLOCATE( InterpolantNES_nuebar( Inte_nPointE, Inte_nPointE, datasize ) )
  ALLOCATE( InterH0i ( Inte_nPointE, Inte_nPointE, datasize ) )
  ALLOCATE( InterH0ii( Inte_nPointE, Inte_nPointE, datasize ) )

  READ( 1, Format3 ) database
  CLOSE( 1, STATUS = 'keep')

  DO i = 1, datasize
    Inte_r(i)   = database(i*4-3)
    Inte_rho(i) = database(i*4-2)
    Inte_T(i)   = database(i*4-1)
    Inte_Ye(i)  = database(i*4)
  END DO

!---------------------------------------
!    read in the reference table
!---------------------------------------
  CALL InitializeHDF( )
  CALL ReadOpacityTableHDF( OpacityTable, &
       FileName_NES_Option = "temp_NES.h5", &
       Verbose_Option = .TRUE. )
  CALL FinalizeHDF( )

  Offset_NES = OpacityTable % Scat_NES % Offsets(1,1:2)
  Offset_cmpe = OpacityTable % EOSTable % DV % Offsets(4)

!--------------------------------------
!   do interpolation
!--------------------------------------
  
  ASSOCIATE &
  ( TableNES_H0i  => OpacityTable % Scat_NES % Kernel(1) % Values(:,:,1,:,:), &
    TableNES_H0ii => OpacityTable % Scat_NES % Kernel(1) % Values(:,:,2,:,:), &
        Tablecmpe => OpacityTable % EOSTable % DV % Variables(4) % Values, &
        Energy    => Inte_E  % Values )

  CALL LogInterpolateSingleVariable &
           ( Inte_rho, Inte_T, Inte_Ye, &
             OpacityTable % EOSTable % TS % States(1) % Values, &
             OpacityTable % EOSTable % TS % States(2) % Values, &
             OpacityTable % EOSTable % TS % States(3) % Values, &
             (/1,1,0/), Offset_cmpe, &
             Tablecmpe, &
             Inte_cmpe )

  CALL LogInterpolateSingleVariable_2D2D_Custom &
          ( LOG10(Energy), LOG10(Inte_T), &
            LOG10(Inte_cmpe / (Inte_T * kMev)), &
            LOG10(OpacityTable % EnergyGrid % Values), &
            LOG10(OpacityTable % EOSTable % TS % States(2) % Values), &
            LOG10(OpacityTable % EtaGrid % Values),    &
            Offset_NES(1), TableNES_H0i, InterH0i )

  CALL LogInterpolateSingleVariable_2D2D_Custom &
          ( LOG10(Energy), LOG10(Inte_T), &
            LOG10(Inte_cmpe / (Inte_T * kMev)), &
            LOG10(OpacityTable % EnergyGrid % Values), &
            LOG10(OpacityTable % EOSTable % TS % States(2) % Values), &
            LOG10(OpacityTable % EtaGrid % Values), &
            Offset_NES(2), TableNES_H0ii, InterH0ii )

  InterpolantNES_nue    = cparp * InterH0i + cparn * InterH0ii
  InterpolantNES_nuebar = cparn * InterH0i + cparp * InterH0ii

  DO i = 1,datasize
    DO ii = 1, Inte_nPointE
      DO jj = ii+1, Inte_nPointE
        InterpolantNES_nue(jj,ii,i) = &
          InterpolantNES_nue(ii,jj,i) * EXP( ( Energy(ii) - Energy(jj) ) / ( Inte_T(i) * kMeV) )
        InterpolantNES_nuebar(jj,ii,i) = &
          InterpolantNES_nuebar(ii,jj,i) * EXP( ( Energy(ii) - Energy(jj) ) / ( Inte_T(i) * kMeV) )
      END DO
    END DO
  END DO

  DO i = 1, datasize  
    DO ii = 1, Inte_nPointE ! e(ii)

      NES0_nue = InterpolantNES_nue(:,ii,i)
      NES0_nuebar = InterpolantNES_nuebar(:,ii,i)

      sum_NES_nue = zero
      sum_NES_nuebar = zero

      DO jj = 2, Inte_nPointE ! ep(jj)

        root2p = roots(jj-1) * roots(jj-1) * widths(jj) * 0.5_dp
        root2n = roots(jj) * roots(jj) * widths(jj) * 0.5_dp

        sum_NES_nue = sum_NES_nue + NES0_nue(jj-1) * root2p &
                                  + NES0_nue(jj)   * root2n
        sum_NES_nuebar = sum_NES_nuebar + NES0_nuebar(jj-1) * root2p &
                                        + NES0_nuebar(jj)   * root2n
      END DO ! jj

      SumNES_nue(ii,i) = sum_NES_nue
      SumNES_nuebar(ii,i) = sum_NES_nuebar

    END DO ! ii
    
  END DO ! i

  END ASSOCIATE ! Table

  CALL DeAllocateOpacityTable( OpacityTable )

!--------------------------------------
!   write out
!--------------------------------------
  CALL InitializeHDF( )
  CALL OpenFileHDF( outfilename, .true., file_id )

  CALL OpenGroupHDF( 'ProfileInfo', .true., file_id, group_id )
  datasize1d(1) = Inte_E % nPoints
  CALL WriteHDF( "Energy", Inte_E % Values(:), group_id, datasize1d )
  datasize1d(1) = datasize
  CALL WriteHDF( "Radius", Inte_r, group_id, datasize1d )
  CALL WriteHDF( "Density", Inte_rho, group_id, datasize1d )
  CALL WriteHDF( "Temperature", Inte_T, group_id, datasize1d )
  CALL WriteHDF( "Electron Fraction", Inte_ye, group_id, datasize1d )
  CALL CloseGroupHDF( group_id )

  CALL OpenGroupHDF( 'OpacitiesIMFP', .true., file_id, group_id )
  datasize2d(2) = datasize
  datasize2d(1) = Inte_E % nPoints
  CALL WriteHDF( "NES_Electron", SumNES_nue, group_id, datasize2d )
  CALL WriteHDF( "NES_ElecAnti", SumNES_nuebar, group_id, datasize2d )
  CALL CloseGroupHDF( group_id )

  CALL OpenGroupHDF( 'Opacities', .true., file_id, group_id )
  datasize3d(3) = datasize
  datasize3d(1:2) = Inte_E % nPoints
  CALL WriteHDF( "NES_Electron", InterpolantNES_nue, group_id, datasize3d )
  CALL WriteHDF( "NES_ElecAnti", InterpolantNES_nuebar, group_id, datasize3d )
  CALL CloseGroupHDF( group_id )

  CALL CloseFileHDF( file_id )
  CALL FinalizeHDF( )

  PRINT*, 'Result was written into ',outfilename

  DEALLOCATE( Inte_r, Inte_rho, Inte_T, Inte_Ye, Inte_TMeV)
  DEALLOCATE( Inte_cmpe, database )
  DEALLOCATE( InterpolantNES_nue, InterpolantNES_nuebar )
  DEALLOCATE( InterH0i, InterH0ii )

END PROGRAM wlInterpolateNES
