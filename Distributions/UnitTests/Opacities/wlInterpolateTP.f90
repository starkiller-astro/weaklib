PROGRAM wlInterpolateTP

  USE wlKindModule, ONLY: dp
  USE wlInterpolationModule, ONLY: &
    LogInterpolateSingleVariable
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
  USE wlExtPhysicalConstantsModule, ONLY: kMeV, cvel, h
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
  REAL(dp)                            :: Offset_cmpe, Offset_TP
  
!-------- variables for output ------------------------
  INTEGER(HID_T)                          :: file_id, group_id
  INTEGER(HSIZE_T)                        :: datasize1d(1)
  INTEGER(HSIZE_T), DIMENSION(2)          :: datasize2d
  INTEGER(HSIZE_T), DIMENSION(3)          :: datasize3d

  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: InterpolantTP

  REAL(dp), DIMENSION(:,:), ALLOCATABLE   :: SumTP

  CHARACTER(LEN=30)                       :: outfilename = &
                                             'IntepolateTPOutput.h5'

!-------- local variables -------------------------
  REAL(dp)                            :: outcome_TP1, root2p, root2n
  REAL(dp), DIMENSION(Inte_nPointE)   :: buffer1, buffer2, buffer3
  REAL(dp), DIMENSION(Inte_nPointE)   :: bufferArr1
  REAL(dp), DIMENSION(Inte_nPointE)   :: roots, widths

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
  ALLOCATE( SumTP  ( Inte_nPointE, datasize ) )
  ALLOCATE( InterpolantTP  ( Inte_nPointE, Inte_nPointE, datasize ) )

  READ( 1, Format3 ) database
  WRITE(*, Format3 ) database
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
  CALL ReadOpacityTableHDF( OpacityTable, "wl-Op-SFHo-15-25-50-Chimera-TP-electronOne-momentTwo.h5" )
  CALL FinalizeHDF( )

  Offset_TP   = OpacityTable % scatt_TP  % Offsets(1,1)
  Offset_cmpe = OpacityTable % EOSTable % DV % Offsets(4)

!--------------------------------------
!   do interpolation
!--------------------------------------
  
  ASSOCIATE( TableTP   => OpacityTable % scatt_TP  % Kernel(1) % Values(:,:,:,:,1), &
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

  DO i = 1, datasize  

    Inte_TMeV(i) = Inte_T(i)* kMeV
    buffer1(:)   = Inte_T(i)
    buffer2(:)   = Inte_cmpe(i) / Inte_TMeV(i)

    DO ii = 1, Inte_nPointE ! e(ii)

    buffer3(:) = Energy(ii) ! e

      CALL LogInterpolateSingleVariable &
           ( Energy, buffer3, buffer1, buffer2, & ! interpolate ep
             OpacityTable % EnergyGrid % Values, &
             OpacityTable % EnergyGrid % Values, &
             OpacityTable % EOSTable % TS % States(2) % Values, &
             OpacityTable % EtaGrid % Values,    &
             (/1,1,1,1/), Offset_TP  , TableTP  , bufferArr1 )

      outcome_TP1  = zero

      DO jj = 2, Inte_nPointE ! ep(jj)

        root2p = roots(jj-1) * roots(jj-1) * widths(jj) * 0.5_dp
        root2n = roots(jj) * roots(jj) * widths(jj) * 0.5_dp
        root2p = root2p * EXP( - (roots(jj-1)+Energy(ii) ) / Inte_TMeV(i) )
        root2n = root2n * EXP( - (roots(jj  )+Energy(ii) ) / Inte_TMeV(i) )

        outcome_TP1  = outcome_TP1  + bufferArr1(jj-1) * root2p &
                                    + bufferArr1(jj)   * root2n
      END DO ! jj

      SumTP  (ii,i) = outcome_TP1
      InterpolantTP(:,ii,i) = bufferArr1 

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
  CALL WriteHDF( "TP_Electron", SumTP, group_id, datasize2d )
  CALL CloseGroupHDF( group_id )

  CALL OpenGroupHDF( 'Opacities', .true., file_id, group_id )
  datasize3d(3) = datasize
  datasize3d(1:2) = Inte_E % nPoints
  CALL WriteHDF( "TP_Electron", InterpolantTP, group_id, datasize3d )
  CALL CloseGroupHDF( group_id )

  CALL CloseFileHDF( file_id )
  CALL FinalizeHDF( )

  PRINT*, 'Result was written into ',outfilename

  DEALLOCATE( Inte_r, Inte_rho, Inte_T, Inte_Ye, Inte_TMeV)
  DEALLOCATE( Inte_cmpe, database )
  DEALLOCATE( InterpolantTP )

END PROGRAM wlInterpolateTP
