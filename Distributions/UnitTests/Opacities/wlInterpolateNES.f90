PROGRAM wlInterpolateNES

  USE wlKindModule, ONLY: dp
  USE wlInterpolationModule, ONLY: &
    LogInterpolateSingleVariable, &
    LogInterpolateSingleVariable_2D2D_Custom
  USE wlOpacityFieldsModule, ONLY: &
    iHi0, iHii0, iHi1, iHii1
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
  USE wlExtNumericalModule, ONLY: pi, half, frpi, zero
  USE HDF5

  IMPLICIT NONE

  !--------- parameters for creating energy grid ----------------------------
  INTEGER, PARAMETER     :: Inte_nPointE = 80
  REAL(dp)               :: Inte_Emin = 1.0d-1
  REAL(dp)               :: Inte_Emax = 3.0d02
  TYPE(GridType)         :: Inte_E

  !-------- variables for reading opacity table -----------------------------
  TYPE(OpacityTableType) :: OpacityTable
  REAL(dp)               :: Offset_cmpe
  REAL(dp), DIMENSION(2) :: Offset_NES

  !-------- variables for reading parameters data ---------------------------
  REAL(dp), DIMENSION(:), ALLOCATABLE     :: Inte_r, Inte_rho, Inte_T, &
                                             Inte_Ye, Inte_cmpe, database
  CHARACTER(LEN=100)                      :: Format1, Format2, Format3
  CHARACTER(LEN=30)                       :: a
  INTEGER                                 :: i, datasize, icmpe

  !-------- variables for output -------------------------------------------
  INTEGER(HID_T)                          :: file_id, group_id
  INTEGER(HSIZE_T)                        :: datasize1d(1)
  INTEGER(HSIZE_T), DIMENSION(2)          :: datasize2d
  INTEGER(HSIZE_T), DIMENSION(3)          :: datasize3d

  REAL(dp), DIMENSION(:,:), ALLOCATABLE   :: SumNES_nue, SumNES_nuebar, &
                                             SumNES_mutau, SumNES_mutaubar

  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: InterpolantNES_nue, &
                                             InterpolantNES_nuebar, &
                                             InterpolantNES_mutau, &
                                             InterpolantNES_mutaubar

  CHARACTER(LEN=128)                      :: EosTableName, OpTableName, &
                                             ProfileName

  CHARACTER(LEN=30)                       :: outfilename = &
                                             'InterpolatedNESOutput.h5'

  !-------- local variables -------------------------------------------------
  INTEGER                                 :: ii, jj
  INTEGER, DIMENSION(4)                   :: LogInterp
  REAL(dp)                                :: cparpe  = (cv+ca)**2
  REAL(dp)                                :: cparne  = (cv-ca)**2
  REAL(dp)                                :: cparpmt = (cv+ca-2.d0)**2
  REAL(dp)                                :: cparnmt = (cv-ca)**2
  REAL(dp)                                :: sum_NES_nue, sum_NES_nuebar, &
                                             sum_NES_mutau, sum_NES_mutaubar, &
                                             root2p, root2n
  REAL(dp), DIMENSION(Inte_nPointE)       :: NES0_nue, NES0_mutau 
  REAL(dp), DIMENSION(Inte_nPointE)       :: NES0_nuebar, NES0_mutaubar
  REAL(dp), DIMENSION(Inte_nPointE)       :: roots, widths

  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: InterH0i, InterH0ii

  !--------------------------------------------------
  !   take in profile and table names
  !--------------------------------------------------
  OPEN(12, FILE="dataList.txt", FORM = "formatted", &
           ACTION = 'read')
  READ(12,*) ProfileName
  READ(12,*) EosTableName
  READ(12,*) OpTableName ! AbEm
  READ(12,*) OpTableName ! Iso
  READ(12,*) OpTableName ! NES
  CLOSE(12, STATUS = 'keep')

  !------------------------------------------------------
  !   interpolated energy 
  !------------------------------------------------------
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

  !--------------------------------------------------------------------------
  !    read in profile ( rho, T, Ye )
  !--------------------------------------------------------------------------
  OPEN(1, FILE = TRIM(ProfileName), FORM = "formatted", &
          ACTION = 'read')
  READ( 1, Format1 ) a, datasize
  READ( 1, Format2 )

  ALLOCATE( database  ( datasize * 4) )
  ALLOCATE( Inte_r    ( datasize ) )
  ALLOCATE( Inte_rho  ( datasize ) )
  ALLOCATE( Inte_T    ( datasize ) )
  ALLOCATE( Inte_Ye   ( datasize ) )
  ALLOCATE( Inte_cmpe ( datasize ) )
  ALLOCATE( SumNES_nue     ( Inte_nPointE, datasize ) )
  ALLOCATE( SumNES_nuebar  ( Inte_nPointE, datasize ) )
  ALLOCATE( SumNES_mutau   ( Inte_nPointE, datasize ) )
  ALLOCATE( SumNES_mutaubar( Inte_nPointE, datasize ) )
  ALLOCATE( InterH0i               ( Inte_nPointE, Inte_nPointE, datasize ) )
  ALLOCATE( InterH0ii              ( Inte_nPointE, Inte_nPointE, datasize ) )
  ALLOCATE( InterpolantNES_nue     ( Inte_nPointE, Inte_nPointE, datasize ) )
  ALLOCATE( InterpolantNES_nuebar  ( Inte_nPointE, Inte_nPointE, datasize ) )
  ALLOCATE( InterpolantNES_mutau   ( Inte_nPointE, Inte_nPointE, datasize ) )
  ALLOCATE( InterpolantNES_mutaubar( Inte_nPointE, Inte_nPointE, datasize ) )

  READ( 1, Format3 ) database
  CLOSE( 1, STATUS = 'keep')

  DO i = 1, datasize
    Inte_r(i)   = database(i*4-3)
    Inte_rho(i) = database(i*4-2)
    Inte_T(i)   = database(i*4-1)
    Inte_Ye(i)  = database(i*4)
  END DO

  !------------------------------------------------------
  !    read in the opacity table
  !------------------------------------------------------
  CALL InitializeHDF( )
  CALL ReadOpacityTableHDF( OpacityTable, &
       FileName_NES_Option                &
       = TRIM(OpTableName),               &
       EquationOfStateTableName_Option    &
       = TRIM(EosTableName),              &
       Verbose_Option = .TRUE. )
  CALL FinalizeHDF( )

  icmpe = OpacityTable % EOSTable % DV % Indices % iElectronChemicalPotential
  Offset_NES = OpacityTable % Scat_NES % Offsets(1,1:2)
  Offset_cmpe = OpacityTable % EOSTable % DV % Offsets(icmpe)

  !----------------------------------------------------------------------------
  !   do interpolation
  !----------------------------------------------------------------------------
  ASSOCIATE &
  ( TableNES_H0i  => OpacityTable % Scat_NES % Kernel(1) % Values(:,:,iHi0,:,:),  &
    TableNES_H0ii => OpacityTable % Scat_NES % Kernel(1) % Values(:,:,iHii0,:,:), &
    Tablecmpe     => OpacityTable % EOSTable % DV % Variables(icmpe) % Values,    &
    Energy        => Inte_E  % Values,                              &
    iEOS_Rho      => OpacityTable % EOSTable % TS % Indices % iRho, &
    iEOS_T        => OpacityTable % EOSTable % TS % Indices % iT,   &
    iEOS_Ye       => OpacityTable % EOSTable % TS % Indices % iYe,  &
    iRho          => OpacityTable % TS % Indices % iRho, &
    iT            => OpacityTable % TS % Indices % iT,   &
    iYe           => OpacityTable % TS % Indices % iYe,  &
    LogInterp     => OpacityTable % EOSTable % TS % LogInterp )

  CALL LogInterpolateSingleVariable   &
         ( Inte_rho, Inte_T, Inte_Ye, &
           OpacityTable % EOSTable % TS % States(iEOS_Rho) % Values, &
           OpacityTable % EOSTable % TS % States(iEOS_T) % Values,   &
           OpacityTable % EOSTable % TS % States(iEOS_Ye) % Values,  &
           Offset_cmpe, Tablecmpe, Inte_cmpe )

  CALL LogInterpolateSingleVariable_2D2D_Custom &
         ( LOG10(Energy), LOG10(Inte_T),        &
           LOG10(Inte_cmpe / (Inte_T * kMev)),  &
           LOG10(OpacityTable % EnergyGrid % Values),      &
           LOG10(OpacityTable % TS % States(iT) % Values), &
           LOG10(OpacityTable % EtaGrid % Values),         &
           Offset_NES(iHi0), TableNES_H0i, InterH0i )

  CALL LogInterpolateSingleVariable_2D2D_Custom &
         ( LOG10(Energy), LOG10(Inte_T),        &
           LOG10(Inte_cmpe / (Inte_T * kMev)),  &
           LOG10(OpacityTable % EnergyGrid % Values),      &
           LOG10(OpacityTable % TS % States(iT) % Values), &
           LOG10(OpacityTable % EtaGrid % Values),         &
           Offset_NES(iHii0), TableNES_H0ii, InterH0ii )

  InterpolantNES_nue      = frpi * ( cparpe  * InterH0i + cparne  * InterH0ii )
  InterpolantNES_nuebar   = frpi * ( cparne  * InterH0i + cparpe  * InterH0ii )

  InterpolantNES_mutau    = frpi * ( cparpmt * InterH0i + cparnmt * InterH0ii )
  InterpolantNES_mutaubar = frpi * ( cparnmt * InterH0i + cparpmt * InterH0ii )

  DO i = 1,datasize
    DO ii = 1, Inte_nPointE
      DO jj = ii+1, Inte_nPointE

        InterpolantNES_nue(jj,ii,i)          &
          = InterpolantNES_nue(ii,jj,i)      &
            * EXP( ( Energy(ii) - Energy(jj) ) / ( Inte_T(i) * kMeV) )

        InterpolantNES_nuebar(jj,ii,i)       &
          = InterpolantNES_nuebar(ii,jj,i)   &
            * EXP( ( Energy(ii) - Energy(jj) ) / ( Inte_T(i) * kMeV) ) 

        InterpolantNES_mutau(jj,ii,i)        &
          = InterpolantNES_mutau(ii,jj,i)    &
            * EXP( ( Energy(ii) - Energy(jj) ) / ( Inte_T(i) * kMeV) )

        InterpolantNES_mutaubar(jj,ii,i)     &
          = InterpolantNES_mutaubar(ii,jj,i) &
            * EXP( ( Energy(ii) - Energy(jj) ) / ( Inte_T(i) * kMeV) )

      END DO
    END DO
  END DO

  DO i = 1, datasize  
    DO ii = 1, Inte_nPointE ! e(ii)

      NES0_nue      = InterpolantNES_nue(:,ii,i)
      NES0_nuebar   = InterpolantNES_nuebar(:,ii,i)
      NES0_mutau    = InterpolantNES_mutau(:,ii,i)
      NES0_mutaubar = InterpolantNES_mutaubar(:,ii,i)

      sum_NES_nue      = zero
      sum_NES_nuebar   = zero
      sum_NES_mutau    = zero
      sum_NES_mutaubar = zero

      DO jj = 2, Inte_nPointE ! ep(jj)

        root2p = roots(jj-1) * roots(jj-1) * widths(jj) * 0.5_dp
        root2n = roots(jj)   * roots(jj)   * widths(jj) * 0.5_dp

        sum_NES_nue = sum_NES_nue + NES0_nue(jj-1) * root2p &
                                  + NES0_nue(jj)   * root2n

        sum_NES_nuebar = sum_NES_nuebar + NES0_nuebar(jj-1) * root2p &
                                        + NES0_nuebar(jj)   * root2n

        sum_NES_mutau = sum_NES_mutau + NES0_mutau(jj-1) * root2p &
                                      + NES0_mutau(jj)   * root2n

        sum_NES_mutaubar = sum_NES_mutaubar + NES0_mutaubar(jj-1) * root2p &
                                            + NES0_mutaubar(jj)   * root2n
      END DO ! jj

      SumNES_nue(ii,i)      = sum_NES_nue  ! (A38) in Bruenn 85
      SumNES_nuebar(ii,i)   = sum_NES_nuebar
      SumNES_mutau(ii,i)    = sum_NES_mutau
      SumNES_mutaubar(ii,i) = sum_NES_mutaubar

    END DO ! ii
    
  END DO ! i

  END ASSOCIATE ! Table

  CALL DeAllocateOpacityTable( OpacityTable )

  !--------------------------------------------------------------------------
  !   write out
  !--------------------------------------------------------------------------
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
  CALL WriteHDF( "NES_MuonTau",  SumNES_mutau, group_id, datasize2d )
  CALL WriteHDF( "NES_MuTaAnti", SumNES_mutaubar, group_id, datasize2d )

  CALL CloseGroupHDF( group_id )

  CALL OpenGroupHDF( 'Opacities', .true., file_id, group_id )

  datasize3d(3) = datasize
  datasize3d(1:2) = Inte_E % nPoints
  CALL WriteHDF( "NES_Electron", InterpolantNES_nue, group_id, datasize3d )
  CALL WriteHDF( "NES_ElecAnti", InterpolantNES_nuebar, group_id, datasize3d )
  CALL WriteHDF( "NES_MuonTau",  InterpolantNES_mutau, group_id, datasize3d )
  CALL WriteHDF( "NES_MuTaAnti", InterpolantNES_mutaubar, group_id, datasize3d)

  CALL CloseGroupHDF( group_id )

  CALL CloseFileHDF( file_id )
  CALL FinalizeHDF( )

  PRINT*, 'Result was written into ',outfilename

  DEALLOCATE( Inte_r, Inte_rho, Inte_T, Inte_Ye )
  DEALLOCATE( Inte_cmpe, database )
  DEALLOCATE( InterpolantNES_nue, InterpolantNES_nuebar )
  DEALLOCATE( InterpolantNES_mutau, InterpolantNES_mutaubar )
  DEALLOCATE( SumNES_nue, SumNES_nuebar, SumNES_mutau, SumNES_mutaubar )
  DEALLOCATE( InterH0i, InterH0ii )

END PROGRAM wlInterpolateNES
