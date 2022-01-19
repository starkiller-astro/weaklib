PROGRAM wlInterpolatePair

  USE wlKindModule, ONLY: dp
  USE wlInterpolationModule, ONLY: &
    LogInterpolateSingleVariable, &
    LogInterpolateSingleVariable_2D2D_Custom
  USE wlOpacityFieldsModule, ONLY: &
    iJi0, iJii0, iJi1, iJii1
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

  !-------- variables for reading opacity table
  TYPE(OpacityTableType) :: OpacityTable
  REAL(dp)               :: Offset_cmpe 
  REAL(dp), DIMENSION(2) :: Offset_TP

  !-------- variables for reading parameters data ---------------------------
  REAL(dp), DIMENSION(:), ALLOCATABLE     :: Inte_r, Inte_rho, Inte_T, &
                                             Inte_Ye, Inte_TMeV, Inte_cmpe, &
                                             database
  CHARACTER(LEN=100)                      :: Format1, Format2, Format3
  CHARACTER(LEN=30)                       :: a
  INTEGER                                 :: i, datasize, icmpe

  !-------- variables for output -------------------------------------------
  INTEGER(HID_T)                          :: file_id, group_id
  INTEGER(HSIZE_T)                        :: datasize1d(1)
  INTEGER(HSIZE_T), DIMENSION(2)          :: datasize2d
  INTEGER(HSIZE_T), DIMENSION(3)          :: datasize3d

  REAL(dp), DIMENSION(:,:), ALLOCATABLE   :: SumPair_nue, SumPair_nuebar, &
                                             SumPair_mutau, SumPair_mutaubar

  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: InterpolantPair_nue, &
                                             InterpolantPair_mutau, &
                                             InterpolantPair_nuebar, &
                                             InterpolantPair_mutaubar

  CHARACTER(LEN=128)                      :: EosTableName, OpTableName, &
                                             ProfileName

  CHARACTER(LEN=30)                       :: outfilename = &
                                             'InterpolatedPairOutput.h5'

  !-------- local variables -------------------------------------------------
  INTEGER                                 :: ii, jj
  INTEGER, DIMENSION(4)                   :: LogInterp
  REAL(dp)                                :: cparpe  = (cv+ca)**2
  REAL(dp)                                :: cparne  = (cv-ca)**2
  REAL(dp)                                :: cparpmt = (cv+ca-2.d0)**2
  REAL(dp)                                :: cparnmt = (cv-ca)**2
  REAL(dp)                                :: root2p, root2n
  REAL(dp)                                :: sum_TP0_nue, sum_TP0_nuebar
  REAL(dp)                                :: sum_TP0_mutau, sum_TP0_mutaubar
  REAL(dp), DIMENSION(Inte_nPointE)       :: TP0_nue, TP0_mutau
  REAL(dp), DIMENSION(Inte_nPointE)       :: TP0_nuebar, TP0_mutaubar
  REAL(dp), DIMENSION(Inte_nPointE)       :: roots, widths
  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: J0i, J0ii

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
  READ(12,*) OpTableName ! Pair
  CLOSE(12, STATUS = 'keep')

  !------------------------------------------------
  !   interpolated energy 
  !------------------------------------------------

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

  !----------------------------------------------------------------------------
  !    read in profile ( rho, T, Ye )
  !----------------------------------------------------------------------------
  OPEN(1, FILE = TRIM(ProfileName), FORM = "formatted", &
          ACTION = 'read')
  READ( 1, Format1 ) a, datasize
  READ( 1, Format2 )

  ALLOCATE( database  ( datasize * 4) )
  ALLOCATE( Inte_r    ( datasize ) )
  ALLOCATE( Inte_rho  ( datasize ) )
  ALLOCATE( Inte_T    ( datasize ) )
  ALLOCATE( Inte_Ye   ( datasize ) )
  ALLOCATE( Inte_TMeV ( datasize ) )
  ALLOCATE( Inte_cmpe ( datasize ) )
  ALLOCATE( SumPair_nue       ( Inte_nPointE, datasize ) )
  ALLOCATE( SumPair_nuebar    ( Inte_nPointE, datasize ) )
  ALLOCATE( SumPair_mutau     ( Inte_nPointE, datasize ) )
  ALLOCATE( SumPair_mutaubar  ( Inte_nPointE, datasize ) )
  ALLOCATE( J0i                     ( Inte_nPointE, Inte_nPointE, datasize ) )
  ALLOCATE( J0ii                    ( Inte_nPointE, Inte_nPointE, datasize ) )
  ALLOCATE( InterpolantPair_nue     ( Inte_nPointE, Inte_nPointE, datasize ) )
  ALLOCATE( InterpolantPair_nuebar  ( Inte_nPointE, Inte_nPointE, datasize ) )
  ALLOCATE( InterpolantPair_mutau   ( Inte_nPointE, Inte_nPointE, datasize ) )
  ALLOCATE( InterpolantPair_mutaubar( Inte_nPointE, Inte_nPointE, datasize ) )

  READ( 1, Format3 ) database
  CLOSE( 1, STATUS = 'keep')

  DO i = 1, datasize
    Inte_r(i)   = database(i*4-3)
    Inte_rho(i) = database(i*4-2)
    Inte_T(i)   = database(i*4-1)
    Inte_Ye(i)  = database(i*4)
    Inte_TMeV(i) = Inte_T(i)* kMeV
  END DO

  !--------------------------------------------------------
  !    read in the reference table
  !--------------------------------------------------------
  CALL InitializeHDF( )
  CALL ReadOpacityTableHDF( OpacityTable, &
       FileName_Pair_Option               &
       = TRIM(OpTableName),               &
       EquationOfStateTableName_Option    &
       = TRIM(EosTableName),              &
       Verbose_Option = .TRUE. )
  CALL FinalizeHDF( )
 
  icmpe = OpacityTable % EOSTable % DV % Indices % iElectronChemicalPotential 
  Offset_TP   = OpacityTable % Scat_Pair  % Offsets(1,1:2)
  Offset_cmpe = OpacityTable % EOSTable % DV % Offsets(icmpe)

  !---------------------------------------------------------------------------
  !   do interpolation
  !---------------------------------------------------------------------------
  ASSOCIATE &
 ( TableTPJ0i  => OpacityTable % Scat_Pair  % Kernel(1) % Values(:,:,iJi0,:,:),  &
   TableTPJ0ii => OpacityTable % Scat_Pair  % Kernel(1) % Values(:,:,iJii0,:,:), &
   Tablecmpe   => OpacityTable % EOSTable % DV % Variables(icmpe) % Values,      &
   Energy      => Inte_E  % Values,                              &
   iEOS_Rho    => OpacityTable % EOSTable % TS % Indices % iRho, &
   iEOS_T      => OpacityTable % EOSTable % TS % Indices % iT,   &
   iEOS_Ye     => OpacityTable % EOSTable % TS % Indices % iYe,  &
   iRho        => OpacityTable % TS % Indices % iRho, &
   iT          => OpacityTable % TS % Indices % iT,   &
   iYe         => OpacityTable % TS % Indices % iYe,  &
   LogInterp   => OpacityTable % EOSTable % TS % LogInterp )

  CALL LogInterpolateSingleVariable &
       ( Inte_rho, Inte_T, Inte_Ye, &
         OpacityTable % EOSTable % TS % States(iEOS_Rho) % Values, &
         OpacityTable % EOSTable % TS % States(iEOS_T) % Values,   &
         OpacityTable % EOSTable % TS % States(iEOS_Ye) % Values,  &
         Offset_cmpe, Tablecmpe, Inte_cmpe )

  CALL LogInterpolateSingleVariable_2D2D_Custom &
       ( LOG10(Energy), LOG10(Inte_T), &
         LOG10(Inte_cmpe / Inte_TMeV), &
         LOG10(OpacityTable % EnergyGrid % Values),      &
         LOG10(OpacityTable % TS % States(iT) % Values), &
         LOG10(OpacityTable % EtaGrid % Values),         &
         Offset_TP(iJi0)  , TableTPJ0i  , J0i )

  CALL LogInterpolateSingleVariable_2D2D_Custom &
       ( LOG10(Energy), LOG10(Inte_T), &
         LOG10(Inte_cmpe / Inte_TMeV), &
         LOG10(OpacityTable % EnergyGrid % Values),      &
         LOG10(OpacityTable % TS % States(iT) % Values), &
         LOG10(OpacityTable % EtaGrid % Values),         &
         Offset_TP(iJii0)  , TableTPJ0ii  , J0ii )

  InterpolantPair_nue      = frpi * ( cparpe  * J0i + cparne  * J0ii )
  InterpolantPair_nuebar   = frpi * ( cparne  * J0i + cparpe  * J0ii )
  InterpolantPair_mutau    = frpi * ( cparpmt * J0i + cparnmt * J0ii )
  InterpolantPair_mutaubar = frpi * ( cparnmt * J0i + cparpmt * J0ii )

  DO i = 1,datasize
    DO ii = 1, Inte_nPointE
      DO jj = ii+1, Inte_nPointE
        InterpolantPair_nue(jj,ii,i)      = InterpolantPair_nuebar(ii,jj,i)
        InterpolantPair_nuebar(jj,ii,i)   = InterpolantPair_nue(ii,jj,i)
        InterpolantPair_mutau(jj,ii,i)    = InterpolantPair_mutaubar(ii,jj,i)
        InterpolantPair_mutaubar(jj,ii,i) = InterpolantPair_mutau(ii,jj,i)
      END DO
    END DO
  END DO

  DO i = 1, datasize  
    DO ii = 1, Inte_nPointE ! e(ii)

      TP0_nue      = InterpolantPair_nue(:,ii,i)
      TP0_nuebar   = InterpolantPair_nuebar(:,ii,i)
      TP0_mutau    = InterpolantPair_mutau(:,ii,i)
      TP0_mutaubar = InterpolantPair_mutaubar(:,ii,i)

      sum_TP0_nue      = zero
      sum_TP0_nuebar   = zero
      sum_TP0_mutau    = zero
      sum_TP0_mutaubar = zero

      DO jj = 2, Inte_nPointE ! ep(jj)

        root2p = roots(jj-1) * roots(jj-1) * widths(jj) * 0.5_dp
        root2n = roots(jj) * roots(jj) * widths(jj) * 0.5_dp
        root2p = root2p * EXP( - (roots(jj-1)+Energy(ii) ) / Inte_TMeV(i) )
        root2n = root2n * EXP( - (roots(jj  )+Energy(ii) ) / Inte_TMeV(i) )

        sum_TP0_nue       = sum_TP0_nue  + TP0_nue(jj-1) * root2p &
                                         + TP0_nue(jj)   * root2n

        sum_TP0_nuebar    = sum_TP0_nuebar  + TP0_nuebar(jj-1) * root2p &
                                            + TP0_nuebar(jj)   * root2n

        sum_TP0_mutau     = sum_TP0_mutau  + TP0_mutau(jj-1) * root2p &
                                           + TP0_mutau(jj)   * root2n

        sum_TP0_mutaubar  = sum_TP0_mutaubar  + TP0_mutaubar(jj-1) * root2p &
                                              + TP0_mutaubar(jj)   * root2n
      END DO ! jj

      SumPair_nue(ii,i)      = sum_TP0_nue    ! (A47) in Bruenn 85
      SumPair_nuebar(ii,i)   = sum_TP0_nuebar
      SumPair_mutau(ii,i)    = sum_TP0_mutau
      SumPair_mutaubar(ii,i) = sum_TP0_mutaubar
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
  CALL WriteHDF( "TP_Electron", SumPair_nue, group_id, datasize2d )
  CALL WriteHDF( "TP_ElecAnti", SumPair_nuebar, group_id, datasize2d )
  CALL WriteHDF( "TP_MuonTau", SumPair_mutau, group_id, datasize2d )
  CALL WriteHDF( "TP_MuTaAnti", SumPair_mutaubar, group_id, datasize2d )

  CALL CloseGroupHDF( group_id )

  CALL OpenGroupHDF( 'Opacities', .true., file_id, group_id )

  datasize3d(3) = datasize
  datasize3d(1:2) = Inte_E % nPoints
  CALL WriteHDF( "TP_Electron", InterpolantPair_nue, group_id, datasize3d )
  CALL WriteHDF( "TP_ElecAnti", InterpolantPair_nuebar, group_id, datasize3d )
  CALL WriteHDF( "TP_MuonTau",  InterpolantPair_mutau, group_id, datasize3d )
  CALL WriteHDF( "TP_MuTaAnti", InterpolantPair_mutaubar, group_id, datasize3d)

  CALL CloseGroupHDF( group_id )

  CALL CloseFileHDF( file_id )
  CALL FinalizeHDF( )

  PRINT*, 'Result was written into ',outfilename

  DEALLOCATE( Inte_r, Inte_rho, Inte_T, Inte_Ye, Inte_TMeV)
  DEALLOCATE( Inte_cmpe, database )
  DEALLOCATE( J0i, J0ii )
  DEALLOCATE( SumPair_nue, SumPair_nuebar )
  DEALLOCATE( SumPair_mutau, SumPair_mutaubar )
  DEALLOCATE( InterpolantPair_nue, InterpolantPair_nuebar )
  DEALLOCATE( InterpolantPair_mutau, InterpolantPair_mutaubar )

END PROGRAM wlInterpolatePair
