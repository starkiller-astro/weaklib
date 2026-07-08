PROGRAM wlTestGeneralKernelVsTable

  USE wlKindModule,         ONLY: dp
  USE wlEosConstantsModule, ONLY: &
    mp, mn, me, mmu, kmev, kmev_inv, ca, cv, pi
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
    FinalizeHDF
  USE wlHelmIOModuleHDF, ONLY: &
    ReadHelmholtzTableHDF
  USE wlLeptonEOSTableModule, ONLY: &
    HelmTableType
  USE wlLeptonPhotonGasEOS, ONLY: &
    LeptonGasType, LeptonGasEOS
  USE wlOpacityTableIOModuleHDF, ONLY: &
    ReadOpacityTableHDF
  USE wlGridModule, ONLY: &
    GridType, &
    AllocateGrid, &
    DescribeGrid, &
    MakeLogGrid
  USE wlGeneralLeptonScatteringModule, ONLY: &
    CalculatePhoutPhin, &
    ProcessIndexFromReactionString
  USE HDF5

  IMPLICIT NONE

  !--------- parameters for creating energy grid ----------------------------
  INTEGER, PARAMETER     :: Inte_nPointE = 40
  REAL(dp)               :: Inte_Emin = 1.0d-1
  REAL(dp)               :: Inte_Emax = 3.0d02
  TYPE(GridType)         :: Inte_E

  !-------- variables for reading opacity table -----------------------------
  TYPE(OpacityTableType) :: OpacityTable
  REAL(dp)               :: Offset_cmpe
  REAL(dp), DIMENSION(4) :: Offset_NES

  !--- Variables for profile read ---
  REAL(DP), ALLOCATABLE :: Inte_T(:) , Inte_Rho(:), &
                           Inte_Ye(:), Inte_Mue(:), Inte_cmpe(:)
  REAL(DP) :: buffer
  INTEGER  :: nThermoPoints
  CHARACTER(LEN=256) :: ThermoConditionsName='ThermoConditions.dat'

  !--- Variables for opacity read ---
  CHARACTER(LEN=256) :: OpTableName
  CHARACTER(LEN=256) :: Eos3DTableName
  CHARACTER(LEN=256) :: EosComposeTableName

  REAL(dp), DIMENSION(:,:), ALLOCATABLE   :: SumNES_nue, SumNES_nuebar, &
                                             SumNES_mutau, SumNES_mutaubar

  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: Phi0Out_Table_nue, &
                                             Phi0Out_Table_nuebar, &
                                             Phi0Out_Table_numt, &
                                             Phi0Out_Table_numtbar

  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: Phi1Out_Table_nue, &
                                             Phi1Out_Table_nuebar, &
                                             Phi1Out_Table_numt, &
                                             Phi1Out_Table_numtbar

  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: InterH0i, InterH0ii
  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: InterH1i, InterH1ii

  INTEGER                                 :: i, ii, jj, icmpe
  INTEGER, DIMENSION(4)                   :: LogInterp
  REAL(dp)                                :: cparpe  = (cv+ca)**2
  REAL(dp)                                :: cparne  = (cv-ca)**2
  REAL(dp)                                :: cparpmt = (cv+ca-2.d0)**2
  REAL(dp)                                :: cparnmt = (cv-ca)**2

  ! Variables to calculate my general scattering       
  TYPE(HelmTableType) :: HelmTableMuons
  TYPE(HelmTableType) :: HelmTableElectrons
  TYPE(LeptonGasType) :: LeptonGasState
  CHARACTER(LEN=64)   :: process_string

  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: Phi0Out_General_nue,    &
                                             Phi0Out_General_nuebar, &
                                             Phi0Out_General_numt,   &
                                             Phi0Out_General_numtbar

  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: Phi1Out_General_nue,    &
                                             Phi1Out_General_nuebar, &
                                             Phi1Out_General_numt,   &
                                             Phi1Out_General_numtbar

  INTEGER , PARAMETER :: nL = 2
  REAL(DP) :: Phout(nL), Phin(nL)
  REAL(DP) :: E1, E3, me_loc
  INTEGER  :: iProcess, ierr
  INTEGER, PARAMETER :: nTheta = 24
  INTEGER  :: count_start, count_interm, count_end, count_rate
  REAL(DP) :: t_Table, t_General
  LOGICAL  :: UseExactTablePoints = .TRUE.

  CALL MPI_INIT( ierr )

  ! Eos3DTableName      = "3DEquationOfState.h5"
  Eos3DTableName      = "wl-EOS-SFHo-15-25-50.h5"
  EosComposeTableName = "BaryonsPlusPhotonsPlusLeptonsEOS.h5"
  OpTableName         = "wl-Op-SFHo-15-25-50-E40-NES.h5"
  CALL InitializeHDF( )
  CALL ReadHelmholtzTableHDF(    &
            HelmTableElectrons,  &
            EosComposeTableName, &
            "HelmTableElectrons" )
  CALL ReadHelmholtzTableHDF(    &
            HelmTableMuons,      &
            EosComposeTableName, &
            "HelmTableMuons" )
  CALL FinalizeHDF( )

  IF (.NOT. UseExactTablePoints) THEN
    ! ------------------------------------------------------
    !    read in thermodynamic conditions
    ! ------------------------------------------------------
    OPEN(UNIT=123, FILE=trim(adjustl(ThermoConditionsName)), STATUS='OLD', ACTION='READ')
    READ(123,*) nThermoPoints
    READ(123,*)

    ALLOCATE(Inte_T (nThermoPoints) , Inte_Rho(nThermoPoints), &
             Inte_Ye(nThermoPoints) , Inte_Mue(nThermoPoints))
    ALLOCATE( Inte_cmpe ( nThermoPoints ) )

    DO i = 1, nThermoPoints
      
      READ(123,*) Inte_T(i), Inte_Rho(i), Inte_Ye(i), buffer, &
        buffer, buffer, buffer, buffer, buffer, buffer, buffer, buffer
      Inte_T(i) = Inte_T(i) * kmev_inv

      LeptonGasState % rho = Inte_Rho(i)
      LeptonGasState % T   = Inte_T(i)
      LeptonGasState % yL  = Inte_Ye(i)
      CALL LeptonGasEOS(HelmTableElectrons, LeptonGasState)
      Inte_Mue(i) = LeptonGasState % mu

    END DO
    CLOSE(123)
  ENDIF

  !------------------------------------------------------
  !    Make energy grid
  !------------------------------------------------------
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

  !------------------------------------------------------
  !    read in the opacity table
  !------------------------------------------------------
  CALL InitializeHDF( )
  CALL ReadOpacityTableHDF( OpacityTable, &
       FileName_NES_Option                &
       = TRIM(OpTableName),               &
       EquationOfStateTableName_Option    &
       = TRIM(Eos3DTableName),              &
       Verbose_Option = .TRUE. )
  CALL FinalizeHDF( )
  
  IF (UseExactTablePoints) THEN
    nThermoPoints = 81
    nThermoPoints = 120
    ALLOCATE(Inte_T (nThermoPoints) , Inte_Rho(nThermoPoints), &
            Inte_Ye(nThermoPoints) , Inte_Mue(nThermoPoints))
    ALLOCATE( Inte_cmpe ( nThermoPoints ) )

    ! ! Fix Mu
    ! Inte_T(:) = OpacityTable % TS % States( OpacityTable % TS % Indices % iT ) % Values
    ! Inte_cmpe (:) = OpacityTable % EtaGrid % Values(80) * Inte_T(:) * kmev
    ! Inte_cmpe (:) = OpacityTable % EtaGrid % Values(40) * Inte_T(:) * kmev

    ! Fix T
    Inte_T(:) = OpacityTable % TS % States( OpacityTable % TS % Indices % iT ) % Values(45)
    Inte_cmpe (:) = OpacityTable % EtaGrid % Values(:) * Inte_T(1) * kmev

    ! USE THE SAME ENERGY GRID
    Inte_E % Values = OpacityTable % EnergyGrid % Values
  ENDIF

  ALLOCATE( InterH0i              ( Inte_nPointE, Inte_nPointE, nThermoPoints ) )
  ALLOCATE( InterH0ii             ( Inte_nPointE, Inte_nPointE, nThermoPoints ) )
  ALLOCATE( InterH1i              ( Inte_nPointE, Inte_nPointE, nThermoPoints ) )
  ALLOCATE( InterH1ii             ( Inte_nPointE, Inte_nPointE, nThermoPoints ) )
  
  ALLOCATE( Phi0Out_Table_nue      ( Inte_nPointE, Inte_nPointE, nThermoPoints ) )
  ALLOCATE( Phi0Out_Table_nuebar   ( Inte_nPointE, Inte_nPointE, nThermoPoints ) )
  ALLOCATE( Phi0Out_Table_numt     ( Inte_nPointE, Inte_nPointE, nThermoPoints ) )
  ALLOCATE( Phi0Out_Table_numtbar  ( Inte_nPointE, Inte_nPointE, nThermoPoints ) )
  ALLOCATE( Phi0Out_General_nue    ( Inte_nPointE, Inte_nPointE, nThermoPoints ) )
  ALLOCATE( Phi0Out_General_nuebar ( Inte_nPointE, Inte_nPointE, nThermoPoints ) )
  ALLOCATE( Phi0Out_General_numt   ( Inte_nPointE, Inte_nPointE, nThermoPoints ) )
  ALLOCATE( Phi0Out_General_numtbar( Inte_nPointE, Inte_nPointE, nThermoPoints ) )
  
  ALLOCATE( Phi1Out_Table_nue      ( Inte_nPointE, Inte_nPointE, nThermoPoints ) )
  ALLOCATE( Phi1Out_Table_nuebar   ( Inte_nPointE, Inte_nPointE, nThermoPoints ) )
  ALLOCATE( Phi1Out_Table_numt     ( Inte_nPointE, Inte_nPointE, nThermoPoints ) )
  ALLOCATE( Phi1Out_Table_numtbar  ( Inte_nPointE, Inte_nPointE, nThermoPoints ) )
  ALLOCATE( Phi1Out_General_nue    ( Inte_nPointE, Inte_nPointE, nThermoPoints ) )
  ALLOCATE( Phi1Out_General_nuebar ( Inte_nPointE, Inte_nPointE, nThermoPoints ) )
  ALLOCATE( Phi1Out_General_numt   ( Inte_nPointE, Inte_nPointE, nThermoPoints ) )
  ALLOCATE( Phi1Out_General_numtbar( Inte_nPointE, Inte_nPointE, nThermoPoints ) )

  icmpe = OpacityTable % EOSTable % DV % Indices % iElectronChemicalPotential
  Offset_NES = OpacityTable % Scat_NES % Offsets(1,1:4)
  Offset_cmpe = OpacityTable % EOSTable % DV % Offsets(icmpe)

  WRITE(*,*) 'Table Interpolation started'

  !----------------------------------------------------------------------------
  !   do interpolation
  !----------------------------------------------------------------------------
  CALL SYSTEM_CLOCK(count_start, count_rate=count_rate)
  ASSOCIATE &
  ( TableNES_H0i  => OpacityTable % Scat_NES % Kernel(1) % Values(:,:,iHi0,:,:),  &
    TableNES_H0ii => OpacityTable % Scat_NES % Kernel(1) % Values(:,:,iHii0,:,:), &
    TableNES_H1i  => OpacityTable % Scat_NES % Kernel(1) % Values(:,:,iHi1,:,:), &
    TableNES_H1ii => OpacityTable % Scat_NES % Kernel(1) % Values(:,:,iHii1,:,:), &
    Tablecmpe     => OpacityTable % EOSTable % DV % Variables(icmpe) % Values,    &
    Energy        => Inte_E  % Values,                              &
    iEOS_Rho      => OpacityTable % EOSTable % TS % Indices % iRho, &
    iEOS_T        => OpacityTable % EOSTable % TS % Indices % iT,   &
    iEOS_Ye       => OpacityTable % EOSTable % TS % Indices % iYe,  &
    iRho          => OpacityTable % TS % Indices % iRho, &
    iT            => OpacityTable % TS % Indices % iT,   &
    iYe           => OpacityTable % TS % Indices % iYe,  &
    LogInterp     => OpacityTable % EOSTable % TS % LogInterp )

  IF (.NOT. UseExactTablePoints) THEN
    CALL LogInterpolateSingleVariable   &
          ( Inte_rho, Inte_T, Inte_Ye, &
            OpacityTable % EOSTable % TS % States(iEOS_Rho) % Values, &
            OpacityTable % EOSTable % TS % States(iEOS_T) % Values,   &
            OpacityTable % EOSTable % TS % States(iEOS_Ye) % Values,  &
            Offset_cmpe, Tablecmpe, Inte_cmpe )
  ENDIF

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

  CALL LogInterpolateSingleVariable_2D2D_Custom &
         ( LOG10(Energy), LOG10(Inte_T),        &
           LOG10(Inte_cmpe / (Inte_T * kMev)),  &
           LOG10(OpacityTable % EnergyGrid % Values),      &
           LOG10(OpacityTable % TS % States(iT) % Values), &
           LOG10(OpacityTable % EtaGrid % Values),         &
           Offset_NES(iHi1), TableNES_H1i, InterH1i )

  CALL LogInterpolateSingleVariable_2D2D_Custom &
         ( LOG10(Energy), LOG10(Inte_T),        &
           LOG10(Inte_cmpe / (Inte_T * kMev)),  &
           LOG10(OpacityTable % EnergyGrid % Values),      &
           LOG10(OpacityTable % TS % States(iT) % Values), &
           LOG10(OpacityTable % EtaGrid % Values),         &
           Offset_NES(iHii1), TableNES_H1ii, InterH1ii )

  Phi0Out_Table_nue     = 4.0d0*pi * ( cparpe  * InterH0i + cparne  * InterH0ii )
  Phi0Out_Table_nuebar  = 4.0d0*pi * ( cparne  * InterH0i + cparpe  * InterH0ii )

  Phi0Out_Table_numt    = 4.0d0*pi * ( cparpmt * InterH0i + cparnmt * InterH0ii )
  Phi0Out_Table_numtbar = 4.0d0*pi * ( cparnmt * InterH0i + cparpmt * InterH0ii )

  Phi1Out_Table_nue     = 4.0d0*pi * ( cparpe  * InterH1i + cparne  * InterH1ii )
  Phi1Out_Table_nuebar  = 4.0d0*pi * ( cparne  * InterH1i + cparpe  * InterH1ii )

  Phi1Out_Table_numt    = 4.0d0*pi * ( cparpmt * InterH1i + cparnmt * InterH1ii )
  Phi1Out_Table_numtbar = 4.0d0*pi * ( cparnmt * InterH1i + cparpmt * InterH1ii )

  DO i = 1, nThermoPoints
    DO ii = 1, Inte_nPointE
      DO jj = ii+1, Inte_nPointE

        Phi0Out_Table_nue(jj,ii,i)          &
          = Phi0Out_Table_nue(ii,jj,i)      &
            * EXP( ( Energy(ii) - Energy(jj) ) / ( Inte_T(i) * kMeV) )

        Phi0Out_Table_nuebar(jj,ii,i)       &
          = Phi0Out_Table_nuebar(ii,jj,i)   &
            * EXP( ( Energy(ii) - Energy(jj) ) / ( Inte_T(i) * kMeV) ) 

        Phi0Out_Table_numt(jj,ii,i)        &
          = Phi0Out_Table_numt(ii,jj,i)    &
            * EXP( ( Energy(ii) - Energy(jj) ) / ( Inte_T(i) * kMeV) )

        Phi0Out_Table_numtbar(jj,ii,i)     &
          = Phi0Out_Table_numtbar(ii,jj,i) &
            * EXP( ( Energy(ii) - Energy(jj) ) / ( Inte_T(i) * kMeV) )

        Phi1Out_Table_nue(jj,ii,i)          &
          = Phi1Out_Table_nue(ii,jj,i)      &
            * EXP( ( Energy(ii) - Energy(jj) ) / ( Inte_T(i) * kMeV) )

        Phi1Out_Table_nuebar(jj,ii,i)       &
          = Phi1Out_Table_nuebar(ii,jj,i)   &
            * EXP( ( Energy(ii) - Energy(jj) ) / ( Inte_T(i) * kMeV) ) 

        Phi1Out_Table_numt(jj,ii,i)        &
          = Phi1Out_Table_numt(ii,jj,i)    &
            * EXP( ( Energy(ii) - Energy(jj) ) / ( Inte_T(i) * kMeV) )

        Phi1Out_Table_numtbar(jj,ii,i)     &
          = Phi1Out_Table_numtbar(ii,jj,i) &
            * EXP( ( Energy(ii) - Energy(jj) ) / ( Inte_T(i) * kMeV) )

      END DO
    END DO
  END DO

  END ASSOCIATE ! Table
  CALL SYSTEM_CLOCK(count_end)
  t_Table = REAL(count_end - count_start) / REAL(count_rate)

  CALL DeAllocateOpacityTable( OpacityTable )

!============================================================================
  ! FIGURE 5: LOOP OVER ENERGY SCAN
  !============================================================================
  WRITE(*,*) 'Table done, now doing General Scattering Kernels'
  CALL SYSTEM_CLOCK(count_start, count_rate=count_rate)
  DO i = 1, nThermoPoints
    ! Make sure you are using the exact same chemical potential
    Inte_Mue(i) = Inte_cmpe(i) + 0.511d0
    WRITE(*,*) i, '/', nThermoPoints
    DO ii = 1, Inte_nPointE
      DO jj = ii, Inte_nPointE

        E1 = Inte_E % Values(ii)
        E3 = Inte_E % Values(jj)
        me_loc = me
        me_loc = 0.0d0

        process_string = 'nu_e + e- -> nu_e + e-'
        CALL ProcessIndexFromReactionString( process_string, iProcess)
        CALL CalculatePhoutPhin( E1, E3, Inte_T(i) * kmev, Inte_Mue(i), Inte_Mue(i), me_loc, me_loc, &
                                iProcess, Phout(:), Phin(:), nL, nTheta_in=nTheta )
        ! Multiply by 4*pi to match what is done for the interpolated table
        Phout(:) = 4.0d0 * pi * Phout(:)
        Phin(:)  = 4.0d0 * pi * Phin(:)
        Phi0Out_General_nue    (ii,jj,i) = Phin (1)
        Phi0Out_General_nue    (jj,ii,i) = Phout(1)
        Phi1Out_General_nue    (ii,jj,i) = Phin (2)
        Phi1Out_General_nue    (jj,ii,i) = Phout(2)

        process_string = 'nu_bar_e + e- -> nu_bar_e + e-'
        CALL ProcessIndexFromReactionString( process_string, iProcess)
        CALL CalculatePhoutPhin( E1, E3, Inte_T(i) * kmev, Inte_Mue(i), Inte_Mue(i), me_loc, me_loc, &
                                iProcess, Phout(:), Phin(:), nL, nTheta_in=nTheta )
        ! Multiply by 4*pi to match what is done for the interpolated table
        Phout(:) = 4.0d0 * pi * Phout(:)
        Phin(:)  = 4.0d0 * pi * Phin(:)
        Phi0Out_General_nuebar (ii,jj,i) = Phin (1)
        Phi0Out_General_nuebar (jj,ii,i) = Phout(1)
        Phi1Out_General_nuebar (ii,jj,i) = Phin (2)
        Phi1Out_General_nuebar (jj,ii,i) = Phout(2)

        process_string = 'nu_mu + e- -> nu_mu + e-'
        CALL ProcessIndexFromReactionString( process_string, iProcess)
        CALL CalculatePhoutPhin( E1, E3, Inte_T(i) * kmev, Inte_Mue(i), Inte_Mue(i), me_loc, me_loc, &
                                iProcess, Phout(:), Phin(:), nL, nTheta_in=nTheta )
        ! Multiply by 4*pi to match what is done for the interpolated table
        Phout(:) = 4.0d0 * pi * Phout(:)
        Phin(:)  = 4.0d0 * pi * Phin(:)
        Phi0Out_General_numt   (ii,jj,i) = Phin (1)
        Phi0Out_General_numt   (jj,ii,i) = Phout(1)
        Phi1Out_General_numt   (ii,jj,i) = Phin (2)
        Phi1Out_General_numt   (jj,ii,i) = Phout(2)

        process_string = 'nu_bar_mu + e- -> nu_bar_mu + e-'
        CALL ProcessIndexFromReactionString( process_string, iProcess)
        CALL CalculatePhoutPhin( E1, E3, Inte_T(i) * kmev, Inte_Mue(i), Inte_Mue(i), me_loc, me_loc, &
                                iProcess, Phout(:), Phin(:), nL, nTheta_in=nTheta )
        ! Multiply by 4*pi to match what is done for the interpolated table
        Phout(:) = 4.0d0 * pi * Phout(:)
        Phin(:)  = 4.0d0 * pi * Phin(:)
        Phi0Out_General_numtbar(ii,jj,i) = Phin (1)
        Phi0Out_General_numtbar(jj,ii,i) = Phout(1)
        Phi1Out_General_numtbar(ii,jj,i) = Phin (2)
        Phi1Out_General_numtbar(jj,ii,i) = Phout(2)

      ENDDO
    ENDDO
  END DO
  CALL SYSTEM_CLOCK(count_end)
  t_General = REAL(count_end - count_start) / REAL(count_rate)

  ! --- Output Routine ---
  OPEN(UNIT=66, FILE='Bruenn_Table_Phi.dat', STATUS='REPLACE', ACTION='WRITE')
  OPEN(UNIT=77, FILE='General_Kernel_Phi.dat', STATUS='REPLACE', ACTION='WRITE')
  WRITE(66,'(A)') '# T, Mue, E1, E3, ' // &
                   'Phi0_nue, Phi1_nue, Phi0_nue_bar, Phi1_nue_bar, ' // &
                   'Phi0_numt, Phi1_numt, Phi0_numt_bar, Phi1_numt_bar'
  WRITE(77,'(A)') '# T, Mue, E1, E3, ' // &
                   'Phi0_nue, Phi1_nue, Phi0_nue_bar, Phi1_nue_bar, ' // &
                   'Phi0_numt, Phi1_numt, Phi0_numt_bar, Phi1_numt_bar'

  DO i = 1, nThermoPoints
    DO ii = 1, Inte_nPointE
      DO jj = 1, Inte_nPointE
        
        ! Retrieve the coordinate/thermodynamic variables for this row
        E1  = Inte_E % Values(ii)
        E3  = Inte_E % Values(jj)
        
        ! Write to Bruenn_Table_Phi.dat (UNIT=66)
        WRITE(66, '(12(ES17.8E3, 1X))') &
          Inte_T(i) * kmev, Inte_Mue(i), E1, E3, &
          Phi0Out_Table_nue(ii,jj,i),     Phi1Out_Table_nue(ii,jj,i),     &
          Phi0Out_Table_nuebar(ii,jj,i),  Phi1Out_Table_nuebar(ii,jj,i),  &
          Phi0Out_Table_numt(ii,jj,i),    Phi1Out_Table_numt(ii,jj,i),    &
          Phi0Out_Table_numtbar(ii,jj,i), Phi1Out_Table_numtbar(ii,jj,i)

        ! Write to General_Kernel_Phi.dat (UNIT=77)
        WRITE(77, '(12(ES17.8E3, 1X))') &
          Inte_T(i) * kmev, Inte_Mue(i), E1, E3, &
          Phi0Out_General_nue(ii,jj,i),     Phi1Out_General_nue(ii,jj,i),     &
          Phi0Out_General_nuebar(ii,jj,i),  Phi1Out_General_nuebar(ii,jj,i),  &
          Phi0Out_General_numt(ii,jj,i),    Phi1Out_General_numt(ii,jj,i),    &
          Phi0Out_General_numtbar(ii,jj,i), Phi1Out_General_numtbar(ii,jj,i)

      END DO
    END DO
  END DO

  ! Close the files when done
  CLOSE(66)
  CLOSE(77)
            
  DEALLOCATE( InterH0i, InterH0ii )
  DEALLOCATE( Phi0Out_Table_nue, Phi0Out_Table_nuebar, Phi0Out_Table_numt, Phi0Out_Table_numtbar)
  DEALLOCATE( Phi0Out_General_nue, Phi0Out_General_nuebar, Phi0Out_General_numt, Phi0Out_General_numtbar)

  DEALLOCATE( InterH1i, InterH1ii )
  DEALLOCATE( Phi1Out_Table_nue, Phi1Out_Table_nuebar, Phi1Out_Table_numt, Phi1Out_Table_numtbar)
  DEALLOCATE( Phi1Out_General_nue, Phi1Out_General_nuebar, Phi1Out_General_numt, Phi1Out_General_numtbar)

  WRITE(*,*) 't_Table, t_General', t_Table, t_General

END PROGRAM wlTestGeneralKernelVsTable
