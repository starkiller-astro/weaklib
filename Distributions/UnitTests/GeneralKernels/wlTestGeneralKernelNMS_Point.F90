PROGRAM wlTestGeneralKernelNMS_Point

  USE wlKindModule,         ONLY: dp
  USE wlEosConstantsModule, ONLY: &
    mp, mn, me, mmu, kmev, kmev_inv, ca, cv, pi
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
    ProcessIndexFromReactionString
  USE wlGeneralLeptonScatteringModuleThornadoInterface, ONLY: &
    CalculateAllPhout, &
    InitGeneralScatteringKernels, &
    FinalizeGeneralScatteringKernels, &
    iProcessMin_Default, &
    iProcessMax_Default
  USE HDF5

  IMPLICIT NONE

  !--------- parameters for creating energy grid ----------------------------
  INTEGER, PARAMETER     :: Inte_nPointE = 32
  REAL(dp)               :: Inte_Emin = 1.0d-1
  REAL(dp)               :: Inte_Emax = 3.0d02
  TYPE(GridType)         :: Inte_E

  REAL(dp) :: ThornadoEgrid(32) = [ &
          0.39623412_dp,   1.47876588_dp,   2.36878696_dp,   3.71783804_dp, &
          4.82698231_dp,   6.50816976_dp,   7.89038536_dp,   9.98548116_dp, &
         11.70799824_dp,  14.31890662_dp,  16.46550740_dp,  19.71922143_dp, &
         22.39431532_dp,  26.44909374_dp,  29.78279576_dp,  34.83586062_dp, &
         38.99032007_dp,  45.28744950_dp,  50.46473687_dp,  58.31221957_dp, &
         64.76415509_dp,  74.54368803_dp,  82.58409007_dp,  94.77134402_dp, &
        104.79129305_dp, 119.97904892_dp, 132.46590940_dp, 151.39289077_dp, &
        166.95401629_dp, 190.54082047_dp, 209.93309510_dp, 239.32697048_dp  &
    ]

  REAL(DP), DIMENSION(:), ALLOCATABLE :: EnergyValues

  !--- Variables for profile read ---
  REAL(DP) :: Rho_point, T_point, Ye_point, Ym_point, Mue_Point, Mum_Point

  !--- Variables for opacity read ---
  CHARACTER(LEN=256) :: EosComposeTableName

  INTEGER             :: i, ii, jj

  ! Variables to calculate my general scattering       
  TYPE(HelmTableType) :: HelmTableMuons
  TYPE(HelmTableType) :: HelmTableElectrons
  TYPE(LeptonGasType) :: LeptonGasState
  CHARACTER(LEN=64)   :: process_string

  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: Phi0Out_General_nue_lep,    &
                                           Phi0Out_General_nueb_lep, &
                                           Phi0Out_General_num_lep,   &
                                           Phi0Out_General_numb_lep
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: Phi1Out_General_nue_lep,    &
                                           Phi1Out_General_nueb_lep, &
                                           Phi1Out_General_num_lep,   &
                                           Phi1Out_General_numb_lep

  INTEGER , PARAMETER :: nL = 2
  INTEGER , PARAMETER :: iProcessMin = 1
  INTEGER , PARAMETER :: iProcessMax = 24
  REAL(DP) :: Phout(iProcessMax - iProcessMin + 1,nL)
  REAL(DP) :: Phin (iProcessMax - iProcessMin + 1,nL)
  REAL(DP) :: E1, E3, Delta_mu, exponent
  INTEGER  :: iProcess, ierr
  INTEGER, PARAMETER :: nTheta = 24
  INTEGER(8) :: count_start, count_interm, count_end, count_rate
  REAL(DP)   :: t_General

  LOGICAL  :: Add_AntiLepton_Contribution = .TRUE.
  INTEGER  :: iProcess_nue_lep, iProcess_nueb_lep, iProcess_num_lep, iProcess_numb_lep
  INTEGER  :: iProcess_nue_alep, iProcess_nueb_alep, iProcess_num_alep, iProcess_numb_alep
  
  EosComposeTableName = "BaryonsPlusPhotonsPlusLeptonsEOS.h5"
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

  ALLOCATE( EnergyValues(Inte_nPointE) )
  EnergyValues = Inte_E % Values  
  EnergyValues = ThornadoEgrid

  ! ------------------------------------------------------
  !    Prescribe thermodynamic condition
  ! ------------------------------------------------------
  Rho_point = 1.135926133227719E+13
  T_point   = 1.020624899736616E+01 * kmev_inv
  Ye_point  = 3.228589290834957E-01
  Ym_point  = 5.839648672898146E-02
  
  ! ------------------------------------------------------
  LeptonGasState % rho = Rho_point
  LeptonGasState % T   = T_point
  LeptonGasState % yL  = Ym_point
  CALL LeptonGasEOS(HelmTableMuons, LeptonGasState)
  Mum_Point = LeptonGasState % mu

  WRITE(*,*) Mum_Point

  LeptonGasState % rho = Rho_point
  LeptonGasState % T   = T_point
  LeptonGasState % yL  = Ye_point
  CALL LeptonGasEOS(HelmTableElectrons, LeptonGasState)
  Mue_Point = LeptonGasState % mu

  CALL InitGeneralScatteringKernels(  &
      EnergyValues,                &
      EnergyValues,                &
      Inte_nPointE,                   &
      Inte_nPointE,                   &
      nTheta,                         &
      iProcessMin_Option=iProcessMin, &
      iProcessMax_Option=iProcessMax  )
      
  ALLOCATE( Phi0Out_General_nue_lep ( Inte_nPointE, Inte_nPointE ) )
  ALLOCATE( Phi0Out_General_nueb_lep( Inte_nPointE, Inte_nPointE ) )
  ALLOCATE( Phi0Out_General_num_lep ( Inte_nPointE, Inte_nPointE ) )
  ALLOCATE( Phi0Out_General_numb_lep( Inte_nPointE, Inte_nPointE ) )

  ALLOCATE( Phi1Out_General_nue_lep ( Inte_nPointE, Inte_nPointE ) )
  ALLOCATE( Phi1Out_General_nueb_lep( Inte_nPointE, Inte_nPointE ) )
  ALLOCATE( Phi1Out_General_num_lep ( Inte_nPointE, Inte_nPointE ) )
  ALLOCATE( Phi1Out_General_numb_lep( Inte_nPointE, Inte_nPointE ) )

  ! Map all indices correctly
  process_string = 'nu_e + mu- -> nu_e + mu-'
  CALL ProcessIndexFromReactionString( process_string, iProcess)
  iProcess_nue_lep = iProcess - iProcessMin + 1
  process_string = 'nu_bar_e + mu- -> nu_bar_e + mu-'
  CALL ProcessIndexFromReactionString( process_string, iProcess)
  iProcess_nueb_lep = iProcess - iProcessMin + 1
  process_string = 'nu_mu + mu- -> nu_mu + mu-'
  CALL ProcessIndexFromReactionString( process_string, iProcess)
  iProcess_num_lep = iProcess - iProcessMin + 1
  process_string = 'nu_bar_mu + mu- -> nu_bar_mu + mu-'
  CALL ProcessIndexFromReactionString( process_string, iProcess)
  iProcess_numb_lep = iProcess - iProcessMin + 1
  process_string = 'nu_e + mu+ -> nu_e + mu+'
  CALL ProcessIndexFromReactionString( process_string, iProcess)
  iProcess_nue_alep = iProcess - iProcessMin + 1
  process_string = 'nu_bar_e + mu+ -> nu_bar_e + mu+'
  CALL ProcessIndexFromReactionString( process_string, iProcess)
  iProcess_nueb_alep = iProcess - iProcessMin + 1
  process_string = 'nu_mu + mu+ -> nu_mu + mu+'
  CALL ProcessIndexFromReactionString( process_string, iProcess)
  iProcess_num_alep = iProcess - iProcessMin + 1
  process_string = 'nu_bar_mu + mu+ -> nu_bar_mu + mu+'
  CALL ProcessIndexFromReactionString( process_string, iProcess)
  iProcess_numb_alep = iProcess - iProcessMin + 1

  !============================================================================
  ! FIGURE 5: LOOP OVER ENERGY SCAN
  !============================================================================
  WRITE(*,*) 'Table done, now doing General Scattering Kernels'
  CALL SYSTEM_CLOCK(count_start, count_rate=count_rate)
  ! Make sure you are using the exact same chemical potential
  DO ii = 1, Inte_nPointE
    DO jj = ii, Inte_nPointE

      E1 = EnergyValues(ii)
      E3 = EnergyValues(jj)

      CALL CalculateAllPhout( ii, jj, E1, E3, T_point * kmev, Mue_Point, &
            Mum_Point, Phout(:,:), nL )

      Delta_mu = 0.0d0
      exponent = MIN( (E3 - E1 - Delta_mu) / (T_point * kmev), 500.0d0 )
      Phin = Phout * EXP(exponent)
            
      Phi0Out_General_nue_lep(ii,jj) = Phin (iProcess_nue_lep,1)
      Phi0Out_General_nue_lep(jj,ii) = Phout(iProcess_nue_lep,1)

      Phi0Out_General_nueb_lep(ii,jj) = Phin (iProcess_nueb_lep,1)
      Phi0Out_General_nueb_lep(jj,ii) = Phout(iProcess_nueb_lep,1)

      Phi0Out_General_num_lep(ii,jj) = Phin (iProcess_num_lep,1)
      Phi0Out_General_num_lep(jj,ii) = Phout(iProcess_num_lep,1)

      Phi0Out_General_numb_lep(ii,jj) = Phin (iProcess_numb_lep,1)
      Phi0Out_General_numb_lep(jj,ii) = Phout(iProcess_numb_lep,1)

      Phi1Out_General_nue_lep(ii,jj)  = Phin (iProcess_nue_lep,2)
      Phi1Out_General_nue_lep(jj,ii)  = Phout(iProcess_nue_lep,2)

      Phi1Out_General_nueb_lep(ii,jj) = Phin (iProcess_nueb_lep,2)
      Phi1Out_General_nueb_lep(jj,ii) = Phout(iProcess_nueb_lep,2)

      Phi1Out_General_num_lep(ii,jj)  = Phin (iProcess_num_lep,2)
      Phi1Out_General_num_lep(jj,ii)  = Phout(iProcess_num_lep,2)

      Phi1Out_General_numb_lep(ii,jj) = Phin (iProcess_numb_lep,2)
      Phi1Out_General_numb_lep(jj,ii) = Phout(iProcess_numb_lep,2)

    ENDDO
  ENDDO

  CALL SYSTEM_CLOCK(count_end)
  t_General = REAL(count_end - count_start) / REAL(count_rate)

  WRITE(*,*) 't_General', t_General

  ! --- Output Routine ---
  OPEN(UNIT=66, FILE='NMS_Phi.dat', STATUS='REPLACE', ACTION='WRITE')
  WRITE(66,'(A)') '# T, Mum, E1, E3, ' // &
                   'Phi0_nue, Phi1_nue, Phi0_nueb_lep, Phi1_nueb_lep, ' // &
                   'Phi0_num, Phi1_num, Phi0_numb_lep, Phi1_numb_lep'

  DO ii = 1, Inte_nPointE
    DO jj = 1, Inte_nPointE
      
      ! Retrieve the coordinate/thermodynamic variables for this row
      E1  = EnergyValues(ii)
      E3  = EnergyValues(jj)
      
      WRITE(66, '(12(ES17.8E3, 1X))') &
        T_point * kmev, Mum_Point, E1, E3, &
        Phi0Out_General_nue_lep(ii,jj),  Phi1Out_General_nue_lep(ii,jj),     &
        Phi0Out_General_nueb_lep(ii,jj), Phi1Out_General_nueb_lep(ii,jj),  &
        Phi0Out_General_num_lep(ii,jj),  Phi1Out_General_num_lep(ii,jj),    &
        Phi0Out_General_numb_lep(ii,jj), Phi1Out_General_numb_lep(ii,jj)

    END DO
  END DO

  CLOSE(66)

  DEALLOCATE( Phi0Out_General_nue_lep,    &
            Phi0Out_General_nueb_lep,   &
            Phi0Out_General_num_lep,    &
            Phi0Out_General_numb_lep,   &
            Phi1Out_General_nue_lep,    &
            Phi1Out_General_nueb_lep,   &
            Phi1Out_General_num_lep,    &
            Phi1Out_General_numb_lep )

  DEALLOCATE( EnergyValues )

END PROGRAM wlTestGeneralKernelNMS_Point
