PROGRAM wlTestNLSOnRealTrajectory

  USE wlKindModule,         ONLY: dp
  USE wlEosConstantsModule, ONLY: &
    mp, mn, me, mmu, kmev, kmev_inv, ca, cv, pi, hbarc
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
    ProcessIndexFromReactionString, gauleg
  USE wlGeneralLeptonScatteringModuleThornadoInterface, ONLY: &
    CalculateAllRoutIntegrated, &
    InitGeneralScatteringKernels, &
    FinalizeGeneralScatteringKernels, &
    iProcessMin_Default, &
    iProcessMax_Default
  USE HDF5

  IMPLICIT NONE

  !--------- parameters for creating energy grid ----------------------------
  REAL(DP), PARAMETER :: conv_fac = 2.0d0*pi / ( (2.0d0 * pi)**3 * hbarc ) 

  INTEGER, PARAMETER :: nTheta = 24    ! costheta GL points for opacity
  INTEGER, PARAMETER :: nE3    = 24     ! E3 GL points for opacity

  REAL(DP) :: E3_min = 0.0d0     ! MeV
  REAL(DP) :: E3_max = 320.0d0  ! MeV
  REAL(DP) :: xa_E3(nE3), wa_E3(nE3)

  INTEGER , PARAMETER :: nE1 = 24    ! number of E_nu points for opacity
  REAL(DP), PARAMETER :: E1_min = 0.1d0     ! MeV
  REAL(DP), PARAMETER :: E1_max = 320.0d0   ! MeV
  REAL(DP), PARAMETER :: log_E1_min = LOG(E1_min)
  REAL(DP), PARAMETER :: log_E1_max = LOG(E1_max)
  REAL(DP)            :: E1_arr(nE1)

  !--- Variables for profile read ---
  REAL(DP), ALLOCATABLE :: Inte_T(:) , Inte_Rho(:), Inte_Ym(:), &
                           Inte_Ye(:), Inte_Mue(:), Inte_Mum(:)
  REAL(DP) :: buffer
  INTEGER  :: nThermoPoints
  CHARACTER(LEN=256) :: ThermoConditionsName='trajectory_example.dat'

  !--- Variables for opacity read ---
  CHARACTER(LEN=256) :: OpTableName
  CHARACTER(LEN=256) :: Eos3DTableName
  CHARACTER(LEN=256) :: EosComposeTableName

  INTEGER                                 :: i, ii, jj

  ! Variables to calculate my general scattering       
  TYPE(HelmTableType) :: HelmTableMuons
  TYPE(HelmTableType) :: HelmTableElectrons
  TYPE(LeptonGasType) :: LeptonGasState
  CHARACTER(LEN=64)   :: process_string

  REAL(DP):: chi_nue_mu, chi_nueb_mu, chi_num_mu, chi_numb_mu
  REAL(DP):: chi_nue_e, chi_nueb_e, chi_num_e, chi_numb_e

  INTEGER , PARAMETER :: iProcessMin = 1
  INTEGER , PARAMETER :: iProcessMax = 24
  REAL(DP)            :: Rout_Int(iProcessMax - iProcessMin + 1)
  INTEGER  :: iE1, iE3
  REAL(DP) :: E1, E3, wE3, Delta_mu, exponent
  INTEGER  :: iProcess, ierr
  INTEGER(8) :: count_start, count_interm, count_end, count_rate
  REAL(DP)   :: t_General

  LOGICAL  :: Add_AntiLepton_Contribution = .TRUE.
  INTEGER  :: iProcess_nue_mum, iProcess_nueb_mum, iProcess_num_mum, iProcess_numb_mum
  INTEGER  :: iProcess_nue_mup, iProcess_nueb_mup, iProcess_num_mup, iProcess_numb_mup
  INTEGER  :: iProcess_nue_em, iProcess_nueb_em, iProcess_num_em, iProcess_numb_em
  INTEGER  :: iProcess_nue_ep, iProcess_nueb_ep, iProcess_num_ep, iProcess_numb_ep
  
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

  !============================================================================
  ! PRE-COMPUTE LIN-SPACED ENERGY GRID
  !============================================================================
  DO iE1 = 1, nE1
    E1_arr(iE1) = E1_min + REAL(iE1-1, dp) * (E1_max - E1_min) / REAL(nE1-1, dp)
  ENDDO

  !============================================================================
  ! PRE-COMPUTE LOG-SPACED ENERGY GRID
  !============================================================================
  DO iE1 = 1, nE1
    E1_arr(iE1) = EXP( log_E1_min + REAL(iE1-1, dp) * (log_E1_max - log_E1_min) / REAL(nE1-1, dp) )
  END DO

  CALL gauleg(  E3_min, E3_max, xa_E3, wa_E3, nE3  )

  OPEN(UNIT=11, FILE='NMS_NES_traj.dat',  STATUS='REPLACE', ACTION='WRITE')

  ! ------------------------------------------------------
  !    read in thermodynamic conditions
  ! ------------------------------------------------------
  ! Radius, Mass, Rho, T, Ye, Ym, Mu_e, Mu_mu, Mu_n, M_n_eff, U_n, Xn_val, Mu_p, M_p_eff, U_p, Xp_val
  OPEN(UNIT=123, FILE=trim(adjustl(ThermoConditionsName)), STATUS='OLD', ACTION='READ')
  READ(123,*) nThermoPoints
  READ(123,*)

  ALLOCATE( Inte_T(nThermoPoints)  , Inte_Rho(nThermoPoints) , &
            Inte_Ye(nThermoPoints) , Inte_Mue(nThermoPoints), &
            Inte_Ym(nThermoPoints) , Inte_Mum(nThermoPoints))

  DO i = 1, nThermoPoints
    
    READ(123,*) buffer, buffer, Inte_Rho(i), Inte_T(i), Inte_Ye(i), Inte_Ym(i), &
      Inte_Mue(i), Inte_Mum(i), buffer, buffer, buffer, buffer, buffer, buffer, buffer, buffer

    ! LeptonGasState % rho = Inte_Rho(i)
    ! LeptonGasState % T   = Inte_T(i)
    ! LeptonGasState % yL  = Inte_Ym(i)
    ! CALL LeptonGasEOS(HelmTableMuons, LeptonGasState)
    ! Inte_Mum(i) = LeptonGasState % mu

    ! LeptonGasState % rho = Inte_Rho(i)
    ! LeptonGasState % T   = Inte_T(i)
    ! LeptonGasState % yL  = Inte_Ye(i)
    ! CALL LeptonGasEOS(HelmTableElectrons, LeptonGasState)
    ! Inte_Mue(i) = LeptonGasState % mu

  END DO
  CLOSE(123)

  CALL InitGeneralScatteringKernels(  &
      E1_arr,                         &
      xa_E3,                          &
      nE1,                            &
      nE3,                            &
      nTheta,                         &
      iProcessMin_Option=iProcessMin, &
      iProcessMax_Option=iProcessMax  )

  ! Map all indices correctly
  process_string = 'nu_e + mu- -> nu_e + mu-'
  CALL ProcessIndexFromReactionString( process_string, iProcess)
  iProcess_nue_mum = iProcess
  process_string = 'nu_bar_e + mu- -> nu_bar_e + mu-'
  CALL ProcessIndexFromReactionString( process_string, iProcess)
  iProcess_nueb_mum = iProcess
  process_string = 'nu_mu + mu- -> nu_mu + mu-'
  CALL ProcessIndexFromReactionString( process_string, iProcess)
  iProcess_num_mum = iProcess
  process_string = 'nu_bar_mu + mu- -> nu_bar_mu + mu-'
  CALL ProcessIndexFromReactionString( process_string, iProcess)
  iProcess_numb_mum = iProcess
  process_string = 'nu_e + mu+ -> nu_e + mu+'
  CALL ProcessIndexFromReactionString( process_string, iProcess)
  iProcess_nue_mup = iProcess
  process_string = 'nu_bar_e + mu+ -> nu_bar_e + mu+'
  CALL ProcessIndexFromReactionString( process_string, iProcess)
  iProcess_nueb_mup = iProcess
  process_string = 'nu_mu + mu+ -> nu_mu + mu+'
  CALL ProcessIndexFromReactionString( process_string, iProcess)
  iProcess_num_mup = iProcess
  process_string = 'nu_bar_mu + mu+ -> nu_bar_mu + mu+'
  CALL ProcessIndexFromReactionString( process_string, iProcess)
  iProcess_numb_mup = iProcess

  process_string = 'nu_e + e- -> nu_e + e-'
  CALL ProcessIndexFromReactionString( process_string, iProcess)
  iProcess_nue_em = iProcess
  process_string = 'nu_bar_e + e- -> nu_bar_e + e-'
  CALL ProcessIndexFromReactionString( process_string, iProcess)
  iProcess_nueb_em = iProcess
  process_string = 'nu_mu + e- -> nu_mu + e-'
  CALL ProcessIndexFromReactionString( process_string, iProcess)
  iProcess_num_em = iProcess
  process_string = 'nu_bar_mu + e- -> nu_bar_mu + e-'
  CALL ProcessIndexFromReactionString( process_string, iProcess)
  iProcess_numb_em = iProcess
  process_string = 'nu_e + e+ -> nu_e + e+'
  CALL ProcessIndexFromReactionString( process_string, iProcess)
  iProcess_nue_ep = iProcess
  process_string = 'nu_bar_e + e+ -> nu_bar_e + e+'
  CALL ProcessIndexFromReactionString( process_string, iProcess)
  iProcess_nueb_ep = iProcess
  process_string = 'nu_mu + e+ -> nu_mu + e+'
  CALL ProcessIndexFromReactionString( process_string, iProcess)
  iProcess_num_ep = iProcess
  process_string = 'nu_bar_mu + e+ -> nu_bar_mu + e+'
  CALL ProcessIndexFromReactionString( process_string, iProcess)
  iProcess_numb_ep = iProcess

  !============================================================================
  ! FIGURE 5: LOOP OVER ENERGY SCAN
  !============================================================================
  WRITE(*,*) 'Table done, now doing General Scattering Kernels'
  CALL SYSTEM_CLOCK(count_start, count_rate=count_rate)
  DO i = 1, nThermoPoints
    ! Make sure you are using the exact same chemical potential
    WRITE(*,*) i, '/', nThermoPoints

    DO iE1 = 1, nE1
      E1 = E1_arr(iE1)

      ! --- CONDITION (A) ---
      chi_nue_mu  = 0.0d0
      chi_nueb_mu = 0.0d0
      chi_num_mu  = 0.0d0
      chi_numb_mu = 0.0d0
      chi_nue_e   = 0.0d0
      chi_nueb_e  = 0.0d0
      chi_num_e   = 0.0d0
      chi_numb_e  = 0.0d0
      DO iE3 = 1, nE3
        E3    = xa_E3(iE3)
        wE3   = wa_E3(iE3)
        
        CALL CalculateAllRoutIntegrated( iE1, iE3, E1, E3, &
          Inte_T(i), Inte_Mue(i), Inte_Mum(i), Rout_Int(:) )

        chi_nue_mu = chi_nue_mu + Rout_Int(iProcess_nue_mum) * E3**2 * conv_fac * wE3
        chi_nueb_mu = chi_nueb_mu + Rout_Int(iProcess_nueb_mum) * E3**2 * conv_fac * wE3
        chi_num_mu = chi_num_mu + Rout_Int(iProcess_num_mum) * E3**2 * conv_fac * wE3
        chi_numb_mu = chi_numb_mu + Rout_Int(iProcess_numb_mum) * E3**2 * conv_fac * wE3

        IF (Add_AntiLepton_Contribution) THEN
          chi_nue_mu = chi_nue_mu + Rout_Int(iProcess_nue_mup) * E3**2 * conv_fac * wE3
          chi_nueb_mu = chi_nueb_mu + Rout_Int(iProcess_nueb_mup) * E3**2 * conv_fac * wE3
          chi_num_mu = chi_num_mu + Rout_Int(iProcess_num_mup) * E3**2 * conv_fac * wE3
          chi_numb_mu = chi_numb_mu + Rout_Int(iProcess_numb_mup) * E3**2 * conv_fac * wE3
        END IF
        
        chi_nue_e = chi_nue_e + Rout_Int(iProcess_nue_em) * E3**2 * conv_fac * wE3
        chi_nueb_e = chi_nueb_e + Rout_Int(iProcess_nueb_em) * E3**2 * conv_fac * wE3
        chi_num_e = chi_num_e + Rout_Int(iProcess_num_em) * E3**2 * conv_fac * wE3
        chi_numb_e = chi_numb_e + Rout_Int(iProcess_numb_em) * E3**2 * conv_fac * wE3

        IF (Add_AntiLepton_Contribution) THEN
          chi_nue_e = chi_nue_e + Rout_Int(iProcess_nue_ep) * E3**2 * conv_fac * wE3
          chi_nueb_e = chi_nueb_e + Rout_Int(iProcess_nueb_ep) * E3**2 * conv_fac * wE3
          chi_num_e = chi_num_e + Rout_Int(iProcess_num_ep) * E3**2 * conv_fac * wE3
          chi_numb_e = chi_numb_e + Rout_Int(iProcess_numb_ep) * E3**2 * conv_fac * wE3
        END IF

      END DO

      WRITE(11,'(I4,9ES18.8E3)') i, E1, chi_nue_mu, chi_nueb_mu, chi_num_mu, chi_numb_mu, &
        chi_nue_e, chi_nueb_e, chi_num_e, chi_numb_e

    END DO

  END DO
  CALL SYSTEM_CLOCK(count_end)
  t_General = REAL(count_end - count_start) / REAL(count_rate)

  WRITE(*,*) 't_General', t_General

END PROGRAM wlTestNLSOnRealTrajectory
