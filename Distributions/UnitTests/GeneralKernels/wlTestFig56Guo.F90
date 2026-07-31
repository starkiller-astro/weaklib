PROGRAM wlTestFig56Guo

  USE wlKindModule,                ONLY: dp
  USE wlEosConstantsModule,        ONLY: mp, mn, me, mmu, pi, hbarc
  USE wlGeneralLeptonScatteringModule, ONLY: &
    ProcessIndexFromReactionString, &
    gauleg, FD
  USE wlGeneralLeptonScatteringModuleThornadoInterface, ONLY: &
    CalculateAllRoutIntegrated, &
    InitGeneralScatteringKernels, &
    FinalizeGeneralScatteringKernels
  USE wlSemiLeptonicOpacityModule2D,   ONLY: &
    Opacity_CC_2D

  IMPLICIT NONE

  ! REAL(DP), PARAMETER :: conv_fac = 1.0d0 / ( (2.0d0 * pi)**3 * hbarc ) 
  REAL(DP), PARAMETER :: conv_fac = 2.0d0*pi / ( (2.0d0 * pi)**3 * hbarc ) 
  REAL(DP), PARAMETER :: tfac = 100.0d0

  !--- Guo 2020 Table II conditions ---
  ! Condition (a)
  REAL(DP), PARAMETER :: T_a      = 38.3d0
  REAL(DP), PARAMETER :: MuE_a    = 83.3d0
  REAL(DP), PARAMETER :: MuMu_a   = 64.1d0
  REAL(DP), PARAMETER :: MuN_a    = 886.0d0 ! based on weaklib's SFHo value
  REAL(DP), PARAMETER :: MuP_a    = 800.7d0 ! based on weaklib's SFHo value
  REAL(DP), PARAMETER :: Un_a     = -24.9d0 ! mean-field neutron potential [MeV] 
  REAL(DP), PARAMETER :: Up_a     = -42.9d0 ! mean-field proton  potential [MeV] 
  REAL(DP), PARAMETER :: mn_eff_a = mn      ! effective neutron  mass [MeV] 
  REAL(DP), PARAMETER :: mp_eff_a = mp      ! effective proton   mass [MeV] 
  REAL(DP), PARAMETER :: MuNue_a  = -2.1d0
  REAL(DP), PARAMETER :: MuNumu_a = -20.0d0

  ! Condition (b)
  REAL(DP), PARAMETER :: T_b      = 15.2d0
  REAL(DP), PARAMETER :: MuE_b    = 44.4d0
  REAL(DP), PARAMETER :: MuMu_b   = 37.3d0
  REAL(DP), PARAMETER :: MuN_b    = 912.7d0 ! based on weaklib's SFHo value
  REAL(DP), PARAMETER :: MuP_b    = 875.4d0 ! based on weaklib's SFHo value
  REAL(DP), PARAMETER :: Un_b     = -6.1d0  ! mean-field neutron potential [MeV] 
  REAL(DP), PARAMETER :: Up_b     = -8.9d0  ! mean-field proton  potential [MeV] 
  REAL(DP), PARAMETER :: mn_eff_b = mn      ! effective neutron  mass [MeV] 
  REAL(DP), PARAMETER :: mp_eff_b = mp      ! effective proton   mass [MeV] 
  REAL(DP), PARAMETER :: MuNue_b  = 7.1d0
  REAL(DP), PARAMETER :: MuNumu_b = 1.3d0

  INTEGER , PARAMETER :: WhichCorrection = 3   ! LO+WM+PS+FF (same convention as Guo)
  INTEGER , PARAMETER :: iProcessMin = 1
  INTEGER , PARAMETER :: iProcessMax = 34
  REAL(DP)            :: Rout_Int(iProcessMax - iProcessMin + 1)

  !--- Quadrature ---
  INTEGER, PARAMETER :: nTheta = 24      ! costheta GL points for opacity
  INTEGER, PARAMETER :: nE3    = 128     ! E3 GL points for opacity

  REAL(DP) :: E3_min = 0.0d0     ! MeV
  REAL(DP) :: E3_max = 250.0d0  ! MeV
  REAL(DP) :: xa_E3(nE3), wa_E3(nE3)

  !--- Scan parameters ---
  INTEGER , PARAMETER :: nE1 = 128    ! number of E_nu points for opacity
  REAL(DP), PARAMETER :: E1_min = 0.1d0     ! MeV
  REAL(DP), PARAMETER :: E1_max = 250.0d0   ! MeV
  REAL(DP), PARAMETER :: log_E1_min = LOG(E1_min)
  REAL(DP), PARAMETER :: log_E1_max = LOG(E1_max)
  REAL(DP)            :: E1_arr(nE1)
  INTEGER             :: iProcess_Fig5_2, iProcess_Fig5_3, iProcess_Fig5_4
  INTEGER             :: iProcess_Fig6_2, iProcess_Fig6_4, iProcess_Fig6_5
  REAL(DP)            :: Rout_avg_2, Rout_avg_3, Rout_avg_4, Rout_avg_5

  !--- Variables ---
  INTEGER  :: iE1, iE3
  REAL(DP) :: E1, E3, wE3
  INTEGER  :: iProcess
  
  ! Opacity arrays for Figure 5 (4 processes)
  REAL(DP) :: chi_numu_n_mu_p, chi_numu_e_numu_e, chi_numu_e_nue_mu, chi_nue_e_numu_mu
  
  ! Opacity arrays for Figure 6 (4 processes)
  REAL(DP) :: chi1, chi2, chi3, chi4, chi5

  CHARACTER(LEN=64), DIMENSION(4) :: process_string_Fig5
  CHARACTER(LEN=64), DIMENSION(5) :: process_string_Fig6

  !--- Open output files ---
  OPEN(UNIT=10, FILE='Guo_Fig5_cond_a.dat', STATUS='REPLACE', ACTION='WRITE')
  OPEN(UNIT=11, FILE='Guo_Fig5_cond_b.dat', STATUS='REPLACE', ACTION='WRITE')
  OPEN(UNIT=12, FILE='Guo_Fig6_cond_a.dat', STATUS='REPLACE', ACTION='WRITE')
  OPEN(UNIT=13, FILE='Guo_Fig6_cond_b.dat', STATUS='REPLACE', ACTION='WRITE')

  ! Initialize process strings
  process_string_Fig5(1) = 'nu_mu + n -> mu- + p'
  process_string_Fig5(2) = 'nu_mu + e- -> nu_mu + e-'
  process_string_Fig5(3) = 'nu_mu + e- -> nu_e + mu-'
  process_string_Fig5(4) = 'nu_bar_e + e- -> nu_bar_mu + mu-'

  process_string_Fig6(1) = 'nu_bar_e + p -> e+ + n'
  process_string_Fig6(2) = 'nu_bar_e + e- -> nu_bar_e + e-'
  process_string_Fig6(2) = 'nu_bar_mu + e- -> nu_bar_mu + e-' ! I think they messed up BIG in Guo...
  process_string_Fig6(3) = 'nu_bar_e + e- + p -> n'
  process_string_Fig6(4) = 'nu_bar_e + e- -> nu_bar_mu + mu-'
  process_string_Fig6(5) = 'nu_bar_e + e- + nu_mu -> mu-'

  !--- Write Headers ---
  ! Fig 5 Condition (a)
  WRITE(10,'(A)') '# Opacities [cm^-1] vs Incoming Neutrino Energy [MeV]'
  WRITE(10,'(A)') '# Condition (a):'
  WRITE(10,'(A)') '# Col 1: E_nu [MeV]'
  WRITE(10,'(A)') '# Col 2: ' // TRIM(ADJUSTL(process_string_Fig5(1)))
  WRITE(10,'(A)') '# Col 3: ' // TRIM(ADJUSTL(process_string_Fig5(2)))
  WRITE(10,'(A)') '# Col 4: ' // TRIM(ADJUSTL(process_string_Fig5(3)))
  WRITE(10,'(A)') '# Col 5: ' // TRIM(ADJUSTL(process_string_Fig5(4)))
  
  ! Fig 5 Condition (b)
  WRITE(11,'(A)') '# Opacities [cm^-1] vs Incoming Neutrino Energy [MeV]'
  WRITE(11,'(A)') '# Condition (b):'
  WRITE(11,'(A)') '# Col 1: E_nu [MeV]'
  WRITE(11,'(A)') '# Col 2: ' // TRIM(ADJUSTL(process_string_Fig5(1)))
  WRITE(11,'(A)') '# Col 3: ' // TRIM(ADJUSTL(process_string_Fig5(2)))
  WRITE(11,'(A)') '# Col 4: ' // TRIM(ADJUSTL(process_string_Fig5(3)))
  WRITE(11,'(A)') '# Col 5: ' // TRIM(ADJUSTL(process_string_Fig5(4)))
  
  ! Fig 6 Condition (a)
  WRITE(12,'(A)') '# Opacities [cm^-1] vs Incoming Neutrino Energy [MeV]'
  WRITE(12,'(A)') '# Condition (a):'
  WRITE(12,'(A)') '# Col 1: E_nu [MeV]'
  WRITE(12,'(A)') '# Col 2: ' // TRIM(ADJUSTL(process_string_Fig6(1)))
  WRITE(12,'(A)') '# Col 3: ' // TRIM(ADJUSTL(process_string_Fig6(2)))
  WRITE(12,'(A)') '# Col 4: ' // TRIM(ADJUSTL(process_string_Fig6(3)))
  WRITE(12,'(A)') '# Col 5: ' // TRIM(ADJUSTL(process_string_Fig6(4)))
  WRITE(12,'(A)') '# Col 6: ' // TRIM(ADJUSTL(process_string_Fig6(5)))

  ! Fig 6 Condition (b)
  WRITE(13,'(A)') '# Opacities [cm^-1] vs Incoming Neutrino Energy [MeV]'
  WRITE(13,'(A)') '# Condition (b):'
  WRITE(13,'(A)') '# Col 1: E_nu [MeV]'
  WRITE(13,'(A)') '# Col 2: ' // TRIM(ADJUSTL(process_string_Fig6(1)))
  WRITE(13,'(A)') '# Col 3: ' // TRIM(ADJUSTL(process_string_Fig6(2)))
  WRITE(13,'(A)') '# Col 4: ' // TRIM(ADJUSTL(process_string_Fig6(3)))
  WRITE(13,'(A)') '# Col 5: ' // TRIM(ADJUSTL(process_string_Fig6(4)))
  WRITE(13,'(A)') '# Col 6: ' // TRIM(ADJUSTL(process_string_Fig6(5)))

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
  CALL InitGeneralScatteringKernels(  &
      E1_arr,                         &
      xa_E3,                          &
      nE1,                            &
      nE3,                            &
      nTheta,                         &
      iProcessMin_Option=iProcessMin, &
      iProcessMax_Option=iProcessMax  )

  ! 2. nu_mu + e- -> nu_mu + e- (e- in, e- out)
  CALL ProcessIndexFromReactionString( process_string_Fig5(2), iProcess)
  iProcess_Fig5_2 = iProcess - iProcessMin + 1

  ! 4. nu_bar_e + e- -> nu_bar_mu + mu- (e- in, mu- out)
  CALL ProcessIndexFromReactionString( process_string_Fig5(3), iProcess)
  iProcess_Fig5_3 = iProcess - iProcessMin + 1

  ! 4. nu_bar_e + e- -> nu_bar_mu + mu- (e- in, mu- out)
  CALL ProcessIndexFromReactionString( process_string_Fig5(4), iProcess)
  iProcess_Fig5_4 = iProcess - iProcessMin + 1

  !============================================================================
  ! FIGURE 5: LOOP OVER ENERGY SCAN
  !============================================================================
  DO iE1 = 1, nE1
    E1 = E1_arr(iE1)

    ! --- CONDITION (A) ---
    ! 1. nu_mu + n -> mu- + p (Produces muon)
    CALL Opacity_CC_2D(WhichCorrection, 1, E1, chi_numu_n_mu_p, &
                        T_a, MuMu_a, MuN_a, MuP_a, mmu,      &
                        mn_eff_a, mp_eff_a, Un_a, Up_a, nE3)
    chi_numu_n_mu_p = chi_numu_n_mu_p * 1.0d5  
    
    Rout_avg_2 = 0.0d0
    Rout_avg_3 = 0.0d0
    Rout_avg_4 = 0.0d0
    DO iE3 = 1, nE3
      E3    = xa_E3(iE3)
      wE3   = wa_E3(iE3)
      
      CALL CalculateAllRoutIntegrated( iE1, iE3, E1, E3, T_a, MuE_a, MuMu_a, Rout_Int(:) )

      Rout_avg_2 = Rout_avg_2 + Rout_Int(iProcess_Fig5_2) * E3**2 * conv_fac * wE3
      Rout_avg_3 = Rout_avg_3 + Rout_Int(iProcess_Fig5_3) * E3**2 * conv_fac * wE3
      Rout_avg_4 = Rout_avg_4 + Rout_Int(iProcess_Fig5_4) * E3**2 * conv_fac * wE3

    ENDDO

    WRITE(10,'(5ES18.8E3)') E1, chi_numu_n_mu_p, Rout_avg_2, &
                          Rout_avg_3, Rout_avg_4

    ! --- CONDITION (B) --- 
    CALL Opacity_CC_2D(WhichCorrection, 1, E1, chi_numu_n_mu_p, &
                        T_b, MuMu_b, MuN_b, MuP_b, mmu,      &
                        mn_eff_b, mp_eff_b, Un_b, Up_b, nE3)
    chi_numu_n_mu_p = chi_numu_n_mu_p * 1.0d5  
    
    Rout_avg_2 = 0.0d0
    Rout_avg_3 = 0.0d0
    Rout_avg_4 = 0.0d0
    DO iE3 = 1, nE3
      E3    = xa_E3(iE3)
      wE3   = wa_E3(iE3)
      
      CALL CalculateAllRoutIntegrated( iE1, iE3, E1, E3, T_b, MuE_b, MuMu_b, Rout_Int(:) )

      Rout_avg_2 = Rout_avg_2 + Rout_Int(iProcess_Fig5_2) * E3**2 * conv_fac * wE3
      Rout_avg_3 = Rout_avg_3 + Rout_Int(iProcess_Fig5_3) * E3**2 * conv_fac * wE3
      Rout_avg_4 = Rout_avg_4 + Rout_Int(iProcess_Fig5_4) * E3**2 * conv_fac * wE3

    ENDDO

    WRITE(11,'(5ES18.8E3)') E1, chi_numu_n_mu_p, Rout_avg_2, &
                          Rout_avg_3, Rout_avg_4

  END DO

    ! 2. nu_mu + e- -> nu_mu + e- (e- in, e- out)
  CALL ProcessIndexFromReactionString( process_string_Fig6(2), iProcess)
  iProcess_Fig6_2 = iProcess - iProcessMin + 1

  ! 4. nu_bar_e + e- -> nu_bar_mu + mu- (e- in, mu- out)
  CALL ProcessIndexFromReactionString( process_string_Fig6(4), iProcess)
  iProcess_Fig6_4 = iProcess - iProcessMin + 1

  ! 5. nu_bar_e + e- + nu_mu -> mu-
  CALL ProcessIndexFromReactionString( process_string_Fig6(5), iProcess)
  iProcess_Fig6_5 = iProcess - iProcessMin + 1

  !============================================================================
  ! FIGURE 6: LOOP OVER ENERGY SCAN
  !============================================================================
  DO iE1 = 1, nE1
    E1 = E1_arr(iE1)

    ! --- CONDITION (A) ---
    ! 1. nu_bar_e + p -> e+ + n (Involves positron)
    CALL Opacity_CC_2D(WhichCorrection, 2, E1, chi1, &
                        T_a, MuE_a, MuN_a, MuP_a, me,      &
                        mn_eff_a, mp_eff_a, Un_a, Up_a, nE3)
    chi1 = chi1 * 1.0d5  
            
    ! 3. nu_bar_e + e- + p -> n (Involves electron)
    CALL Opacity_CC_2D(WhichCorrection, 3, E1, chi3, &
                        T_a, MuE_a, MuN_a, MuP_a, me,      &
                        mn_eff_a, mp_eff_a, Un_a, Up_a, nE3)
    chi3 = chi3 * 1.0d5  
            
    Rout_avg_2 = 0.0d0
    Rout_avg_4 = 0.0d0
    Rout_avg_5 = 0.0d0
    DO iE3 = 1, nE3
      E3    = xa_E3(iE3)
      wE3   = wa_E3(iE3)
      
      CALL CalculateAllRoutIntegrated( iE1, iE3, E1, E3, T_a, MuE_a, MuMu_a, Rout_Int(:) )

      Rout_avg_2 = Rout_avg_2 + Rout_Int(iProcess_Fig6_2) * E3**2 * conv_fac * wE3
      Rout_avg_4 = Rout_avg_4 + Rout_Int(iProcess_Fig6_4) * E3**2 * conv_fac * wE3
      Rout_avg_5 = Rout_avg_5 + Rout_Int(iProcess_Fig6_5) * E3**2 * conv_fac * wE3 * FD(E3, MuNumu_a, T_a)

    ENDDO

    WRITE(12,'(6ES18.8E3)') E1, chi1, Rout_avg_2, chi3, Rout_avg_4, Rout_avg_5

    ! --- CONDITION (B) ---
    ! 1. nu_bar_e + p -> e+ + n (Involves positron)
    CALL Opacity_CC_2D(WhichCorrection, 2, E1, chi1, &
                        T_b, MuE_b, MuN_b, MuP_b, me,      &
                        mn_eff_b, mp_eff_b, Un_b, Up_b, nE3)
    chi1 = chi1 * 1.0d5  
            
    ! 3. nu_bar_e + e- + p -> n (Involves electron)
    CALL Opacity_CC_2D(WhichCorrection, 3, E1, chi3, &
                        T_b, MuE_b, MuN_b, MuP_b, me,      &
                        mn_eff_b, mp_eff_b, Un_b, Up_b, nE3)
    chi3 = chi3 * 1.0d5  
            
    Rout_avg_2 = 0.0d0
    Rout_avg_4 = 0.0d0
    Rout_avg_5 = 0.0d0
    DO iE3 = 1, nE3
      E3    = xa_E3(iE3)
      wE3   = wa_E3(iE3)
      
      CALL CalculateAllRoutIntegrated( iE1, iE3, E1, E3, T_b, MuE_b, MuMu_b, Rout_Int(:) )

      ! Notice the blocking factor
      Rout_avg_2 = Rout_avg_2 + Rout_Int(iProcess_Fig6_2) * E3**2 * conv_fac * wE3
      Rout_avg_4 = Rout_avg_4 + Rout_Int(iProcess_Fig6_4) * E3**2 * conv_fac * wE3
      Rout_avg_5 = Rout_avg_5 + Rout_Int(iProcess_Fig6_5) * E3**2 * conv_fac * wE3 * FD(E3, MuNumu_b, T_b)

    ENDDO

    WRITE(13,'(6ES18.8E3)') E1, chi1, Rout_avg_2, chi3, Rout_avg_4, Rout_avg_5
   
  END DO

  CLOSE(10)
  CLOSE(11)
  CLOSE(12)
  CLOSE(13)
  CALL FinalizeGeneralScatteringKernels()

  WRITE(*,'(A)') '=== test_Fig56_Guo completed ==='
  WRITE(*,'(A)') 'Output files written: Guo_Fig5_cond_a.dat, Guo_Fig5_cond_b.dat'
  WRITE(*,'(A)') '                      Guo_Fig6_cond_a.dat, Guo_Fig6_cond_b.dat'

END PROGRAM wlTestFig56Guo