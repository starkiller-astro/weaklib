PROGRAM wlTestFig56Guo

  USE wlKindModule,                ONLY: dp
  USE wlEosConstantsModule,        ONLY: mp, mn, me, mmu
  USE wlGeneralLeptonScatteringModule, ONLY: &
    GeneralScatteringOpacity, &
    ProcessIndexFromReactionString
  USE wlSemiLeptonicOpacityModule2D,   ONLY: &
    Opacity_CC_2D

  IMPLICIT NONE

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

  ! Set final state neutrino chemical potential to 0 (free final state approximation)
  REAL(DP), PARAMETER :: MuNu_Final = 0.0d0
  INTEGER , PARAMETER :: WhichCorrection = 3   ! LO+WM+PS+FF (same convention as Guo)

  !--- Quadrature ---
  INTEGER, PARAMETER :: nPhi   = 28      ! Not used if FreeFinalState
  INTEGER, PARAMETER :: nTheta = 64      ! costheta GL points for opacity
  INTEGER, PARAMETER :: nE3pts = 128     ! E3 GL points for opacity

  !--- Scan parameters ---
  INTEGER , PARAMETER :: nE1_scan = 150    ! number of E_nu points for opacity
  REAL(DP), PARAMETER :: E1_min_scan = 0.5d0     ! MeV
  REAL(DP), PARAMETER :: E1_max_scan = 250.0d0   ! MeV
  LOGICAL, PARAMETER  :: FreeFinalState = .FALSE.
  LOGICAL, PARAMETER  :: DoFullIntegration = .FALSE.
  ! LOGICAL, PARAMETER  :: FreeFinalState = .FALSE.

  !--- Variables ---
  INTEGER  :: iE1
  REAL(DP) :: E1
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
  ! FIGURE 5: LOOP OVER ENERGY SCAN
  !============================================================================
  DO iE1 = 1, nE1_scan
    E1 = E1_min_scan + REAL(iE1-1, dp) * (E1_max_scan - E1_min_scan) / REAL(nE1_scan-1, dp)

    ! --- CONDITION (A) ---
    ! 1. nu_mu + n -> mu- + p (Produces muon)
    CALL Opacity_CC_2D(WhichCorrection, 1, E1, chi_numu_n_mu_p, &
                        T_a, MuMu_a, MuN_a, MuP_a, mmu,      &
                        mn_eff_a, mp_eff_a, Un_a, Up_a, nE3pts)
    chi_numu_n_mu_p = chi_numu_n_mu_p * 1.0d5  
                        
    ! 2. nu_mu + e- -> nu_mu + e- (e- in, e- out)
    CALL ProcessIndexFromReactionString( process_string_Fig5(2), iProcess)
    CALL GeneralScatteringOpacity( E1, T_a, MuE_a, MuE_a, MuNumu_a, me, me, &
                             iProcess, chi_numu_e_numu_e, nE3pts, nTheta, nPhi, &
                             IsFinalStateFree_in=FreeFinalState, &
                             DoFullIntegration_in=DoFullIntegration )

    ! 3. nu_mu + e- -> nu_e + mu- (e- in, mu- out)
    CALL ProcessIndexFromReactionString( process_string_Fig5(3), iProcess)
    CALL GeneralScatteringOpacity( E1, T_a, MuE_a, MuMu_a, MuNue_a, me, mmu, &
                             iProcess, chi_numu_e_nue_mu, nE3pts, nTheta, nPhi, &
                             IsFinalStateFree_in=FreeFinalState, &
                             DoFullIntegration_in=DoFullIntegration )

    ! 4. nu_bar_e + e- -> nu_bar_mu + mu- (e- in, mu- out)
    CALL ProcessIndexFromReactionString( process_string_Fig5(4), iProcess)
    CALL GeneralScatteringOpacity( E1, T_a, MuE_a, MuMu_a, -MuNumu_a, me, mmu, &
                             iProcess, chi_nue_e_numu_mu, nE3pts, nTheta, nPhi, &
                             IsFinalStateFree_in=FreeFinalState, &
                             DoFullIntegration_in=DoFullIntegration )

    WRITE(10,'(5ES18.8E3)') E1, chi_numu_n_mu_p, chi_numu_e_numu_e, &
                          chi_numu_e_nue_mu, chi_nue_e_numu_mu

    ! --- CONDITION (B) --- 
    ! (Same variable mapping as A, just using _b parameters)
    CALL Opacity_CC_2D(WhichCorrection, 1, E1, chi_numu_n_mu_p, &
                        T_b, MuMu_b, MuN_b, MuP_b, mmu,      &
                        mn_eff_b, mp_eff_b, Un_b, Up_b, nE3pts)
    chi_numu_n_mu_p = chi_numu_n_mu_p * 1.0d5 

    CALL ProcessIndexFromReactionString( process_string_Fig5(2), iProcess)
    CALL GeneralScatteringOpacity( E1, T_b, MuE_b, MuE_b, MuNumu_b, me, me, &
                             iProcess, chi_numu_e_numu_e, nE3pts, nTheta, nPhi, &
                             IsFinalStateFree_in=FreeFinalState, &
                             DoFullIntegration_in=DoFullIntegration )

    CALL ProcessIndexFromReactionString( process_string_Fig5(3), iProcess)
    CALL GeneralScatteringOpacity( E1, T_b, MuE_b, MuMu_b, -MuNue_b, me, mmu, &
                             iProcess, chi_numu_e_nue_mu, nE3pts, nTheta, nPhi, &
                             IsFinalStateFree_in=FreeFinalState, &
                             DoFullIntegration_in=DoFullIntegration )

    CALL ProcessIndexFromReactionString( process_string_Fig5(4), iProcess)
    CALL GeneralScatteringOpacity( E1, T_b, MuE_b, MuMu_b, MuNumu_b, me, mmu, &
                             iProcess, chi_nue_e_numu_mu, nE3pts, nTheta, nPhi, &
                             IsFinalStateFree_in=FreeFinalState, &
                             DoFullIntegration_in=DoFullIntegration )

    WRITE(11,'(5ES18.8E3)') E1, chi_numu_n_mu_p, chi_numu_e_numu_e, &
                          chi_numu_e_nue_mu, chi_nue_e_numu_mu
  END DO

  !============================================================================
  ! FIGURE 6: LOOP OVER ENERGY SCAN
  !============================================================================
  DO iE1 = 1, nE1_scan
    E1 = E1_min_scan + REAL(iE1-1, dp) * (E1_max_scan - E1_min_scan) / REAL(nE1_scan-1, dp)

    ! --- CONDITION (A) ---
    ! 1. nu_bar_e + p -> e+ + n (Involves positron)
    CALL Opacity_CC_2D(WhichCorrection, 2, E1, chi1, &
                        T_a, MuE_a, MuN_a, MuP_a, me,      &
                        mn_eff_a, mp_eff_a, Un_a, Up_a, nE3pts)
    chi1 = chi1 * 1.0d5  
                        
    ! 2. nu_bar_e + e- -> nu_bar_e + e- (e- in, e- out)
    CALL ProcessIndexFromReactionString( process_string_Fig6(2), iProcess)
    CALL GeneralScatteringOpacity( E1, T_a, MuE_a, MuE_a, MuNue_a, me, me, &
                             iProcess, chi2, nE3pts, nTheta, nPhi, &
                             IsFinalStateFree_in=FreeFinalState, &
                             DoFullIntegration_in=DoFullIntegration )

    ! 3. nu_bar_e + e- + p -> n (Involves electron)
    CALL Opacity_CC_2D(WhichCorrection, 3, E1, chi3, &
                        T_a, MuE_a, MuN_a, MuP_a, me,      &
                        mn_eff_a, mp_eff_a, Un_a, Up_a, nE3pts)
    chi3 = chi3 * 1.0d5  

    ! 4. nu_bar_e + e- -> nu_bar_mu + mu-
    CALL ProcessIndexFromReactionString( process_string_Fig6(4), iProcess)
    CALL GeneralScatteringOpacity( E1, T_a, MuE_a, MuMu_a, -MuNue_a, me, mmu, &
                             iProcess, chi4, nE3pts, nTheta, nPhi, &
                             IsFinalStateFree_in=FreeFinalState, &
                             DoFullIntegration_in=DoFullIntegration )

    ! 5. 'nu_bar_e + e- + nu_mu -> mu-'
    iProcess = 33
    CALL GeneralScatteringOpacity( E1, T_a, MuE_a, MuMu_a, MuNuMu_a, me, mmu, &
                                   iProcess, chi5, nE3pts, nTheta )

    WRITE(12,'(6ES18.8E3)') E1, chi1, chi2, chi3, chi4, chi5

    ! --- CONDITION (B) ---
    ! 1. nu_bar_e + p -> e+ + n (Involves positron)
    CALL Opacity_CC_2D(WhichCorrection, 2, E1, chi1, &
                        T_b, MuE_b, MuN_b, MuP_b, me,      &
                        mn_eff_b, mp_eff_b, Un_b, Up_b, nE3pts)
    chi1 = chi1 * 1.0d5  
                        
    ! 2. nu_bar_e + e- -> nu_bar_e + e- (e- in, e- out)
    CALL ProcessIndexFromReactionString( process_string_Fig6(2), iProcess)
    CALL GeneralScatteringOpacity( E1, T_b, MuE_b, MuE_b, -MuNue_b, me, me, &
                             iProcess, chi2, nE3pts, nTheta, nPhi, &
                             IsFinalStateFree_in=FreeFinalState, &
                             DoFullIntegration_in=DoFullIntegration )

    ! 3. nu_bar_e + e- + p -> n (Involves electron)
    CALL Opacity_CC_2D(WhichCorrection, 3, E1, chi3, &
                        T_b, MuE_b, MuN_b, MuP_b, me,      &
                        mn_eff_b, mp_eff_b, Un_b, Up_b, nE3pts)
    chi3 = chi3 * 1.0d5  

    ! 4. nu_bar_e + e- -> nu_bar_mu + mu-
    CALL ProcessIndexFromReactionString( process_string_Fig6(4), iProcess)
    CALL GeneralScatteringOpacity( E1, T_b, MuE_b, MuMu_b, -MuNue_b, me, mmu, &
                             iProcess, chi4, nE3pts, nTheta, nPhi, &
                             IsFinalStateFree_in=FreeFinalState, &
                             DoFullIntegration_in=DoFullIntegration )

    ! 5. nu_bar_e + e- + nu_mu -> mu-
    iProcess = 33
    CALL GeneralScatteringOpacity( E1, T_b, MuE_b, MuMu_b, MuNuMu_b, me, mmu, &
                                   iProcess, chi5, nE3pts, nTheta )

    WRITE(13,'(6ES18.8E3)') E1, chi1, chi2, chi3, chi4, chi5
  END DO

  CLOSE(10)
  CLOSE(11)
  CLOSE(12)
  CLOSE(13)

  WRITE(*,'(A)') '=== test_Fig56_Guo completed ==='
  WRITE(*,'(A)') 'Output files written: Guo_Fig5_cond_a.dat, Guo_Fig5_cond_b.dat'
  WRITE(*,'(A)') '                      Guo_Fig6_cond_a.dat, Guo_Fig6_cond_b.dat'

END PROGRAM wlTestFig56Guo