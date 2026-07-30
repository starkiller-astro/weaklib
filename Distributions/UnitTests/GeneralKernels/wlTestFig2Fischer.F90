PROGRAM wlTestFig2Fischer

  USE wlKindModule,                ONLY: dp
  USE wlEosConstantsModule,        ONLY: mp, mn, me, mmu, pi, hbarc
  USE wlGeneralLeptonScatteringModule, ONLY: &
    ProcessIndexFromReactionString, &
    gauleg
  USE wlGeneralLeptonScatteringModuleThornadoInterface, ONLY: &
    CalculateAllRoutIntegrated, &
    CalculateAllPhout, &
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
  REAL(DP), PARAMETER :: T_a     = 10.0d0      ! MeV
  REAL(DP), PARAMETER :: MuE_a   = 108.1d0     ! MeV
  REAL(DP), PARAMETER :: MuMu_a  = 51.7d0      ! MeV
  ! Condition (b)
  REAL(DP), PARAMETER :: T_b     = 25.0d0      ! MeV
  REAL(DP), PARAMETER :: MuE_b   = 147.4d0     ! MeV
  REAL(DP), PARAMETER :: MuMu_b  = 132.8d0     ! MeV

  INTEGER , PARAMETER :: WhichCorrection = 3   ! LO+WM+PS+FF (same convention as Guo)
  INTEGER , PARAMETER :: iProcessMin = 1
  INTEGER , PARAMETER :: iProcessMax = 32
  INTEGER , PARAMETER :: nL = 1
  REAL(DP)            :: Rout_Int(iProcessMax - iProcessMin + 1)
  REAL(DP)            :: Phout(iProcessMax - iProcessMin + 1, nL)

  !--- Quadrature ---
  INTEGER, PARAMETER :: nTheta = 24    ! costheta GL points for opacity
  INTEGER, PARAMETER :: nE3    = 128     ! E3 GL points for opacity

  REAL(DP) :: E3_min = 0.0d0     ! MeV
  REAL(DP) :: E3_max = 320.0d0  ! MeV
  REAL(DP) :: xa_E3(nE3), wa_E3(nE3)

  !--- Scan parameters ---
  INTEGER , PARAMETER :: nE1 = 128    ! number of E_nu points for opacity
  REAL(DP), PARAMETER :: E1_min = 0.1d0     ! MeV
  REAL(DP), PARAMETER :: E1_max = 320.0d0   ! MeV
  REAL(DP), PARAMETER :: log_E1_min = LOG(E1_min)
  REAL(DP), PARAMETER :: log_E1_max = LOG(E1_max)
  REAL(DP)            :: E1_arr(nE1)
  INTEGER             :: iProcess_Fig5_2, iProcess_Fig5_3, iProcess_Fig5_4
  INTEGER             :: iProcess_Fig6_2, iProcess_Fig6_4, iProcess_Fig6_5
  REAL(DP)            :: Rout_avg_2, Rout_avg_3, Rout_avg_4, Rout_avg_5

  !--- Variables ---
  INTEGER  :: iE1, iE3
  REAL(DP) :: E1, E3, wE3
  INTEGER  :: iProcess, iProcess_1, iProcess_2, iProcess_3, iProcess_4
  
  ! Opacity arrays for Figure 6 (4 processes)
  REAL(DP) :: chi1, chi2, chi3, chi4

  CHARACTER(LEN=64), DIMENSION(4) :: process_string

  !--- Open output files ---
  OPEN(UNIT=10, FILE='Fischer_Fig2_cond_a.dat',  STATUS='REPLACE', ACTION='WRITE')
  OPEN(UNIT=11, FILE='Fischer_Fig2_cond_b.dat',  STATUS='REPLACE', ACTION='WRITE')

  ! Initialize process strings
  process_string(1) = 'nu_mu + mu- -> nu_mu + mu-'
  CALL ProcessIndexFromReactionString( process_string(1), iProcess)
  iProcess_1 = iProcess - iProcessMin + 1

  process_string(2) = 'nu_mu + e- -> nu_mu + e-'
  CALL ProcessIndexFromReactionString( process_string(2), iProcess)
  iProcess_2 = iProcess - iProcessMin + 1

  process_string(3) = 'nu_bar_mu + mu- -> nu_bar_mu + mu-'
  CALL ProcessIndexFromReactionString( process_string(3), iProcess)
  iProcess_3 = iProcess - iProcessMin + 1
  
  process_string(4) = 'nu_bar_mu + e- -> nu_bar_mu + e-'
  CALL ProcessIndexFromReactionString( process_string(4), iProcess)
  iProcess_4 = iProcess - iProcessMin + 1

  ! Write Headers
  WRITE(10,'(A)') '# Opacities [cm^-1] vs Incoming Neutrino Energy [MeV]'
  WRITE(10,'(A)') '# Condition (a): T=10 MeV, mu_e=108.1 MeV, mu_mu=51.7 MeV'
  WRITE(10,'(A)') '# Col 1: E_nu [MeV]'
  WRITE(10,'(A)') '# Col 2: ' // TRIM(ADJUSTL(process_string(1)))
  WRITE(10,'(A)') '# Col 3: ' // TRIM(ADJUSTL(process_string(2)))
  WRITE(10,'(A)') '# Col 4: ' // TRIM(ADJUSTL(process_string(3)))
  WRITE(10,'(A)') '# Col 5: ' // TRIM(ADJUSTL(process_string(4)))


  WRITE(11,'(A)') '# Opacities [cm^-1] vs Incoming Neutrino Energy [MeV]'
  WRITE(11,'(A)') '# Condition (b): T=25 MeV, mu_e=147.4 MeV, mu_mu=132.8 MeV'
  WRITE(11,'(A)') '# Col 1: E_nu [MeV]'
  WRITE(11,'(A)') '# Col 2: ' // TRIM(ADJUSTL(process_string(1)))
  WRITE(11,'(A)') '# Col 3: ' // TRIM(ADJUSTL(process_string(2)))
  WRITE(11,'(A)') '# Col 4: ' // TRIM(ADJUSTL(process_string(3)))
  WRITE(11,'(A)') '# Col 5: ' // TRIM(ADJUSTL(process_string(4)))
  
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

  !============================================================================
  ! FIGURE 2: LOOP OVER ENERGY SCAN
  !============================================================================
  DO iE1 = 1, nE1
    E1 = E1_arr(iE1)

    ! --- CONDITION (A) ---
    chi1 = 0.0d0
    chi2 = 0.0d0
    chi3 = 0.0d0
    chi4 = 0.0d0
    DO iE3 = 1, nE3
      E3    = xa_E3(iE3)
      wE3   = wa_E3(iE3)
      
      ! CALL CalculateAllRoutIntegrated( iE1, iE3, E1, E3, T_a, MuE_a, MuMu_a, Rout_Int(:) )

      ! chi1 = chi1 + Rout_Int(iProcess_1) * E3**2 * conv_fac * wE3
      ! chi2 = chi2 + Rout_Int(iProcess_2) * E3**2 * conv_fac * wE3
      ! chi3 = chi3 + Rout_Int(iProcess_3) * E3**2 * conv_fac * wE3
      ! chi4 = chi4 + Rout_Int(iProcess_4) * E3**2 * conv_fac * wE3
      
      CALL CalculateAllPhout( iE1, iE3, E1, E3, T_a, MuE_a, MuMu_a, Phout(:,:), nL )

      chi1 = chi1 + Phout(iProcess_1, 1) * E3**2 * 4.0d0*pi * wE3
      chi2 = chi2 + Phout(iProcess_2, 1) * E3**2 * 4.0d0*pi * wE3
      chi3 = chi3 + Phout(iProcess_3, 1) * E3**2 * 4.0d0*pi * wE3
      chi4 = chi4 + Phout(iProcess_4, 1) * E3**2 * 4.0d0*pi * wE3

    ENDDO

    WRITE(10,'(5ES18.8E3)') E1, chi1, chi2, chi3, chi4

    ! --- CONDITION (B) --- 
    chi1 = 0.0d0
    chi2 = 0.0d0
    chi3 = 0.0d0
    chi4 = 0.0d0
    DO iE3 = 1, nE3
      E3    = xa_E3(iE3)
      wE3   = wa_E3(iE3)
      
      ! CALL CalculateAllRoutIntegrated( iE1, iE3, E1, E3, T_b, MuE_b, MuMu_b, Rout_Int(:) )

      ! chi1 = chi1 + Rout_Int(iProcess_1) * E3**2 * conv_fac * wE3
      ! chi2 = chi2 + Rout_Int(iProcess_2) * E3**2 * conv_fac * wE3
      ! chi3 = chi3 + Rout_Int(iProcess_3) * E3**2 * conv_fac * wE3
      ! chi4 = chi4 + Rout_Int(iProcess_4) * E3**2 * conv_fac * wE3

      CALL CalculateAllPhout( iE1, iE3, E1, E3, T_b, MuE_b, MuMu_b, Phout(:,:), nL )

      chi1 = chi1 + Phout(iProcess_1, 1) * E3**2 * 4.0d0*pi * wE3
      chi2 = chi2 + Phout(iProcess_2, 1) * E3**2 * 4.0d0*pi * wE3
      chi3 = chi3 + Phout(iProcess_3, 1) * E3**2 * 4.0d0*pi * wE3
      chi4 = chi4 + Phout(iProcess_4, 1) * E3**2 * 4.0d0*pi * wE3

    ENDDO

    WRITE(11,'(5ES18.8E3)') E1, chi1, chi2, chi3, chi4

  END DO

  CLOSE(10)
  CLOSE(11)
  CALL FinalizeGeneralScatteringKernels()

  WRITE(*,'(A)') '=== test_Fig2_Fischer completed ==='
  WRITE(*,'(A)') 'Output files written: Fischer_Fig2_cond_a.dat, Fischer_Fig2_cond_b.dat'

END PROGRAM wlTestFig2Fischer