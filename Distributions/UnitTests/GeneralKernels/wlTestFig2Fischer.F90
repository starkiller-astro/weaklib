PROGRAM wlTestFig2Fischer

  ! Fischer et al. (2020) Fig. 2: nu_mu / nu_bar_mu scattering opacities on
  ! e- and mu- (optionally adding the e+ and mu+ contributions), by direct
  ! integration of the kernel with GeneralScatteringOpacity
  ! (wlGeneralLeptonScatteringModule).

  USE wlKindModule,                 ONLY: dp
  USE wlEosConstantsModule,          ONLY: &
    me, mmu
  USE wlGeneralLeptonScatteringModule, ONLY: &
    GeneralScatteringOpacity, &
    ProcessIndexFromReactionString

  IMPLICIT NONE

  !--- Fischer 2020 Table II conditions ---
  ! Condition (a)
  REAL(DP), PARAMETER :: T_a     = 10.0d0      ! MeV
  REAL(DP), PARAMETER :: MuE_a   = 108.1d0     ! MeV
  REAL(DP), PARAMETER :: MuMu_a  = 51.7d0      ! MeV
  ! Condition (b)
  REAL(DP), PARAMETER :: T_b     = 25.0d0      ! MeV
  REAL(DP), PARAMETER :: MuE_b   = 147.4d0     ! MeV
  REAL(DP), PARAMETER :: MuMu_b  = 132.8d0     ! MeV

  ! Set final state neutrino chemical potential to 0 (free final state approximation)
  REAL(DP), PARAMETER :: MuNu_Final = 0.0d0

  !--- Quadrature ---
  INTEGER, PARAMETER :: nPhi   = 24     ! Phi GL points for opacity
  INTEGER, PARAMETER :: nTheta = 256     ! costheta GL points for opacity
  INTEGER, PARAMETER :: nE3pts = 256     ! E3 GL points for opacity

  !--- Scan parameters ---
  INTEGER , PARAMETER :: nE1_scan = 150    ! number of E_nu points for opacity
  REAL(DP), PARAMETER :: E1_min_scan = 1.0d0     ! MeV
  REAL(DP), PARAMETER :: E1_max_scan = 320.0d0   ! MeV
  LOGICAL, PARAMETER  :: FreeFinalState = .TRUE.
  LOGICAL, PARAMETER  :: DoFullIntegration = .FALSE.

  ! Add scattering on anti-leptons (e+, mu+) to each column. Anti-lepton
  ! chemical potentials are mu_{l+} = -mu_{l-} (as in the Thornado interface).
  LOGICAL, PARAMETER  :: Include_AntiLeptons = .TRUE.

  !--- Variables ---
  INTEGER  :: iE1
  REAL(DP) :: E1

  ! Opacities for the 4 columns of Fig. 2
  REAL(DP) :: chi(4)

  CHARACTER(LEN=64), DIMENSION(4) :: process_string, process_anti_string
  INTEGER :: iProcess(4), iProcess_Anti(4), i

  !--- Open output files ---
  OPEN(UNIT=10, FILE='Fischer_Fig2_cond_a.dat',  STATUS='REPLACE', ACTION='WRITE')
  OPEN(UNIT=11, FILE='Fischer_Fig2_cond_b.dat',  STATUS='REPLACE', ACTION='WRITE')

    ! Initialize process strings
  process_string(1) = 'nu_mu + mu- -> nu_mu + mu-'
  process_string(2) = 'nu_mu + e- -> nu_mu + e-'
  process_string(3) = 'nu_bar_mu + mu- -> nu_bar_mu + mu-'
  process_string(4) = 'nu_bar_mu + e- -> nu_bar_mu + e-'

  process_anti_string(1) = 'nu_mu + mu+ -> nu_mu + mu+'
  process_anti_string(2) = 'nu_mu + e+ -> nu_mu + e+'
  process_anti_string(3) = 'nu_bar_mu + mu+ -> nu_bar_mu + mu+'
  process_anti_string(4) = 'nu_bar_mu + e+ -> nu_bar_mu + e+'

  ! Write Headers
  CALL WriteHeader( 10, '# Condition (a): T=10 MeV, mu_e=108.1 MeV, mu_mu=51.7 MeV' )
  CALL WriteHeader( 11, '# Condition (b): T=25 MeV, mu_e=147.4 MeV, mu_mu=132.8 MeV' )

  DO i=1,SIZE(process_string)
    CALL ProcessIndexFromReactionString( process_string(i), iProcess(i) )
    WRITE(*,*) process_string(i), '->', iProcess(i)
    IF (Include_AntiLeptons) THEN
      CALL ProcessIndexFromReactionString( process_anti_string(i), iProcess_Anti(i) )
      WRITE(*,*) process_anti_string(i), '->', iProcess_Anti(i)
    ENDIF
  ENDDO

  !============================================================================
  ! LOOP OVER ENERGY SCAN
  !============================================================================
  DO iE1 = 1, nE1_scan
    E1 = E1_min_scan + REAL(iE1-1, dp) * (E1_max_scan - E1_min_scan) / REAL(nE1_scan-1, dp)

    ! --- CONDITION (A) ---
    CALL Fig2Opacities( E1, T_a, MuE_a, MuMu_a, chi )
    WRITE(10,'(5ES18.8E3)') E1, chi

    ! --- CONDITION (B) ---
    CALL Fig2Opacities( E1, T_b, MuE_b, MuMu_b, chi )
    WRITE(11,'(5ES18.8E3)') E1, chi

  END DO

  CLOSE(10)
  CLOSE(11)

  WRITE(*,'(A)') '=== wlTestFig2Fischer completed ==='
  WRITE(*,'(A)') 'Output files written: Fischer_Fig2_cond_a.dat, Fischer_Fig2_cond_b.dat'

CONTAINS

  ! Opacities of the 4 Fig. 2 columns (nu_mu mu, nu_mu e, nu_bar_mu mu,
  ! nu_bar_mu e), plus the anti-lepton contributions if requested
  SUBROUTINE Fig2Opacities( E1, T, MuE, MuMu, chi )

    REAL(DP), INTENT(IN)  :: E1, T, MuE, MuMu
    REAL(DP), INTENT(OUT) :: chi(4)

    REAL(DP) :: Mu_l, m_l, chi_anti
    INTEGER  :: i

    DO i = 1, 4

      ! Columns 1,3: target muons; columns 2,4: target electrons
      IF (MOD(i,2) == 1) THEN
        Mu_l = MuMu; m_l = mmu
      ELSE
        Mu_l = MuE;  m_l = me
      ENDIF

      CALL GeneralScatteringOpacity( E1, T, Mu_l, Mu_l, MuNu_Final, m_l, m_l, &
                    iProcess(i), chi(i), nE3pts, nTheta, nPhi, &
                    IsFinalStateFree_in=FreeFinalState, DoFullIntegration_in=DoFullIntegration )

      IF (Include_AntiLeptons) THEN
        CALL GeneralScatteringOpacity( E1, T, -Mu_l, -Mu_l, MuNu_Final, m_l, m_l, &
                      iProcess_Anti(i), chi_anti, nE3pts, nTheta, nPhi, &
                      IsFinalStateFree_in=FreeFinalState, DoFullIntegration_in=DoFullIntegration )
        chi(i) = chi(i) + chi_anti
      ENDIF

    END DO

  END SUBROUTINE Fig2Opacities

  SUBROUTINE WriteHeader( iUnit, Condition )

    INTEGER,          INTENT(IN) :: iUnit
    CHARACTER(LEN=*), INTENT(IN) :: Condition

    INTEGER :: i

    IF (Include_AntiLeptons) THEN
      WRITE(iUnit,'(A)') '# Opacities [cm^-1] (Including AntiLeptons) vs Incoming Neutrino Energy [MeV]'
    ELSE
      WRITE(iUnit,'(A)') '# Opacities [cm^-1] vs Incoming Neutrino Energy [MeV]'
    ENDIF
    WRITE(iUnit,'(A)') Condition
    WRITE(iUnit,'(A)') '# Col 1: E_nu [MeV]'
    DO i = 1, 4
      IF (Include_AntiLeptons) THEN
        WRITE(iUnit,'(A,I1,A)') '# Col ', i+1, ': ' // TRIM(ADJUSTL(process_string(i))) &
          // '  +  ' // TRIM(ADJUSTL(process_anti_string(i)))
      ELSE
        WRITE(iUnit,'(A,I1,A)') '# Col ', i+1, ': ' // TRIM(ADJUSTL(process_string(i)))
      ENDIF
    END DO

  END SUBROUTINE WriteHeader

END PROGRAM wlTestFig2Fischer
