PROGRAM wlTestFig5Fischer

  ! Fischer et al. (2020) Fig. 5: opacities of the purely leptonic
  ! flavour-changing reactions (nu_e <-> nu_mu with e <-> mu), by direct
  ! integration of the kernel with GeneralScatteringOpacity
  ! (wlGeneralLeptonScatteringModule). Reactions on anti-leptons (e+, mu+)
  ! can be switched off with Include_AntiLeptons (their columns are then 0).

  USE wlKindModule,                ONLY: dp
  USE wlGeneralLeptonScatteringModule, ONLY: &
    GeneralScatteringOpacity, &
    ProcessIndexFromReactionString

  IMPLICIT NONE

  !--- Particle masses [MeV] ---
  REAL(DP), PARAMETER :: me_loc  = 0.510998950d0     ! electron mass [MeV]
  REAL(DP), PARAMETER :: mmu_loc = 105.6583755d0     ! muon mass [MeV]

  !--- Fischer 2020 Table II conditions ---
  ! Condition (a)
  REAL(DP), PARAMETER :: T_a     = 10.0d0      ! MeV
  REAL(DP), PARAMETER :: MuE_a   = 108.1d0     ! MeV
  REAL(DP), PARAMETER :: MuMu_a  = 51.3d0      ! MeV
  ! Condition (b)
  REAL(DP), PARAMETER :: T_b     = 25.0d0      ! MeV
  REAL(DP), PARAMETER :: MuE_b   = 147.4d0     ! MeV
  REAL(DP), PARAMETER :: MuMu_b  = 131.8d0     ! MeV

  ! Set final state neutrino chemical potential to 0 (free final state approximation)
  REAL(DP), PARAMETER :: MuNu_Final = 0.0d0

  !--- Quadrature ---
  INTEGER, PARAMETER :: nPhi   = 32     ! Not used if FreeFinalState
  INTEGER, PARAMETER :: nTheta = 64     ! costheta GL points for opacity
  INTEGER, PARAMETER :: nE3pts = 64     ! E3 GL points for opacity

  !--- Scan parameters ---
  INTEGER , PARAMETER :: nE1_scan = 30           ! number of E_nu points for opacity
  REAL(DP), PARAMETER :: E1_min_scan = 1.0d0     ! MeV
  REAL(DP), PARAMETER :: E1_max_scan = 260.0d0   ! MeV
  LOGICAL, PARAMETER  :: FreeFinalState = .TRUE.
  LOGICAL, PARAMETER  :: DoFullIntegration = .TRUE.

  ! Compute the reactions on anti-leptons (e+, mu+). Anti-lepton chemical
  ! potentials are mu_{l+} = -mu_{l-} (as in the Thornado interface).
  LOGICAL, PARAMETER  :: Include_AntiLeptons = .TRUE.

  !--- Reactions (one per output column) ---
  ! Particle 2 is the incoming charged lepton, particle 4 the outgoing one;
  ! they are always of different flavour (e <-> mu) and of the same charge.
  INTEGER, PARAMETER :: nProc = 8
  CHARACTER(LEN=64), PARAMETER :: process_string(nProc) = [ CHARACTER(LEN=64) :: &
    'nu_e + mu- -> nu_mu + e-',          &
    'nu_e + e+ -> nu_mu + mu+',          &
    'nu_bar_e + mu+ -> nu_bar_mu + e+',  &
    'nu_bar_e + e- -> nu_bar_mu + mu-',  &
    'nu_mu + mu+ -> nu_e + e+',          &
    'nu_mu + e- -> nu_e + mu-',          &
    'nu_bar_mu + e+ -> nu_bar_e + mu+',  &
    'nu_bar_mu + mu- -> nu_bar_e + e-' ]
  LOGICAL, PARAMETER :: l2IsMuon(nProc) = &
    [ .TRUE., .FALSE., .TRUE., .FALSE., .TRUE., .FALSE., .FALSE., .TRUE. ]
  LOGICAL, PARAMETER :: IsAntiLepton(nProc) = &
    [ .FALSE., .TRUE., .TRUE., .FALSE., .TRUE., .FALSE., .TRUE., .FALSE. ]

  !--- Variables ---
  INTEGER  :: iE1, i
  REAL(DP) :: E1
  INTEGER  :: iProcess(nProc)

  !--- Arrays to accumulate results for parallel execution ---
  ! Dimensions: (number of energy points, E_nu + nProc columns)
  REAL(DP) :: res_a(nE1_scan, nProc+1)
  REAL(DP) :: res_b(nE1_scan, nProc+1)

  !--- Open output files ---
  OPEN(UNIT=10, FILE='Fischer_Fig5_cond_a.dat',  STATUS='REPLACE', ACTION='WRITE')
  OPEN(UNIT=11, FILE='Fischer_Fig5_cond_b.dat',  STATUS='REPLACE', ACTION='WRITE')

  ! Write Headers
  CALL WriteHeader( 10, 'a', T_a, MuE_a, MuMu_a )
  CALL WriteHeader( 11, 'b', T_b, MuE_b, MuMu_b )

  DO i = 1, nProc
    CALL ProcessIndexFromReactionString( process_string(i), iProcess(i) )
    WRITE(*,*) process_string(i), '->', iProcess(i)
  END DO

  !============================================================================
  ! LOOP OVER ENERGY SCAN (PARALLELIZED)
  !============================================================================
  !$OMP PARALLEL DO PRIVATE(iE1, E1)
  DO iE1 = 1, nE1_scan
    E1 = E1_min_scan + REAL(iE1-1, dp) * (E1_max_scan - E1_min_scan) / REAL(nE1_scan-1, dp)

    ! --- CONDITION (A) ---
    res_a(iE1, 1) = E1
    CALL Fig5Opacities( E1, T_a, MuE_a, MuMu_a, res_a(iE1, 2:) )

    ! --- CONDITION (B) ---
    res_b(iE1, 1) = E1
    CALL Fig5Opacities( E1, T_b, MuE_b, MuMu_b, res_b(iE1, 2:) )

  END DO
  !$OMP END PARALLEL DO

  !============================================================================
  ! SEQUENTIAL WRITE TO FILES
  !============================================================================
  DO iE1 = 1, nE1_scan
    WRITE(10,'(9ES18.8E3)') res_a(iE1, :)
    WRITE(11,'(9ES18.8E3)') res_b(iE1, :)
  END DO

  CLOSE(10)
  CLOSE(11)

  WRITE(*,'(A)') '=== wlTestFig5Fischer completed ==='
  WRITE(*,'(A)') 'Output files written: Fischer_Fig5_cond_a.dat, Fischer_Fig5_cond_b.dat'

CONTAINS

  ! Opacities of all reactions at one energy and thermodynamic condition
  SUBROUTINE Fig5Opacities( E1, T, MuE, MuMu, chi )

    REAL(DP), INTENT(IN)  :: E1, T, MuE, MuMu
    REAL(DP), INTENT(OUT) :: chi(nProc)

    REAL(DP) :: Mu2, Mu4, m2, m4, ChemSign
    INTEGER  :: i

    DO i = 1, nProc

      IF (IsAntiLepton(i) .AND. .NOT. Include_AntiLeptons) THEN
        chi(i) = 0.0d0
        CYCLE
      ENDIF

      ChemSign = MERGE( -1.0d0, 1.0d0, IsAntiLepton(i) )
      IF (l2IsMuon(i)) THEN
        Mu2 = ChemSign * MuMu; m2 = mmu_loc
        Mu4 = ChemSign * MuE;  m4 = me_loc
      ELSE
        Mu2 = ChemSign * MuE;  m2 = me_loc
        Mu4 = ChemSign * MuMu; m4 = mmu_loc
      ENDIF

      CALL GeneralScatteringOpacity( E1, T, Mu2, Mu4, MuNu_Final, m2, m4, &
                               iProcess(i), chi(i), nE3pts, nTheta, nPhi, &
                               IsFinalStateFree_in=FreeFinalState, &
                               DoFullIntegration_in=DoFullIntegration )

    END DO

  END SUBROUTINE Fig5Opacities

  SUBROUTINE WriteHeader( iUnit, Cond, T, MuE, MuMu )

    INTEGER,          INTENT(IN) :: iUnit
    CHARACTER(LEN=1), INTENT(IN) :: Cond
    REAL(DP),         INTENT(IN) :: T, MuE, MuMu

    INTEGER :: i

    WRITE(iUnit,'(A)') '# Opacities [cm^-1] vs Incoming Neutrino Energy [MeV]'
    WRITE(iUnit,'(A,A,A,F5.1,A,F6.1,A,F6.1,A)') '# Condition (', Cond, '): T=', T, &
      ' MeV, mu_e=', MuE, ' MeV, mu_mu=', MuMu, ' MeV'
    IF (.NOT. Include_AntiLeptons) &
      WRITE(iUnit,'(A)') '# Reactions on anti-leptons (e+, mu+) not computed: set to 0'
    WRITE(iUnit,'(A)') '# Col 1: E_nu [MeV]'
    DO i = 1, nProc
      WRITE(iUnit,'(A,I1,A)') '# Col ', i+1, ': ' // TRIM(process_string(i))
    END DO

  END SUBROUTINE WriteHeader

END PROGRAM wlTestFig5Fischer
