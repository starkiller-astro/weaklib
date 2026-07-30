!================================================================================
!!! Authors:  Implemented following Guo et al. (2020) PRD 102, 023037,
!!!           Section II.A and Appendix A [primary source],
!!!           and Fischer et al. (2020) PRD 102, 123001,
!!!           Bollig (2018) PhD thesis Section 6.3.1 for additional notation
!!!           as well as for the neutral current coefficients

!!! Date:     April 2026
!!!
!!!   General process (Guo Table I, process (c)):
!!!     nu_l2(E1) + l2(E2)  ->  nu_l4(E3) + l4(E4)
!!!   and its reverse via detailed balance:
!!!     nu_l4(E3) + l4(E4)  ->  nu_l2(E1) + l2(E2)
!!!
!!!   where l2 and l4 are any charged leptons (e, mu, tau) with masses m2, m4.
!!!   The neutrinos nu_l2 and nu_l4 are massless.
!!! Physics summary:
!!!   - Purely leptonic neutral-current (Z-boson exchange only) 1-24
!!!   - Purely leptonic charged-current (W-boson exchange only) 25-29
!!!   - Spin-summed matrix element (Guo Eq. A3 / Fischer Eq. B21):
!!!       |M|^2 = 64 G_F^2 (p1.p2)(p3.p4)
!!!     where 1=nu_l2, 2=l2, 3=nu_l4, 4=l4
!!!   - Both neutrinos massless (E=|p|); charged leptons fully relativistic
!!!   - Fermi-Dirac statistics for l2 and l4; no NR approximation
!!!   - Full kinematics: exact E_minus kinematic lower bound (Guo Eq. 4)
!!!
!!! R_out kernel (Guo Eqs. 2-6, outscattering direction):
!!!   Forward process nu_l2(E1)+l2(E2)->nu_l4(E3)+l4(E4), occupation f2*(1-f4):
!!!
!!!   R_out(E1,E3,mu) = lambda_1 * (1/16*pi*Delta^5) * (A1*I_2 + B1*I_1 + C1*I_0)
!!!                                                                    (Guo Eq. 3)
!!!
!!!   I_s = integral_{E-}^inf dE2 E2^s f2(E2) [1 - f4(E1+E2-E3)]    (Guo Eq. 6)
!!!
!!!   Delta = sqrt(E1^2 - 2*E1*E3*mu + E3^2)                         (Guo Eq. 4)
!!!   k     = Q / [E1*E3*(1-mu)],  Q = (m4^2 - m2^2)/2               (Guo Eq. 5)
!!!   E-    = (1/2)*{(E3-E1)*(1+k) + Delta*sqrt[(1+k)^2 + 2*m2^2/(E1*E3*(1-mu))]}
!!!           max(E-, m2)                                             (Guo Eq. 4)
!!!
!!! Conjugate channel for charged-current reactions: 
!!! replace xMu_l2 -> -xMu_l2, xMu_l4 -> -xMu_l4 (Guo after Eq. A5)
!!!
!!! I_s integrals can be computed using GL integration, polylogarithm, or Fukushima (2015) formulae
!!!
!!! Units: R_out/R_in in MeV^-4 (natural units hbar=c=1 absorbed throughout).
!!!   Opacity: chi [cm^-1] = (1/(2*pi*hbarc)^3) * integral R_out*E3^2 dE3 d(costh)
!!!
!!! nE3:      GL points for E3 energy integral in GeneralScatteringOpacity (default 32)
!!! nTheta: GL points for costheta integral in GeneralScatteringOpacity (default 8)
!================================================================================

MODULE wlGeneralLeptonScatteringModule

  USE wlKindModule,         ONLY: dp
  USE wlEosConstantsModule, ONLY: &
    pi, Gw_MeV, hbarc, sin2W, me, cvel
  USE wlPolylogModule, ONLY: &
    Fermi_Dirac_Integrals, &
    Init_Polylogarithms 
  USE wlFukushimaFermiDiracIntegralsModule, ONLY: &
    F0z_Fukushima, &
    F1z_Fukushima, &
    F2z_Fukushima

  USE, INTRINSIC :: iso_c_binding, ONLY: c_double

  IMPLICIT NONE
  PRIVATE

  ! 1 is direct GL quadarture integration (slowest by far)
  ! 2 is polylogarithm-based evaluation (Bollig 2018, Chapter V, Section 8)
  ! 3 is Fukushima (2015) rational approximation formulae (fastest)
  INTEGER, PARAMETER  :: CompleteFDIntegralsEvaluationMethod = 3
  ! Method 1: If you want to integrate directly only call gauleg once
  LOGICAL                :: GauLegInitialized = .FALSE.
  INTEGER , PARAMETER    :: nGL_FI = 64
  REAL(DP), DIMENSION(:) :: xa_FI(nGL_FI), wa_FI(nGL_FI)

  !====================================================================
  ! INELASTIC FLAVOR EXCHANGE & CONVERSION (Guo Eq. A3 & A4)
  !====================================================================
  ! LFE: Lepton Flavor Exchange (t-channel W-exchange)
  ! Example: nu_mu + e- <-> nu_e + mu-  (Guo Eq. A3, process c)
  REAL(DP), PARAMETER :: lam1_LFE = 64.0d0*Gw_MeV**2
  REAL(DP), PARAMETER :: lam2_LFE = 0.0d0
  REAL(DP), PARAMETER :: lam3_LFE = 0.0d0

  ! LFC: Lepton Flavor Conversion (s-channel W-exchange)
  ! Example: nu_bar_e + e- <-> nu_bar_mu + mu- (Guo Eq. A4, process d)
  REAL(DP), PARAMETER :: lam1_LFC = 0.0d0
  REAL(DP), PARAMETER :: lam2_LFC = 64.0d0*Gw_MeV**2
  REAL(DP), PARAMETER :: lam3_LFC = 0.0d0

  ! These are not in Guo 2020 but are NMS and can be found in Fischer 2020
  ! and Bollig's thesis. In principle NMS has less terms than all these processes,
  ! but maybe we can use the same machinery to compute them, it's just that all terms
  ! that depends on a mass difference will be zero. Not sure if this slows down the code (?)
  ! Bollig's thesis has alphas instead of lambdas. lambda = Gw_MeV**2/(2.0d0*pi**2) * alpha

  !=============================================================================
  ! Bollig Thesis Table 6.2: Neutrino and antineutrino scattering on electrons (e-)
  !=============================================================================
  ! nu_e + e-
  REAL(DP), PARAMETER :: CV_ee_m  = +0.5d0 + 2.0d0*sin2W
  REAL(DP), PARAMETER :: CA_ee_m  = +0.5d0
  ! anti-nu_e + e-
  REAL(DP), PARAMETER :: CV_aee_m = +0.5d0 + 2.0d0*sin2W
  REAL(DP), PARAMETER :: CA_aee_m = -0.5d0
  ! nu_mu + e-
  REAL(DP), PARAMETER :: CV_me_m  = -0.5d0 + 2.0d0*sin2W
  REAL(DP), PARAMETER :: CA_me_m  = -0.5d0
  ! anti-nu_mu + e-
  REAL(DP), PARAMETER :: CV_ame_m = -0.5d0 + 2.0d0*sin2W
  REAL(DP), PARAMETER :: CA_ame_m = +0.5d0
  ! nu_tau + e-
  REAL(DP), PARAMETER :: CV_te_m  = -0.5d0 + 2.0d0*sin2W
  REAL(DP), PARAMETER :: CA_te_m  = -0.5d0
  ! anti-nu_tau + e-
  REAL(DP), PARAMETER :: CV_ate_m = -0.5d0 + 2.0d0*sin2W
  REAL(DP), PARAMETER :: CA_ate_m = +0.5d0

  !=============================================================================
  !Bollig Thesis Table 6.3: Neutrino and antineutrino scattering on positrons (e+)
  !=============================================================================
  ! nu_e + e+
  REAL(DP), PARAMETER :: CV_ee_p  = +0.5d0 + 2.0d0*sin2W
  REAL(DP), PARAMETER :: CA_ee_p  = -0.5d0
  ! anti-nu_e + e+
  REAL(DP), PARAMETER :: CV_aee_p = +0.5d0 + 2.0d0*sin2W
  REAL(DP), PARAMETER :: CA_aee_p = +0.5d0
  ! nu_mu + e+
  REAL(DP), PARAMETER :: CV_me_p  = -0.5d0 + 2.0d0*sin2W
  REAL(DP), PARAMETER :: CA_me_p  = +0.5d0
  ! anti-nu_mu + e+
  REAL(DP), PARAMETER :: CV_ame_p = -0.5d0 + 2.0d0*sin2W
  REAL(DP), PARAMETER :: CA_ame_p = -0.5d0
  ! nu_tau + e+
  REAL(DP), PARAMETER :: CV_te_p  = -0.5d0 + 2.0d0*sin2W
  REAL(DP), PARAMETER :: CA_te_p  = +0.5d0
  ! anti-nu_tau + e+
  REAL(DP), PARAMETER :: CV_ate_p = -0.5d0 + 2.0d0*sin2W
  REAL(DP), PARAMETER :: CA_ate_p = -0.5d0

  !=============================================================================
  ! Bollig Thesis Table 6.4: Neutrino and antineutrino scattering on muons (mu-)
  !=============================================================================
  ! nu_e + mu-
  REAL(DP), PARAMETER :: CV_em_m  = -0.5d0 + 2.0d0*sin2W
  REAL(DP), PARAMETER :: CA_em_m  = -0.5d0
  ! anti-nu_e + mu-
  REAL(DP), PARAMETER :: CV_aem_m = -0.5d0 + 2.0d0*sin2W
  REAL(DP), PARAMETER :: CA_aem_m = +0.5d0
  ! nu_mu + mu-
  REAL(DP), PARAMETER :: CV_mm_m  = +0.5d0 + 2.0d0*sin2W
  REAL(DP), PARAMETER :: CA_mm_m  = +0.5d0
  ! anti-nu_mu + mu-
  REAL(DP), PARAMETER :: CV_amm_m = +0.5d0 + 2.0d0*sin2W
  REAL(DP), PARAMETER :: CA_amm_m = -0.5d0
  ! nu_tau + mu-
  REAL(DP), PARAMETER :: CV_tm_m  = -0.5d0 + 2.0d0*sin2W
  REAL(DP), PARAMETER :: CA_tm_m  = -0.5d0
  ! anti-nu_tau + mu-
  REAL(DP), PARAMETER :: CV_atm_m = -0.5d0 + 2.0d0*sin2W
  REAL(DP), PARAMETER :: CA_atm_m = +0.5d0

  !=============================================================================
  ! Bollig Thesis Table 6.5: Neutrino and antineutrino scattering on antimuons (mu+)
  !=============================================================================
  ! nu_e + mu+
  REAL(DP), PARAMETER :: CV_em_p  = -0.5d0 + 2.0d0*sin2W
  REAL(DP), PARAMETER :: CA_em_p  = +0.5d0
  ! anti-nu_e + mu+
  REAL(DP), PARAMETER :: CV_aem_p = -0.5d0 + 2.0d0*sin2W
  REAL(DP), PARAMETER :: CA_aem_p = -0.5d0
  ! nu_mu + mu+
  REAL(DP), PARAMETER :: CV_mm_p  = +0.5d0 + 2.0d0*sin2W
  REAL(DP), PARAMETER :: CA_mm_p  = -0.5d0
  ! anti-nu_mu + mu+
  REAL(DP), PARAMETER :: CV_amm_p = +0.5d0 + 2.0d0*sin2W
  REAL(DP), PARAMETER :: CA_amm_p = +0.5d0
  ! nu_tau + mu+
  REAL(DP), PARAMETER :: CV_tm_p  = -0.5d0 + 2.0d0*sin2W
  REAL(DP), PARAMETER :: CA_tm_p  = +0.5d0
  ! anti-nu_tau + mu+
  REAL(DP), PARAMETER :: CV_atm_p = -0.5d0 + 2.0d0*sin2W
  REAL(DP), PARAMETER :: CA_atm_p = -0.5d0
  
  ! Integration parameters
  !--- Default quadrature resolution ---
  REAL(DP), PARAMETER :: tfac           = 100.0d0
  REAL(DP), PARAMETER :: MaxExponent    = 100.0d0
  INTEGER, PARAMETER  :: nE3_default    = 32   ! GL points for E3 integral in opacity
  INTEGER, PARAMETER  :: nTheta_default = 16   ! GL points for costheta in opacity
  INTEGER, PARAMETER  :: nPhi_default   = 16   ! GL points for costheta in opacity

  PUBLIC :: CalculatePhoutPhin
  PUBLIC :: GeneralScatteringKernel         ! R_out and R_in at single (E1,E3,costheta) point
  PUBLIC :: GeneralScatteringOpacity        ! angle-integrated outscattering opacity [cm^-1]
  PUBLIC :: ProcessIndexFromReactionString
  PUBLIC :: ReactionStringFromProcessIndex
  PUBLIC :: GeneralAngleAveragedOutgoingScattering_Simple
  PUBLIC :: InverseMuonDecayEmissivity
  PUBLIC :: GeneralAngleAveragedOutgoingScattering_Moments
  PUBLIC :: SelectLambdaFromProcessIndex
  PUBLIC :: Compute_Is_Integrals
  PUBLIC :: A1_f
  PUBLIC :: B1_f
  PUBLIC :: C1_f
  PUBLIC :: C3_f
  PUBLIC :: FD
  PUBLIC :: CalculateCollinearH0
  PUBLIC :: gaulag, gauleg, LegendrePolynomial

CONTAINS

  SUBROUTINE ProcessIndexFromReactionString(reac, ProcessIndex)

    CHARACTER(*), INTENT(IN)  :: reac
    INTEGER     , INTENT(OUT) :: ProcessIndex

    CHARACTER(:), ALLOCATABLE :: s, lhs, rhs
    CHARACTER(32) :: nu_in, l_in, nu_out, l_out
    INTEGER :: pArrow, ios

    ProcessIndex = -1
    s = NormalizeSpaces(ADJUSTL(TRIM(reac)))

    !------------------------------------------------------------
    ! Special 3->1 IMD channels (Table I, 5a & 5b)
    !------------------------------------------------------------
    IF ( s == 'nu_bar_e + e- + nu_mu -> mu-' ) THEN
      ProcessIndex = 33
      RETURN
    ELSE IF ( s == 'nu_mu + e- + nu_bar_e -> mu-' ) THEN
      ProcessIndex = 34
      RETURN
    END IF

    pArrow = INDEX(s, ' -> ')
    IF (pArrow <= 0) THEN
      WRITE(*,*) 'ProcessIndexFromReactionString: missing "->" in: ', TRIM(reac)
      STOP
    END IF

    lhs = TRIM(s(1:pArrow-1))
    rhs = TRIM(s(pArrow+4:LEN(s)))

    CALL ParseSide(lhs, nu_in, l_in, ios)
    IF (ios /= 0) THEN
      WRITE(*,*) 'ProcessIndexFromReactionString: bad LHS: ', TRIM(reac)
      STOP
    END IF

    CALL ParseSide(rhs, nu_out, l_out, ios)
    IF (ios /= 0) THEN
      WRITE(*,*) 'ProcessIndexFromReactionString: bad RHS: ', TRIM(reac)
      STOP
    END IF

    !------------------------------------------------------------
    ! First: Explicit LFE and LFC mappings (25..32)
    !------------------------------------------------------------
    ! LFE: nu_mu + e- <-> nu_e + mu-
    IF (TRIM(nu_in)=='nu_mu' .AND. TRIM(l_in)=='e-' .AND. TRIM(nu_out)=='nu_e' .AND. TRIM(l_out)=='mu-') THEN
      ProcessIndex = 25; RETURN
    ELSE IF (TRIM(nu_in)=='nu_e' .AND. TRIM(l_in)=='mu-' .AND. TRIM(nu_out)=='nu_mu' .AND. TRIM(l_out)=='e-') THEN
      ProcessIndex = 26; RETURN
    END IF

    ! LFE: nu_bar_mu + e+ <-> nu_bar_e + mu+
    IF (TRIM(nu_in)=='nu_bar_mu' .AND. TRIM(l_in)=='e+' .AND. TRIM(nu_out)=='nu_bar_e' .AND. TRIM(l_out)=='mu+') THEN
      ProcessIndex = 27; RETURN
    ELSE IF (TRIM(nu_in)=='nu_bar_e' .AND. TRIM(l_in)=='mu+' .AND. TRIM(nu_out)=='nu_bar_mu' .AND. TRIM(l_out)=='e+') THEN
      ProcessIndex = 28; RETURN
    END IF

    ! LFC: nu_mu + mu+ <-> nu_e + e+
    IF (TRIM(nu_in)=='nu_mu' .AND. TRIM(l_in)=='mu+' .AND. TRIM(nu_out)=='nu_e' .AND. TRIM(l_out)=='e+') THEN
      ProcessIndex = 29; RETURN
    ELSE IF (TRIM(nu_in)=='nu_e' .AND. TRIM(l_in)=='e+' .AND. TRIM(nu_out)=='nu_mu' .AND. TRIM(l_out)=='mu+') THEN
      ProcessIndex = 30; RETURN
    END IF

    ! LFC: nu_bar_mu + mu- <-> nu_bar_e + e-
    IF (TRIM(nu_in)=='nu_bar_mu' .AND. TRIM(l_in)=='mu-' .AND. TRIM(nu_out)=='nu_bar_e' .AND. TRIM(l_out)=='e-') THEN
      ProcessIndex = 31; RETURN
    ELSE IF (TRIM(nu_in)=='nu_bar_e' .AND. TRIM(l_in)=='e-' .AND. TRIM(nu_out)=='nu_bar_mu' .AND. TRIM(l_out)=='mu-') THEN
      ProcessIndex = 32; RETURN
    END IF

    !------------------------------------------------------------
    ! Then: elastic-only map (1..24)
    !------------------------------------------------------------
    IF (TRIM(nu_in) /= TRIM(nu_out) .OR. TRIM(l_in) /= TRIM(l_out)) THEN
      WRITE(*,*) 'ProcessIndexFromReactionString: non-elastic reaction not mapped: ', TRIM(reac)
      STOP
    END IF

    !------------------------------------------------------------
    ! Then: elastic-only map (1..24)
    !------------------------------------------------------------
    IF (TRIM(nu_in) /= TRIM(nu_out) .OR. TRIM(l_in) /= TRIM(l_out)) THEN
      WRITE(*,*) 'ProcessIndexFromReactionString: non-elastic reaction not mapped: ', TRIM(reac)
      WRITE(*,*) 'Supported non-elastic 2->2: Guo (a)-(d). Process 29 is 3->1 (not parsed here).'
      STOP
    END IF

    SELECT CASE (TRIM(nu_in))
    CASE ('nu_e')
      SELECT CASE (TRIM(l_in))
      CASE ('e-');  ProcessIndex = 1
      CASE ('e+');  ProcessIndex = 2
      CASE ('mu-'); ProcessIndex = 13
      CASE ('mu+'); ProcessIndex = 14
      END SELECT

    CASE ('nu_bar_e')
      SELECT CASE (TRIM(l_in))
      CASE ('e-');  ProcessIndex = 3
      CASE ('e+');  ProcessIndex = 4
      CASE ('mu-'); ProcessIndex = 15
      CASE ('mu+'); ProcessIndex = 16
      END SELECT

    CASE ('nu_mu')
      SELECT CASE (TRIM(l_in))
      CASE ('e-');  ProcessIndex = 5
      CASE ('e+');  ProcessIndex = 6
      CASE ('mu-'); ProcessIndex = 17
      CASE ('mu+'); ProcessIndex = 18
      END SELECT

    CASE ('nu_bar_mu')
      SELECT CASE (TRIM(l_in))
      CASE ('e-');  ProcessIndex = 7
      CASE ('e+');  ProcessIndex = 8
      CASE ('mu-'); ProcessIndex = 19
      CASE ('mu+'); ProcessIndex = 20
      END SELECT

    CASE ('nu_tau')
      SELECT CASE (TRIM(l_in))
      CASE ('e-');  ProcessIndex = 9
      CASE ('e+');  ProcessIndex = 10
      CASE ('mu-'); ProcessIndex = 21
      CASE ('mu+'); ProcessIndex = 22
      END SELECT

    CASE ('nu_bar_tau')
      SELECT CASE (TRIM(l_in))
      CASE ('e-');  ProcessIndex = 11
      CASE ('e+');  ProcessIndex = 12
      CASE ('mu-'); ProcessIndex = 23
      CASE ('mu+'); ProcessIndex = 24
      END SELECT
    END SELECT

    IF (ProcessIndex < 0) THEN
      WRITE(*,*) 'ProcessIndexFromReactionString: unsupported reaction: ', TRIM(reac)
      STOP
    END IF

  CONTAINS

    SUBROUTINE ParseSide(side, nu, lep, ios)
      CHARACTER(*), INTENT(IN)  :: side
      CHARACTER(*), INTENT(OUT) :: nu, lep
      INTEGER,      INTENT(OUT) :: ios
      CHARACTER(:), ALLOCATABLE :: t
      INTEGER :: pPlus

      ios = 0
      nu  = ''
      lep = ''

      t = NormalizeSpaces(ADJUSTL(TRIM(side)))
      pPlus = INDEX(t, ' + ')
      IF (pPlus <= 0) THEN
        ios = 1
        RETURN
      END IF

      nu  = ADJUSTL(TRIM(t(1:pPlus-1)))
      lep = ADJUSTL(TRIM(t(pPlus+3:LEN(t))))

      IF (.NOT. IsValidNu(TRIM(nu)))   ios = 2
      IF (.NOT. IsValidLep(TRIM(lep))) ios = 3
    END SUBROUTINE ParseSide

    LOGICAL FUNCTION IsValidNu(x)
      CHARACTER(*), INTENT(IN) :: x
      IsValidNu = (x == 'nu_e'      .OR. x == 'nu_bar_e'  .OR. &
                  x == 'nu_mu'     .OR. x == 'nu_bar_mu' .OR. &
                  x == 'nu_tau'    .OR. x == 'nu_bar_tau')
    END FUNCTION IsValidNu

    LOGICAL FUNCTION IsValidLep(x)
      CHARACTER(*), INTENT(IN) :: x
      IsValidLep = (x == 'e-' .OR. x == 'e+' .OR. x == 'mu-' .OR. x == 'mu+')
    END FUNCTION IsValidLep

    FUNCTION NormalizeSpaces(str) RESULT(out)
      CHARACTER(*), INTENT(IN) :: str
      CHARACTER(:), ALLOCATABLE :: out
      CHARACTER(LEN(str)) :: tmp
      INTEGER :: i, j
      LOGICAL :: prevSpace

      j = 0
      prevSpace = .FALSE.
      DO i = 1, LEN_TRIM(str)
        IF (str(i:i) == ' ' .OR. str(i:i) == CHAR(9)) THEN
          IF (.NOT. prevSpace) THEN
            j = j + 1
            tmp(j:j) = ' '
            prevSpace = .TRUE.
          END IF
        ELSE
          j = j + 1
          tmp(j:j) = str(i:i)
          prevSpace = .FALSE.
        END IF
      END DO

      IF (j == 0) THEN
        out = ''
      ELSE
        out = TRIM(tmp(1:j))
      END IF
    END FUNCTION NormalizeSpaces

  END SUBROUTINE ProcessIndexFromReactionString

  SUBROUTINE ReactionStringFromProcessIndex(ProcessIndex, reac)

    INTEGER,      INTENT(IN)  :: ProcessIndex
    CHARACTER(*), INTENT(OUT) :: reac

    CHARACTER(80) :: s

    SELECT CASE (ProcessIndex)
    CASE (1);  s = 'nu_e       + e-  -> nu_e       + e- '
    CASE (2);  s = 'nu_e       + e+  -> nu_e       + e+ '
    CASE (3);  s = 'nu_bar_e   + e-  -> nu_bar_e   + e- '
    CASE (4);  s = 'nu_bar_e   + e+  -> nu_bar_e   + e+ '

    CASE (5);  s = 'nu_mu      + e-  -> nu_mu      + e- '
    CASE (6);  s = 'nu_mu      + e+  -> nu_mu      + e+ '
    CASE (7);  s = 'nu_bar_mu  + e-  -> nu_bar_mu  + e- '
    CASE (8);  s = 'nu_bar_mu  + e+  -> nu_bar_mu  + e+ '

    CASE (9);  s = 'nu_tau     + e-  -> nu_tau     + e- '
    CASE (10); s = 'nu_tau     + e+  -> nu_tau     + e+ '
    CASE (11); s = 'nu_bar_tau + e-  -> nu_bar_tau + e- '
    CASE (12); s = 'nu_bar_tau + e+  -> nu_bar_tau + e+ '

    CASE (13); s = 'nu_e       + mu- -> nu_e       + mu-'
    CASE (14); s = 'nu_e       + mu+ -> nu_e       + mu+'
    CASE (15); s = 'nu_bar_e   + mu- -> nu_bar_e   + mu-'
    CASE (16); s = 'nu_bar_e   + mu+ -> nu_bar_e   + mu+'

    CASE (17); s = 'nu_mu      + mu- -> nu_mu      + mu-'
    CASE (18); s = 'nu_mu      + mu+ -> nu_mu      + mu+'
    CASE (19); s = 'nu_bar_mu  + mu- -> nu_bar_mu  + mu-'
    CASE (20); s = 'nu_bar_mu  + mu+ -> nu_bar_mu  + mu+'

    CASE (21); s = 'nu_tau     + mu- -> nu_tau     + mu-'
    CASE (22); s = 'nu_tau     + mu+ -> nu_tau     + mu+'
    CASE (23); s = 'nu_bar_tau + mu- -> nu_bar_tau + mu-'
    CASE (24); s = 'nu_bar_tau + mu+ -> nu_bar_tau + mu+'

  ! --- LFE: Lepton Flavor Exchange (Table I, 3a & 3b) ---
    CASE (25); s = 'nu_mu      + e-  -> nu_e       + mu- '  ! 3a forward
    CASE (26); s = 'nu_e       + mu- -> nu_mu      + e-  '  ! 3a reverse
    CASE (27); s = 'nu_bar_mu  + e+  -> nu_bar_e   + mu+ '  ! 3b forward
    CASE (28); s = 'nu_bar_e   + mu+ -> nu_bar_mu  + e+  '  ! 3b reverse

    ! --- LFC: Lepton Flavor Conversion (Table I, 4a & 4b) ---
    CASE (29); s = 'nu_mu      + mu+ -> nu_e       + e+  '  ! 4a forward
    CASE (30); s = 'nu_e       + e+  -> nu_mu      + mu+ '  ! 4a reverse
    CASE (31); s = 'nu_bar_mu  + mu- -> nu_bar_e   + e-  '  ! 4b forward
    CASE (32); s = 'nu_bar_e   + e-  -> nu_bar_mu  + mu- '  ! 4b reverse

    ! --- IMD: Inverse Muon Decay (IMD) (Table I, 5a & 5b) ---
    CASE (33); s = 'nu_bar_e   + e- + nu_mu -> mu-  '       ! inverse muon decay (3->1)
    CASE (34); s = 'nu_mu      + e- + nu_bar_e -> mu-  '    ! inverse muon decay (3->1)

    CASE DEFAULT
      WRITE(*,*) 'ReactionStringFromProcessIndex: unsupported ProcessIndex = ', ProcessIndex
      STOP
    END SELECT

    reac = TRIM(s)
  END SUBROUTINE ReactionStringFromProcessIndex

  SUBROUTINE SelectLambdaFromProcessIndex(ProcessIndex, lam1, lam2, lam3)

    INTEGER , INTENT(IN)  :: ProcessIndex
    REAL(DP), INTENT(OUT) :: lam1, lam2, lam3

    SELECT CASE(ProcessIndex)
    ! ============================================================
    ! ELASTIC NC SCATTERING ONLY (same in/out flavor): 1..24
    ! ============================================================

    ! 1) nu_e + e-
    CASE(1);  CALL CVCA_to_lambda(CV_ee_m,  CA_ee_m, lam1, lam2, lam3)
    ! 2) nu_e + e+
    CASE(2);  CALL CVCA_to_lambda(CV_ee_p,  CA_ee_p, lam1, lam2, lam3)
    ! 3) nu_bar_e + e-
    CASE(3);  CALL CVCA_to_lambda(CV_aee_m, CA_aee_m, lam1, lam2, lam3)
    ! 4) nu_bar_e + e+
    CASE(4);  CALL CVCA_to_lambda(CV_aee_p, CA_aee_p, lam1, lam2, lam3)

    ! 5) nu_mu + e-
    CASE(5);  CALL CVCA_to_lambda(CV_me_m,  CA_me_m, lam1, lam2, lam3)
    ! 6) nu_mu + e+
    CASE(6);  CALL CVCA_to_lambda(CV_me_p,  CA_me_p, lam1, lam2, lam3)
    ! 7) nu_bar_mu + e-
    CASE(7);  CALL CVCA_to_lambda(CV_ame_m, CA_ame_m, lam1, lam2, lam3)
    ! 8) nu_bar_mu + e+
    CASE(8);  CALL CVCA_to_lambda(CV_ame_p, CA_ame_p, lam1, lam2, lam3)

    ! 9) nu_tau + e-
    CASE(9);  CALL CVCA_to_lambda(CV_te_m,  CA_te_m, lam1, lam2, lam3)
    ! 10) nu_tau + e+
    CASE(10); CALL CVCA_to_lambda(CV_te_p,  CA_te_p, lam1, lam2, lam3)
    ! 11) nu_bar_tau + e-
    CASE(11); CALL CVCA_to_lambda(CV_ate_m, CA_ate_m, lam1, lam2, lam3)
    ! 12) nu_bar_tau + e+
    CASE(12); CALL CVCA_to_lambda(CV_ate_p, CA_ate_p, lam1, lam2, lam3)

    ! 13) nu_e + mu-
    CASE(13); CALL CVCA_to_lambda(CV_em_m,  CA_em_m, lam1, lam2, lam3)
    ! 14) nu_e + mu+
    CASE(14); CALL CVCA_to_lambda(CV_em_p,  CA_em_p, lam1, lam2, lam3)
    ! 15) nu_bar_e + mu-
    CASE(15); CALL CVCA_to_lambda(CV_aem_m, CA_aem_m, lam1, lam2, lam3)
    ! 16) nu_bar_e + mu+
    CASE(16); CALL CVCA_to_lambda(CV_aem_p, CA_aem_p, lam1, lam2, lam3)

    ! 17) nu_mu + mu-
    CASE(17); CALL CVCA_to_lambda(CV_mm_m,  CA_mm_m, lam1, lam2, lam3)
    ! 18) nu_mu + mu+
    CASE(18); CALL CVCA_to_lambda(CV_mm_p,  CA_mm_p, lam1, lam2, lam3)
    ! 19) nu_bar_mu + mu-
    CASE(19); CALL CVCA_to_lambda(CV_amm_m, CA_amm_m, lam1, lam2, lam3)
    ! 20) nu_bar_mu + mu+
    CASE(20); CALL CVCA_to_lambda(CV_amm_p, CA_amm_p, lam1, lam2, lam3)

    ! 21) nu_tau + mu-
    CASE(21); CALL CVCA_to_lambda(CV_tm_m,  CA_tm_m, lam1, lam2, lam3)
    ! 22) nu_tau + mu+
    CASE(22); CALL CVCA_to_lambda(CV_tm_p,  CA_tm_p, lam1, lam2, lam3)
    ! 23) nu_bar_tau + mu-
    CASE(23); CALL CVCA_to_lambda(CV_atm_m, CA_atm_m, lam1, lam2, lam3)
    ! 24) nu_bar_tau + mu+
    CASE(24); CALL CVCA_to_lambda(CV_atm_p, CA_atm_p, lam1, lam2, lam3)

    ! ============================================================
    ! INELASTIC FLAVOR EXCHANGE (LFE) & CONVERSION (LFC): 25..32
    ! ============================================================
    
    ! LFE (t-channel W-exchange) 
    ! nu_mu + e- <-> nu_e + mu- AND nu_bar_mu + e+ <-> nu_bar_e + mu+
    CASE(25, 26, 27, 28)
      lam1 = lam1_LFE
      lam2 = lam2_LFE
      lam3 = lam3_LFE

    ! LFC (s-channel W-exchange)
    ! nu_mu + mu+ <-> nu_e + e+ AND nu_bar_mu + mu- <-> nu_bar_e + e-
    CASE(29, 30, 31, 32)
      lam1 = lam1_LFC
      lam2 = lam2_LFC
      lam3 = lam3_LFC

    ! IMD: Inverse Muon Decay (IMD) (Table I, 5a & 5b)
    CASE(33)  ! inverse muon decay (3->1)
      lam1 = lam1_LFE ! LFE and IMD have the same lambdas
      lam2 = lam2_LFE ! LFE and IMD have the same lambdas
      lam3 = lam3_LFE ! LFE and IMD have the same lambdas

    CASE(34)  ! inverse muon decay (3->1)
      lam1 = lam1_LFC ! Same as above but now you switched places for neutrinos
      lam2 = lam2_LFC ! Same as above but now you switched places for neutrinos
      lam3 = lam3_LFC ! Same as above but now you switched places for neutrinos

    CASE DEFAULT
      WRITE(*,*) 'Error: Unrecognized ProcessIndex: ', ProcessIndex
      STOP
    END SELECT

  END SUBROUTINE SelectLambdaFromProcessIndex

  !================================================================================
  !
  !  GeneralScatteringKernel
  !
  !  Compute outscattering kernel R_out and inscattering kernel R_in at a single
  !  phase-space point (E1, E3, costheta) for the general Scattering process:
  !    nu_l2(E1) + l2(E2)  ->  nu_l4(E3) + l4(E4)
  !
  !  Transport convention (Guo Eq. 2 and Fischer Eqs. B17-B19):
  !
  !   R_out(E1,E3,costheta): kernel for the FORWARD process
  !     nu_l2(E1) + l2(E2) -> nu_l4(E3) + l4(E4)
  !   Occupation factors: f2(E2) * [1 - f4(E4)]
  !   Appears in collision integral reducing F_{nu_l2}(E1).
  !
  !   R_in(E1,E3,costheta): kernel for the REVERSE process
  !     nu_l4(E3) + l4(E4) -> nu_l2(E1) + l2(E2)
  !   Occupation factors: f4(E4) * [1 - f2(E2)]
  !   Appears in collision integral increasing F_{nu_l2}(E1).
  !
  !   Detailed balance (Guo Eq. 13 / Fischer Eq. B29):
  !     R_in(E1,E3) = R_out(E1,E3) * exp((E1-E3+xMu_l2-xMu_l4)/T)
  !
  !  Inputs:
  !    E1        -- incoming nu_l2 energy [MeV]
  !    E3        -- outgoing nu_l4 energy [MeV]
  !    costheta  -- cosine of angle between p_{nu_l2} and p_{nu_l4}
  !    xT        -- temperature [MeV]
  !    xMu_l2      -- l2 relativistic chemical potential [MeV]
  !    xMu_l4      -- l4 relativistic chemical potential [MeV]
  !    m2        -- l2 mass [MeV]
  !    m4        -- l4 mass [MeV]
  !    ProcessIndex -- 1: particle channel (nu_l2 + l2 -> nu_l4 + l4)
  !                     2: antiparticle channel (nu_bar_l2 + l2bar -> nu_bar_l4 + l4bar)
  !                        implemented by negating both chemical potentials
  !
  !  Outputs:
  !    Notice that for some reason Guo 2020 swaps the definitions of R_out and R_in 
  !    compared to Fischer 2020 and Bollig's thesis, so we DO NOT follow Guo's convention here:
  !    Rout  -- R_out(E1, E3, costheta) [MeV^-4]
  !    Rin   -- R_in (E1, E3, costheta) [MeV^-4]
  !
  !================================================================================
  SUBROUTINE CalculatePhoutPhin( E1, E3, xT, xMu_l2, xMu_l4, m2, m4, &
                                ProcessIndex, Phout, Phin, nL, nTheta_in )

    REAL(DP), INTENT(IN)  :: E1, E3, xT, xMu_l2, xMu_l4, m2, m4
    INTEGER , INTENT(IN)  :: ProcessIndex
    REAL(DP), INTENT(OUT) :: Phout(nL), Phin(nL)
    INTEGER , INTENT(IN)  :: nL
    INTEGER , INTENT(IN), OPTIONAL :: nTheta_in

    REAL(DP) :: Rout, Rin, Pl_mu
    REAL(DP) :: costheta, Delta_mu, exponent, conv_fac
    INTEGER  :: iTh, nTheta, iL
    REAL(DP), ALLOCATABLE :: xa_cos_theta(:), wa_cos_theta(:)

    nTheta = nTheta_default
    IF (PRESENT(nTheta_in)) nTheta = nTheta_in

    ALLOCATE( xa_cos_theta(nTheta), wa_cos_theta(nTheta) )

    CALL gauleg( -1.0d0, 1.0d0, xa_cos_theta, wa_cos_theta, nTheta )
    
    conv_fac = 1.0d0 / ( (2.0d0 * pi)**3 * hbarc ) 

    Phout = 0.0d0
    Phin = 0.0d0
    DO iTh=1, nTheta
      costheta = xa_cos_theta(iTh)
      CALL GeneralScatteringKernel( E1, E3, costheta, xT, xMu_l2, xMu_l4, m2, m4, &
                                ProcessIndex, Rout, Rin )

      ! --- Accumulate the integral for each moment ---
      DO iL = 1, nL
        Pl_mu = LegendrePolynomial(iL-1, costheta)
        Phout(iL) = Phout(iL) + Rout * Pl_mu * wa_cos_theta(iTh)
      END DO
      
    ENDDO
    ! factors of 1/2 and 3/2 to match Bruenn
    Phout(1) = Phout(1) * conv_fac * 0.5d0
    Phout(2) = Phout(2) * conv_fac * 1.5d0

    ! --- Phin from detailed balance 
    Delta_mu = xMu_l2 - xMu_l4
    exponent = MIN( (E3 - E1 - Delta_mu) / xT, MaxExponent )
    Phin = Phout * EXP(exponent)

  END SUBROUTINE CalculatePhoutPhin

  SUBROUTINE GeneralScatteringKernel( E1, E3, costheta, xT, xMu_l2, xMu_l4, m2, m4, &
                                ProcessIndex, Rout, Rin )

    REAL(DP), INTENT(IN)  :: E1, E3, costheta, xT, xMu_l2, xMu_l4, m2, m4
    INTEGER , INTENT(IN)  :: ProcessIndex
    REAL(DP), INTENT(OUT) :: Rout, Rin

    REAL(DP) :: exponent, Delta_mu
    REAL(DP) :: lam1, lam2, lam3

    CALL SelectLambdaFromProcessIndex(ProcessIndex, lam1, lam2, lam3)

    ! --- compute R_out with f2*(1-f4) factors ---
    CALL GeneralScatteringRout( E1, E3, costheta, xT, xMu_l2, xMu_l4, m2, m4, lam1, lam2, lam3, Rout )

    ! --- R_in from detailed balance 
    Delta_mu = xMu_l2 - xMu_l4
    exponent = MIN( (E3 - E1 - Delta_mu) / xT, MaxExponent )   ! Guo Eq. 13

    Rin = Rout * EXP(exponent)   ! Guo Eq. 13

  END SUBROUTINE GeneralScatteringKernel

  SUBROUTINE InverseMuonDecayOpacity( E1, xT, xMu_l2, xMu_l4, xMu_nu3, m2, m4, &
                                      ProcessIndex, Opacity, nE3, nTheta )

    ! --- Dummy Arguments ---
    REAL(DP), INTENT(IN)  :: E1, xT, xMu_l2, xMu_l4, xMu_nu3, m2, m4
    INTEGER , INTENT(IN)  :: ProcessIndex
    REAL(DP), INTENT(OUT) :: Opacity
    INTEGER , INTENT(IN)  :: nE3, nTheta

    ! --- Local Variables ---
    INTEGER  :: iE3, iTheta
    REAL(DP) :: E3_min, E3_max, E3_range
    REAL(DP) :: E3, wE3, cos_theta, wcos_theta
    REAL(DP) :: Rout, bin_opacity, f_nu3
    REAL(DP), ALLOCATABLE :: xa_E3(:), wa_E3(:)
    REAL(DP), ALLOCATABLE :: xa_cos_theta(:), wa_cos_theta(:)
    REAL(DP) :: conv_fac, lam1, lam2, lam3

    IF (ProcessIndex /= 33) THEN
      WRITE(*,*) 'Error: InverseMuonDecayOpacity only valid for ProcessIndex=33 (nu_bar_e + e- + nu_mu -> mu-).'
      STOP
    END IF
    CALL SelectLambdaFromProcessIndex(ProcessIndex, lam1, lam2, lam3)

    ! --- E3 integration limits ---
    E3_min = 0.0d0
    E3_max = E1 + ABS(xMu_l2) + tfac*xT
    E3_range = E3_max - E3_min

    ! --- Allocations & Quadrature ---
    ALLOCATE( xa_E3(nE3), wa_E3(nE3) )
    ALLOCATE( xa_cos_theta(nTheta), wa_cos_theta(nTheta) )
    
    CALL gauleg(  E3_min, E3_max,  xa_E3,    wa_E3,    nE3  )
    CALL gauleg( -1.0d0, 1.0d0, xa_cos_theta, wa_cos_theta, nTheta )

    Opacity = 0.0d0
    
    conv_fac = 1.0d0 / ( (2.0d0 * pi)**3 * hbarc ) 

    ! --- Integration Loops ---
    DO iE3 = 1, nE3
      E3    = xa_E3(iE3)
      wE3   = wa_E3(iE3)
      f_nu3 = FD(E3, xMu_nu3, xT)

      ! Replaced the triple loop with a single loop over relative angle mu
      DO iTheta = 1, nTheta
        cos_theta  = xa_cos_theta(iTheta)
        wcos_theta = wa_cos_theta(iTheta)

        CALL InverseMuonDecayRout( E1, E3, cos_theta, xT, xMu_l2, xMu_l4, m2, m4, lam1, lam2, lam3, Rout )
        Opacity = Opacity + Rout * wcos_theta * E3**2 * wE3 * conv_fac * f_nu3  ! include nu3 occupation factor
      END DO
    ENDDO

    DEALLOCATE( xa_E3, wa_E3, xa_cos_theta, wa_cos_theta )

  END SUBROUTINE InverseMuonDecayOpacity

  SUBROUTINE InverseMuonDecayEmissivity( E1, xT, xMu_l2, xMu_l4, xMu_nu3, m2, m4, &
                    ProcessIndex, Emissivity, nE3, nTheta, IsFinalStateFree_in )

    ! --- Dummy Arguments ---
    REAL(DP), INTENT(IN)  :: E1, xT, xMu_l2, xMu_l4, xMu_nu3, m2, m4
    INTEGER , INTENT(IN)  :: ProcessIndex
    REAL(DP), INTENT(OUT) :: Emissivity
    INTEGER , INTENT(IN)  :: nE3, nTheta
    LOGICAL , INTENT(IN), OPTIONAL :: IsFinalStateFree_in

    ! --- Local Variables ---
    INTEGER  :: iE3, iTheta
    REAL(DP) :: E3_min, E3_max, E3_range
    REAL(DP) :: E3, wE3, cos_theta, wcos_theta
    REAL(DP) :: Rout, f_nu3
    REAL(DP), ALLOCATABLE :: xa_E3(:), wa_E3(:)
    REAL(DP), ALLOCATABLE :: xa_cos_theta(:), wa_cos_theta(:)
    REAL(DP) :: conv_fac, lam1, lam2, lam3
    LOGICAL  :: IsFinalStateFree

    IsFinalStateFree = .FALSE.
    IF (PRESENT(IsFinalStateFree_in)) IsFinalStateFree = IsFinalStateFree_in

    IF (ProcessIndex /= 33) THEN
      WRITE(*,*) 'Error: InverseMuonDecayEmissivity only valid for ProcessIndex=33 (nu_bar_e + e- + nu_mu -> mu-).'
      STOP
    END IF
    CALL SelectLambdaFromProcessIndex(ProcessIndex, lam1, lam2, lam3)

    ! --- E3 integration limits ---
    E3_min = 0.0d0
    E3_max = E1 + ABS(xMu_l2) + tfac*xT
    E3_range = E3_max - E3_min

    ! --- Allocations & Quadrature ---
    ALLOCATE( xa_E3(nE3), wa_E3(nE3) )
    ALLOCATE( xa_cos_theta(nTheta), wa_cos_theta(nTheta) )
    
    CALL gauleg(  E3_min, E3_max,  xa_E3,    wa_E3,    nE3  )
    CALL gauleg( -1.0d0, 1.0d0, xa_cos_theta, wa_cos_theta, nTheta )

    Emissivity = 0.0d0
    
    conv_fac = 1.0d0 / ( (2.0d0 * pi)**3 * hbarc ) 

    ! --- Integration Loops ---
    DO iE3 = 1, nE3
      E3    = xa_E3(iE3)
      wE3   = wa_E3(iE3)
      IF (IsFinalStateFree) THEN
        f_nu3 = 0.0d0
      ELSE
        f_nu3 = FD(E3, xMu_nu3, xT)
      ENDIF
      ! Replaced the triple loop with a single loop over relative angle mu
      DO iTheta = 1, nTheta
        cos_theta  = xa_cos_theta(iTheta)
        wcos_theta = wa_cos_theta(iTheta)

        CALL InverseMuonDecayRout( E1, E3, cos_theta, xT, xMu_l2, xMu_l4, m2, m4, lam1, lam2, lam3, Rout )
        Emissivity = Emissivity + Rout * wcos_theta * E3**2 * wE3 * conv_fac * (1.0d0 - f_nu3)  ! include nu3 occupation factor
      END DO
    ENDDO

    DEALLOCATE( xa_E3, wa_E3, xa_cos_theta, wa_cos_theta )

  END SUBROUTINE InverseMuonDecayEmissivity

  SUBROUTINE GeneralScatteringOpacity( E1, xT, xMu_l2, xMu_l4, xMu_nu3, m2, m4, &
                                      ProcessIndex, Opacity_tot, &
                                      nE3, nTheta, nPhi_in, &
                                      IsFinalStateFree_in, DoFullIntegration_in )

    ! --- Dummy Arguments ---
    REAL(DP), INTENT(IN)  :: E1, xT, xMu_l2, xMu_l4, xMu_nu3, m2, m4
    INTEGER , INTENT(IN)  :: ProcessIndex
    REAL(DP), INTENT(OUT) :: Opacity_tot
    INTEGER , INTENT(IN)  :: nE3, nTheta
    INTEGER , INTENT(IN), OPTIONAL :: nPhi_in
    LOGICAL , INTENT(IN), OPTIONAL :: IsFinalStateFree_in, DoFullIntegration_in

    ! --- Local Variables ---
    LOGICAL :: IsFinalStateFree, DoFullIntegration
    REAL(DP), DIMENSION(3,nE3) :: Opacity_array

    IsFinalStateFree = .FALSE.
    IF (PRESENT(IsFinalStateFree_in)) IsFinalStateFree = IsFinalStateFree_in
    DoFullIntegration = .FALSE.
    IF (PRESENT(DoFullIntegration_in)) DoFullIntegration = DoFullIntegration_in

    IF (ProcessIndex == 33) THEN
      CALL InverseMuonDecayOpacity( E1, xT, xMu_l2, xMu_l4, xMu_nu3, m2, m4, &
                                   ProcessIndex, Opacity_tot, nE3, nTheta )
      RETURN
    END IF

    IF (DoFullIntegration) THEN
      CALL GeneralAngleAveragedOutgoingScattering_Full( &
                    E1, xT, xMu_l2, xMu_l4, xMu_nu3, m2, m4, &
                    ProcessIndex, Opacity_array, &
                    nE3, nTheta, nPhi_in, IsFinalStateFree )
      ! Only valid if IsFinalStateFree_in
    ELSE
      CALL GeneralAngleAveragedOutgoingScattering_Simple( &
                    E1, xT, xMu_l2, xMu_l4, xMu_nu3, m2, m4, &
                    ProcessIndex, Opacity_array, &
                    nE3, nTheta, IsFinalStateFree )
    ENDIF

    Opacity_tot = SUM(Opacity_array(2,:) * Opacity_array(3,:))

  END SUBROUTINE GeneralScatteringOpacity

  SUBROUTINE GeneralAngleAveragedOutgoingScattering_Full( E1, xT, xMu_l2, xMu_l4, xMu_nu3, m2, m4, &
                                                    ProcessIndex, Opacity_array, &
                                                    nE3, nTheta, nPhi_in, IsFinalStateFree_in )

    ! --- Dummy Arguments ---
    REAL(DP), INTENT(IN)  :: E1, xT, xMu_l2, xMu_l4, xMu_nu3, m2, m4
    INTEGER , INTENT(IN)  :: ProcessIndex
    REAL(DP), DIMENSION(3,nE3), INTENT(OUT) :: Opacity_array(:,:)
    INTEGER , INTENT(IN) :: nE3, nTheta
    INTEGER , INTENT(IN), OPTIONAL :: nPhi_in
    LOGICAL , INTENT(IN), OPTIONAL :: IsFinalStateFree_in

    ! --- Local Variables ---
    LOGICAL  :: IsFinalStateFree
    INTEGER  :: nPhi, iE3, iTheta, iTheta_p, iPhi
    REAL(DP) :: E3_min, E3_max, E3_range
    REAL(DP) :: E3, wE3, cos_theta, wcos_theta, cos_theta_p, wcos_theta_p, Phi, wPhi
    REAL(DP) :: Rout, bin_opacity, f_nu3
    REAL(DP), ALLOCATABLE :: xa_E3(:), wa_E3(:)
    REAL(DP), ALLOCATABLE :: xa_cos_theta(:), wa_cos_theta(:)
    REAL(DP), ALLOCATABLE :: xa_Phi(:), wa_Phi(:)
    REAL(DP) :: conv_fac, lam1, lam2, lam3, mu

    nPhi   = nPhi_default
    IsFinalStateFree = .FALSE.
    IF (PRESENT(nPhi_in))             nPhi             = nPhi_in
    IF (PRESENT(IsFinalStateFree_in)) IsFinalStateFree = IsFinalStateFree_in

    CALL SelectLambdaFromProcessIndex(ProcessIndex, lam1, lam2, lam3)
    
    ! --- E3 integration limits ---
    E3_min = 0.0d0
    E3_max = E1 + ABS(xMu_l2) + tfac*xT
    E3_range = E3_max - E3_min

    ! --- Allocations & Quadrature ---
    ALLOCATE( xa_E3(nE3)    , wa_E3(nE3)     )
    ALLOCATE( xa_cos_theta(nTheta), wa_cos_theta(nTheta) )
    ALLOCATE( xa_Phi(nPhi), wa_Phi(nPhi) )
    
    CALL gauleg(  E3_min, E3_max,  xa_E3,    wa_E3,    nE3  )
    CALL gauleg( -1.0d0 , 1.0d0 ,  xa_cos_theta, wa_cos_theta, nTheta )
    CALL gauleg(  0.0d0 , 1.0d0 ,  xa_Phi, wa_Phi, nPhi )

    Opacity_array = 0.0d0
    
    ! This is correct but I need to rederive why. Everyone uses different conventions... 
    ! Start with Lohs and Mezzacappa and Bruenn 1993 perhaps
    conv_fac = 1.0d0 / ( (2.0d0 * pi)**3 * hbarc ) 

    ! --- Integration Loops ---
    DO iE3 = 1, nE3
      E3    = xa_E3(iE3)
      wE3   = wa_E3(iE3)
      IF (IsFinalStateFree) THEN
        f_nu3 = 0.0d0
      ELSE
        f_nu3 = FD( E3, xMu_nu3, xT )
      ENDIF
      
      bin_opacity = 0.0d0

      DO iTheta = 1, nTheta
        cos_theta  = xa_cos_theta(iTheta)
        wcos_theta = wa_cos_theta(iTheta)
        DO iTheta_p = 1, nTheta
          cos_theta_p  = xa_cos_theta(iTheta_p)
          wcos_theta_p = wa_cos_theta(iTheta_p)

          DO iPhi = 1, nPhi
            Phi  = xa_Phi(iPhi) * 2.0d0*pi
            wPhi = wa_Phi(iPhi) * 2.0d0*pi
            mu = cos_theta*cos_theta_p + SQRT((1.0d0-cos_theta**2)*(1.0d0-cos_theta_p**2))*COS(Phi)

            CALL GeneralScatteringRout( E1, E3, mu, xT, xMu_l2, xMu_l4, m2, m4, lam1, lam2, lam3, Rout )
            bin_opacity = bin_opacity + Rout * wPhi * wcos_theta_p * wcos_theta   
          ENDDO
        ENDDO
      END DO
      
      ! Multiply by phase space factors: E3^2 * dE3 * final_state_blocking
      Opacity_array(1,iE3) = E3
      Opacity_array(2,iE3) = wE3
      Opacity_array(3,iE3) = bin_opacity * (1.0d0 - f_nu3) * E3**2 * conv_fac
      
    ENDDO

    DEALLOCATE( xa_E3, wa_E3, xa_cos_theta, wa_cos_theta, xa_Phi, wa_Phi )

  END SUBROUTINE GeneralAngleAveragedOutgoingScattering_Full

  SUBROUTINE GeneralAngleAveragedOutgoingScattering_Simple( E1, xT, xMu_l2, xMu_l4, xMu_nu3, m2, m4, &
                                                    ProcessIndex, Opacity_array, &
                                                    nE3, nTheta, IsFinalStateFree_in )

    ! --- Dummy Arguments ---
    REAL(DP), INTENT(IN)  :: E1, xT, xMu_l2, xMu_l4, xMu_nu3, m2, m4
    INTEGER , INTENT(IN)  :: ProcessIndex
    REAL(DP), DIMENSION(3,nE3), INTENT(OUT) :: Opacity_array(:,:)
    INTEGER , INTENT(IN) :: nE3, nTheta
    LOGICAL , INTENT(IN), OPTIONAL :: IsFinalStateFree_in

    ! --- Local Variables ---
    LOGICAL  :: IsFinalStateFree
    INTEGER  :: iE3, iMu
    REAL(DP) :: E3_min, E3_max, E3_range
    REAL(DP) :: E3, wE3, mu, wmu
    REAL(DP) :: Rout, bin_opacity, f_nu3
    REAL(DP), ALLOCATABLE :: xa_E3(:), wa_E3(:)
    REAL(DP), ALLOCATABLE :: xa_mu(:), wa_mu(:)
    REAL(DP) :: conv_fac, lam1, lam2, lam3

    IsFinalStateFree = .FALSE.
    IF (PRESENT(IsFinalStateFree_in)) IsFinalStateFree = IsFinalStateFree_in

    CALL SelectLambdaFromProcessIndex(ProcessIndex, lam1, lam2, lam3)

    ! --- E3 integration limits ---
    E3_min = 0.0d0
    E3_max = E1 + ABS(xMu_l2) + tfac*xT
    E3_range = E3_max - E3_min

    ! --- Allocations & Quadrature ---
    ALLOCATE( xa_E3(nE3), wa_E3(nE3) )
    ALLOCATE( xa_mu(nTheta), wa_mu(nTheta) )
    
    CALL gauleg(  E3_min, E3_max,  xa_E3,    wa_E3,    nE3  )
    CALL gauleg( -1.0d0, 1.0d0, xa_mu, wa_mu, nTheta )

    Opacity_array = 0.0d0
    
    conv_fac = 1.0d0 / ( (2.0d0 * pi)**3 * hbarc ) 

    ! --- Integration Loops ---
    DO iE3 = 1, nE3
      E3    = xa_E3(iE3)
      wE3   = wa_E3(iE3)
      
      IF (IsFinalStateFree) THEN
        f_nu3 = 0.0d0
      ELSE
        f_nu3 = FD( E3, xMu_nu3, xT )
      ENDIF
      bin_opacity = 0.0d0

      ! Replaced the triple loop with a single loop over relative angle mu
      DO iMu = 1, nTheta
        mu  = xa_mu(iMu)
        wmu = wa_mu(iMu)

        CALL GeneralScatteringRout( E1, E3, mu, xT, xMu_l2, xMu_l4, m2, m4, lam1, lam2, lam3, Rout )
        bin_opacity = bin_opacity + Rout * wmu   
      END DO
      
      ! Multiply by phase space factors: E3^2 * dE3 * final_state_blocking
      Opacity_array(1,iE3) = E3
      Opacity_array(2,iE3) = wE3
      Opacity_array(3,iE3) = bin_opacity * E3**2 * conv_fac * (1.0d0 - f_nu3)  ! include final state neutrino blocking factor
      
    ENDDO

    DEALLOCATE( xa_E3, wa_E3, xa_mu, wa_mu )

  END SUBROUTINE GeneralAngleAveragedOutgoingScattering_Simple

  SUBROUTINE GeneralAngleAveragedOutgoingScattering_Moments( E1, xT, xMu_l2, xMu_l4, xMu_nu3, m2, m4, &
                                                              ProcessIndex, Opacity_array, &
                                                              nE3, nTheta, nMoments, IsFinalStateFree_in )

    ! --- Dummy Arguments ---
    REAL(DP), INTENT(IN)  :: E1, xT, xMu_l2, xMu_l4, xMu_nu3, m2, m4
    INTEGER , INTENT(IN)  :: ProcessIndex
    INTEGER , INTENT(IN)  :: nE3, nTheta, nMoments  ! nMoments = number of Legendre moments to compute
    REAL(DP), DIMENSION(2+nMoments, nE3), INTENT(OUT) :: Opacity_array
    LOGICAL , INTENT(IN), OPTIONAL :: IsFinalStateFree_in

    ! --- Local Variables ---
    LOGICAL  :: IsFinalStateFree
    INTEGER  :: iE3, iMu, l
    REAL(DP) :: E3_min, E3_max
    REAL(DP) :: E3, wE3, mu, wmu
    REAL(DP) :: Rout, f_nu3
    REAL(DP), ALLOCATABLE :: xa_E3(:), wa_E3(:)
    REAL(DP), ALLOCATABLE :: xa_mu(:), wa_mu(:)
    REAL(DP), ALLOCATABLE :: Phi_out(:)  ! Legendre moments for current E3
    REAL(DP) :: conv_fac, lam1, lam2, lam3
    REAL(DP) :: Pl_mu  ! Legendre polynomial P_l(mu)

    IsFinalStateFree = .FALSE.
    IF (PRESENT(IsFinalStateFree_in)) IsFinalStateFree = IsFinalStateFree_in

    CALL SelectLambdaFromProcessIndex(ProcessIndex, lam1, lam2, lam3)

    ! --- E3 integration limits ---
    E3_min = 0.0d0
    E3_max = E1 + ABS(xMu_l2) + tfac*xT

    ! --- Allocations & Quadrature ---
    ALLOCATE( xa_E3(nE3), wa_E3(nE3) )
    ALLOCATE( xa_mu(nTheta), wa_mu(nTheta) )
    ALLOCATE( Phi_out(0:nMoments-1) )
    
    CALL gauleg( E3_min, E3_max, xa_E3, wa_E3, nE3 )
    CALL gauleg( -1.0d0, 1.0d0, xa_mu, wa_mu, nTheta )

    Opacity_array = 0.0d0
    
    conv_fac = 1.0d0 / ( (2.0d0 * pi)**3 * hbarc ) 

    ! --- Integration Loops ---
    DO iE3 = 1, nE3
      E3  = xa_E3(iE3)
      wE3 = wa_E3(iE3)
      
      IF (IsFinalStateFree) THEN
        f_nu3 = 0.0d0
      ELSE
        f_nu3 = FD( E3, xMu_nu3, xT )
      ENDIF

      ! Initialize moments for this E3
      Phi_out(:) = 0.0d0

      ! Loop over angular quadrature points
      DO iMu = 1, nTheta
        mu  = xa_mu(iMu)
        wmu = wa_mu(iMu)

        ! Get scattering rate at this angle
        CALL GeneralScatteringRout( E1, E3, mu, xT, xMu_l2, xMu_l4, m2, m4, lam1, lam2, lam3, Rout )
        
        ! Calculate Legendre moments:
        ! Phi_out(l) = ∫ R_out(E1, E3, mu) * P_l(mu) d(mu)
        DO l = 0, nMoments-1
          Pl_mu = LegendrePolynomial(l, mu)
          Phi_out(l) = Phi_out(l) + Rout * Pl_mu * wmu
        END DO
        
      END DO
      
      ! Store results
      ! Opacity_array(1, iE3) = E3
      ! Opacity_array(2, iE3) = wE3
      ! Opacity_array(3:2+nMoments, iE3) = Phi_out(0:nMoments-1) with phase space and blocking
      
      Opacity_array(1, iE3) = E3
      Opacity_array(2, iE3) = wE3
      
      DO l = 0, nMoments-1
        ! Include phase space E3^2, conversion factor, and final state blocking
        Opacity_array(3+l, iE3) = Phi_out(l) * conv_fac * (1.0d0 - f_nu3)
      END DO
      
    END DO

    DEALLOCATE( xa_E3, wa_E3, xa_mu, wa_mu, Phi_out )

  END SUBROUTINE GeneralAngleAveragedOutgoingScattering_Moments


  ! ========================================================================
  ! Helper function to compute Legendre polynomials P_l(x)
  ! ========================================================================
  FUNCTION LegendrePolynomial(l, x) RESULT(Pl)
    INTEGER,  INTENT(IN) :: l
    REAL(DP), INTENT(IN) :: x
    REAL(DP) :: Pl
    
    REAL(DP) :: P0, P1, Pn
    INTEGER  :: n
    
    ! Recurrence relation: (n+1)P_{n+1} = (2n+1)xP_n - nP_{n-1}
    
    IF (l == 0) THEN
      Pl = 1.0d0
    ELSE IF (l == 1) THEN
      Pl = x
    ELSE
      P0 = 1.0d0
      P1 = x
      DO n = 1, l-1
        Pn = ( (2.0d0*n + 1.0d0)*x*P1 - n*P0 ) / (n + 1.0d0)
        P0 = P1
        P1 = Pn
      END DO
      Pl = Pn
    END IF
    
  END FUNCTION LegendrePolynomial

  !================================================================================
  !
  !  GeneralScatteringRout
  !
  !  Compute R_out(E1,E3,costheta) for the forward outscattering process
  !    nu_l2(E1) + l2(E2) -> nu_l4(E3) + l4(E4)
  !  using the analytic Fermi-integral form of Guo Eqs. 7-11.
  !
  !  R_out(E1,E3,mu) = lambda_1 * (1/16*pi*Delta^5) * (A1*I_2 + B1*I_1 + C1*I_0)
  !                                                                    (Guo Eq. 3)
  !
  !  I_s = integral_{E-}^{inf} dE2  E2^s  f2(E2)  [1-f4(E1+E2-E3)]
  !                                                                    (Guo Eq. 6)
  !
  !  Kinematic variables (Guo Eqs. 4-5):
  !    Delta = sqrt(E1^2 - 2*E1*E3*mu + E3^2)
  !    Q     = (m4^2 - m2^2) / 2
  !    k     = Q / [E1*E3*(1-mu)]
  !    E-    = (1/2)*{(E3-E1)*(1+k) + Delta*sqrt[(1+k)^2 + 2*m2^2/(E1*E3*(1-mu))]}
  !            max(E-, m2)
  !
  !  Coefficients A1, B1, C1: Guo Eqs. A6a-c (see ABC_Scattering_Coefficients).
  !  I_s computed analytically via Compute_Is_Integrals (Guo Eqs. 7-11).
  !
  !================================================================================
  SUBROUTINE GeneralScatteringRout( E1, E3, costheta, xT, xMu_l2, xMu_l4, m2, m4, lam1, lam2, lam3, Rout )

    REAL(DP), INTENT(IN)  :: E1, E3, costheta, xT, xMu_l2, xMu_l4, m2, m4, lam1, lam2, lam3
    REAL(DP), INTENT(OUT) :: Rout

    REAL(DP) :: mu              ! = costheta
    REAL(DP) :: Delta, Delta5
    REAL(DP) :: Q, k_val, E_minus
    REAL(DP) :: A1, B1, C1, A2, B2, C2, A3, B3, C3
    REAL(DP) :: I0_val, I1_val, I2_val
    REAL(DP) :: R1, R2, R3
    REAL(DP) :: disc

    Rout = 0.0d0
    mu  = costheta

    ! --- Q = (m4^2 - m2^2)/2 (Guo Eq. 5, generalized) ---
    Q = 0.5d0 * (m4**2 - m2**2)   ! Guo Eq. 5

    ! --- Delta (Guo Eq. 4) ---
    Delta = SQRT(MAX(E1**2 - 2.0d0*E1*E3*mu + E3**2, 0.0d0))   ! Guo Eq. 4
    IF (Delta < 1.0d-10) RETURN   ! pathological collinear case

    ! --- k (Guo Eq. 5): diverges at forward scattering mu=1 ---
    IF (1.0d0 - mu < 1.0d-10) RETURN   ! mu -> 1: E_minus -> inf, no phase space

    k_val = Q / (E1 * E3 * (1.0d0 - mu))   ! Guo Eq. 5

    ! --- E_minus (Guo Eq. 4): lower kinematic bound on l2 energy ---
    !   E- = (1/2)*{(E3-E1)*(1+k) + Delta*sqrt[(1+k)^2 + 2*m2^2/(E1*E3*(1-mu))]}
    disc = (1.0d0 + k_val)**2 + 2.0d0 * m2**2 / (E1 * E3 * (1.0d0 - mu)) ! Guo Eq. 4
    IF (disc < 0.0d0) RETURN
    E_minus = 0.5d0 * ( (E3 - E1)*(1.0d0 + k_val) + Delta*SQRT(disc) )   ! Guo Eq. 4
    E_minus = MAX(E_minus, m2)                                           ! Guo Eq. 4: floor at m2

    ! --- Coefficients A1, B1, C1 (Guo Eqs. A6a-c) ---
    CALL ABC_Scattering_Coefficients( E1, E3, mu, Q, m2, m4, A1, B1, C1, A2, B2, C2, A3, B3, C3 )

    ! --- I_s integrals (Guo Eqs. 7-11) ---
    CALL Compute_Is_Integrals( E1, E3, E_minus, xT, xMu_l2, xMu_l4, I0_val, I1_val, I2_val )
    ! This one does an explicit integration to check that our implementation is correct
    ! CALL Compute_Is_Integrals_explicitly( E1, E3, E_minus, xT, xMu_l2, xMu_l4, m4, I0_val, I1_val, I2_val )

    ! --- Assemble R_out (Guo Eq. 3 * lambda_1) ---
    Delta5 = Delta**5
    IF (Delta5 < 1.0d-50) RETURN

    R1 = 1.0d0 / (16.0d0 * pi * Delta5) &
        * (A1*I2_val + B1*I1_val + C1*I0_val)   ! Guo Eq. 3
    R2 = 1.0d0 / (16.0d0 * pi * Delta5) &
        * (A2*I2_val + B2*I1_val + C2*I0_val)   ! Guo Eq. 3
    R3 = 1.0d0 / (16.0d0 * pi * Delta5) &
        * (A3*I2_val + B3*I1_val + C3*I0_val)   ! Guo Eq. 3

    Rout = lam1*R1 + lam2*R2 + lam3*R3

    Rout = MAX(Rout, 0.0d0)

  END SUBROUTINE GeneralScatteringRout

  SUBROUTINE InverseMuonDecayRout( E1, E3, costheta, xT, xMu_l2, xMu_l4, m2, m4, lam1, lam2, lam3, Rout )

    REAL(DP), INTENT(IN)  :: E1, E3, costheta, xT, xMu_l2, xMu_l4, m2, m4, lam1, lam2, lam3
    REAL(DP), INTENT(OUT) :: Rout

    REAL(DP) :: mu              ! = costheta
    REAL(DP) :: DeltaA, Delta5
    REAL(DP) :: Q, k_val, k0, E_minus, E_plus
    REAL(DP) :: A1, B1, C1, A2, B2, C2, A3, B3, C3
    REAL(DP) :: I0_m, I1_m, I2_m, I0_p, I1_p, I2_p, I0_val, I1_val, I2_val
    REAL(DP) :: R1, R2, R3
    REAL(DP) :: disc

    Rout = 0.0d0
    mu  = costheta

    k0 = (m4 + m2) / (m4 - m2)

    ! --- Q = (m4^2 - m2^2)/2 (Guo Eq. 5, generalized) ---
    Q = 0.5d0 * (m4**2 - m2**2)   ! Guo Eq. 5

    ! --- Delta (Guo Eq. 17) ---
    DeltaA = SQRT(MAX(E1**2 + 2.0d0*E1*E3*mu + E3**2, 0.0d0))   ! Guo Eq. 4
    !IF (DeltaA < 1.0d-10) RETURN   ! pathological collinear case

    ! --- k (Guo Eq. 5): diverges at forward scattering mu=1 ---
    !IF (1.0d0 - mu < 1.0d-10) RETURN   ! mu -> 1: E_minus -> inf, no phase space

    k_val = Q / (E1 * E3 * (1.0d0 - mu))   ! Guo Eq. 5

    IF (k_val < k0) THEN
       Rout = 0.0d0   ! Theta function eq 16
       RETURN
    ENDIF

    ! --- E_minus (Guo Eq. 17): lower and upper kinematic bound on l2 energy ---
    disc = (1.0d0 - k_val)**2 - 2.0d0 * m2**2 / (E1 * E3 * (1.0d0 - mu)) ! Guo Eq. 4
    IF (disc < 0.0d0) RETURN
    E_minus = 0.5d0 * ( (E3 + E1)*(k_val - 1.0d0) - DeltaA*SQRT(disc) )   ! Guo Eq. 4
    E_plus  = 0.5d0 * ( (E3 + E1)*(k_val - 1.0d0) + DeltaA*SQRT(disc) )   ! Guo Eq. 4
    E_minus = MAX(E_minus, m2)
    E_plus  = MAX(E_plus, m2)
    IF (E_minus >= E_plus) THEN
      Rout = 0.0d0
      RETURN
    END IF

    ! --- Coefficients A1, B1, C1 (Guo Eqs. A7) --- 
    CALL ABC_Scattering_Coefficients( E1, -E3, mu, Q, m2, m4, A1, B1, C1, A2, B2, C2, A3, B3, C3 )

    ! --- I_s integrals (Guo Eqs. 18-19) ---
    CALL Compute_Is_Integrals( E1, -E3, E_minus, xT, xMu_l2, xMu_l4, I0_m, I1_m, I2_m )
    CALL Compute_Is_Integrals( E1, -E3, E_plus , xT, xMu_l2, xMu_l4, I0_p, I1_p, I2_p )
    ! Eq. 18
    I0_val = I0_m - I0_p
    I1_val = I1_m - I1_p
    I2_val = I2_m - I2_p

    ! --- Assemble R_out (Guo Eq. 3 * lambda_1) ---
    Delta5 = DeltaA**5
    IF (Delta5 < 1.0d-50) RETURN

    R1 = 1.0d0 / (16.0d0 * pi * Delta5) &
        * (A1*I2_val + B1*I1_val + C1*I0_val)   ! Guo Eq. 3
    R2 = 1.0d0 / (16.0d0 * pi * Delta5) &
        * (A2*I2_val + B2*I1_val + C2*I0_val)   ! Guo Eq. 3
    R3 = 1.0d0 / (16.0d0 * pi * Delta5) &
        * (A3*I2_val + B3*I1_val + C3*I0_val)   ! Guo Eq. 3

    Rout = lam1*R1 + lam2*R2 + lam3*R3
    Rout = MAX(Rout, 0.0d0)

  END SUBROUTINE InverseMuonDecayRout

  !================================================================================
  !
  !  Compute_Is_Integrals
  !
  !  Compute the I_s integrals (s=0,1,2) analytically using Guo Eqs. 7-11.
  !
  !  Starting from Guo Eq. 6:
  !    I_s = integral_{E-}^{inf} dE2  E2^s  f2(E2)  [1 - f4(E2 + delta)]
  !  where delta = E1 - E3 (energy transfer), and f4(E4) = FD(E4, xMu_l4, T)
  !  with E4 = E2 + delta.
  !
  !  Using the identity (Guo Eq. 7 preamble):
  !    f2(E2) [1 - f4(E2+delta)] = [f2(E2) - f4(E2+delta)] / (1 - exp(-xi))
  !  where xi = (delta + xMu_l2 - xMu_l4) / T = (E1-E3+xMu_l2-xMu_l4)/T,
  !  the integral splits into two incomplete Fermi-Dirac integrals:
  !
  !    I_s = T^(s+1) / (1 - exp(-xi))
  !          * [ F_s(eta2, y) - F_s(eta4', y) ]              (Guo Eqs. 7a-c)
  !
  !  where:
  !    y       = E_minus / T                    (lower limit in units of T)
  !    eta2    = xMu_l2 / T                        (l2 degeneracy)
  !    eta4'   = (xMu_l4 - delta) / T             (shifted l4 degeneracy)
  !    F_s(eta, y) = integral_y^inf x^s / (exp(x-eta)+1) dx    (Guo Eqs. 8-10)
  !
  !  On-shell guard: the l4 particle must satisfy E4 = E2 + delta >= m4.
  !  This is enforced by raising the lower limit:
  !    y_eff = max(y, (m4 - delta)/T)  when delta < m4
  !  so that E2 >= m4 - delta is guaranteed.
  !
  !  Degenerate case xi -> 0 (chemical equilibrium):
  !    lim_{xi->0} [F_s(eta2,y) - F_s(eta4',y)] / (1-exp(-xi))
  !    = F_{s-1}(eta2, y)   [via L'Hopital]
  !  handled when |xi| < xi_tol.
  !
  !  Inputs:
  !    E1, E3    -- neutrino energies [MeV]
  !    E_minus   -- lower kinematic bound on E2 [MeV] (from Guo Eq. 4)
  !    xT        -- temperature [MeV]
  !    xMu_l2      -- l2 chemical potential [MeV]
  !    xMu_l4      -- l4 chemical potential [MeV]
  !    m4        -- l4 mass [MeV] (for on-shell guard)
  !
  !  Outputs:
  !    I0, I1, I2 -- integrals [MeV^1, MeV^2, MeV^3]
  !
  !================================================================================
  !================================================================================
  ! SUBROUTINE Compute_Is_Integrals
  ! Evaluates I_s (s=0, 1, 2) defined in Guo 2020 Eq. 6:
  ! I_s = \int_{E_-}^{\infty} dE_2 E_2^s f_2(E_2) [1 - f_4(E_2 + \Delta)]
  !================================================================================
  SUBROUTINE Compute_Is_Integrals( E1, E3, E_minus, T, xMu_l2, xMu_l4, I0, I1, I2 )

    REAL(DP), INTENT(IN)  :: E1, E3, E_minus, T, xMu_l2, xMu_l4
    REAL(DP), INTENT(OUT) :: I0, I1, I2

    REAL(DP), PARAMETER :: xi_tol = 1.0d-8   ! Threshold for \xi -> 0 branch
    INTEGER , PARAMETER :: n_pts  = 32       ! GL quadrature points

    REAL(DP) :: xi, y, eta, eta_prime
    REAL(DP) :: F0_eta, F1_eta, F2_eta
    REAL(DP) :: F0_etap, F1_etap, F2_etap
    REAL(DP) :: fg, df0_dz

    I0 = 0.0d0; I1 = 0.0d0; I2 = 0.0d0

    ! --- Kinematic parameters ---
    ! Notice that eta_prime = eta - (E1 - E3 + xMu_l2 - xMu_l4) / T
    eta       = xMu_l2 / T
    eta_prime = (E3 - E1 + xMu_l4) / T
    xi        = eta_prime - eta

    ! --- Lower Integration Limit (y = E_- / T) ---
    y = E_minus / T

    SELECT CASE (CompleteFDIntegralsEvaluationMethod)
    CASE (1)
      IF (.NOT. GauLegInitialized) THEN
        CALL gauleg(-1.0d0, 1.0d0, xa_FI, wa_FI, nGL_FI)
        GauLegInitialized = .TRUE.
      ENDIF

      ! --- Compute raw complete Fermi-Dirac integrals F_s(z) = F_s(\eta - y) ---
      F0_eta = FI0   (    eta - y )
      F1_eta = Num_Fn( 1, eta - y, xa_FI, wa_FI, nGL_FI)
      F2_eta = Num_Fn( 2, eta - y, xa_FI, wa_FI, nGL_FI)

      F0_etap = FI0   (    eta_prime - y )
      F1_etap = Num_Fn( 1, eta_prime - y, xa_FI, wa_FI, nGL_FI)
      F2_etap = Num_Fn( 2, eta_prime - y, xa_FI, wa_FI, nGL_FI)

    CASE (2)

      CALL Init_Polylogarithms()
      ! --- Evaluate F_s(eta, y) and F_s(eta_prime, y) using polylog-based functions ---
      CALL Fermi_Dirac_Integrals( eta       - y, F0_eta , F1_eta , F2_eta )
      CALL Fermi_Dirac_Integrals( eta_prime - y, F0_etap, F1_etap, F2_etap )

    CASE (3)

      ! --- Evaluate F_s(eta, y) and F_s(eta_prime, y) using Fukushima (2015) routines ---
      F0_eta = F0z_Fukushima( eta - y )
      F1_eta = F1z_Fukushima( eta - y )
      F2_eta = F2z_Fukushima( eta - y )

      F0_etap = F0z_Fukushima( eta_prime - y )
      F1_etap = F1z_Fukushima( eta_prime - y )
      F2_etap = F2z_Fukushima( eta_prime - y )

    END SELECT

    IF (ABS(xi) < xi_tol) THEN
      ! WRITE(*,*) 'Warning: Compute_Is_Integrals: |xi| < xi_tol; using degenerate limit for I_s integrals.'
      ! =====================================================================
      ! SPECIAL CASE: \xi -> 0 (Guo 2020 Eqs. 11a, 11b, 11c)
      ! =====================================================================
      ! df0_dz = 1.0d0 / (EXP(MIN(y - eta, MaxExponent)) + 1.0d0)
      ! df0_dz = 1.0d0 / (EXP(y - eta) + 1.0d0)
      df0_dz = EXP(MIN(eta - y, 30.0d0)) / (1.0_dp + EXP(MIN(eta - y, 30.0d0)))

      ! Using F_n'(z) relations from the text (Eq 11):
      I0 = (T)    * df0_dz                                                  ! Eq. 11a
      I1 = (T**2) * ( F0_eta + y * df0_dz )                                 ! Eq. 11b
      I2 = (T**3) * ( 2.0d0 * F1_eta + 2.0d0 * y * F0_eta + y**2 * df0_dz ) ! Eq. 11c

    ELSE
      ! =====================================================================
      ! GENERAL CASE: Guo 2020 Eqs. 7a, 7b, and 7c
      ! =====================================================================

      fg = f_gamma( xi )
      ! Guo Eq. 7a
      I0 = T * fg * ( F0_etap - F0_eta )

      ! Guo Eq. 7b
      I1 = (T**2) * fg * ( (F1_etap - F1_eta) + &
            y * (F0_etap - F0_eta) )

      ! Guo Eq. 7c
      I2 = (T**3) * fg * ( (F2_etap - F2_eta) + &
          2.0d0 * y * (F1_etap - F1_eta) + (y**2) * (F0_etap - F0_eta) )

    END IF

  CONTAINS

  PURE FUNCTION f_gamma( z ) RESULT(res)
    REAL(DP), INTENT(IN) :: z
    REAL(DP) :: res

    IF (z > MaxExponent) THEN
      res = EXP(-z)                  ! Avoid overflow in exp(z)
    ELSE IF (z < -MaxExponent) THEN
      res = - 1.0d0                  ! Avoid overflow in exp(-z)
    ELSE
      res = 1.0d0 / (EXP(z) - 1.0d0) ! Exact definition
    END IF

  END FUNCTION f_gamma

  END SUBROUTINE Compute_Is_Integrals

  !================================================================================
  ! Helper Functions for Fermi-Dirac Integrals (Eq. 8, Eq. 9, Eq. 10)
  !================================================================================

  !--- F_0(\eta, y) = \ln(1 + \exp(\eta - y))  (Guo Eq. 8) ---
  PURE FUNCTION FI0( z ) RESULT(res)
    REAL(DP), INTENT(IN) :: z
    REAL(DP) :: res

    IF (z > MaxExponent) THEN
      res = z                          ! \ln(1 + e^arg) \approx arg
    ELSE IF (z < -MaxExponent) THEN
      res = EXP(z)                     ! \ln(1 + e^arg) \approx e^arg
    ELSE
      res = LOG(1.0d0 + EXP(z))        ! Exact Eq. 8
    END IF
  END FUNCTION FI0

  !--- Core Numerical Integrator for mapped FD integrand (Eq. 9 shifted by y) ---
  PURE FUNCTION Num_Fn( n, z, xa, wa, n_pts ) RESULT(res)
    INTEGER,  INTENT(IN) :: n 
    REAL(DP), INTENT(IN) :: z 
    REAL(DP), INTENT(IN) :: xa(:), wa(:) 
    INTEGER,  INTENT(IN) :: n_pts
    REAL(DP) :: res, x, w, arg, integrand
    INTEGER  :: i
    REAL(DP) :: xm1, xl1, xm2, xl2, t
    
    res = 0.0d0

    IF (z > 5.0d0) THEN
      ! Split GL quadrature domain at the Fermi surface (x=z) to resolve sharp drops
      xm1 = 0.5d0 * z
      xl1 = 0.5d0 * z
      xm2 = z + 15.0d0
      xl2 = 15.0d0

      DO i = 1, n_pts
        ! Domain 1: [0, z]
        t = xa(i)
        x = xm1 + xl1 * t
        w = wa(i) * xl1
        arg = x - z
        integrand = (x**n) / (EXP(arg) + 1.0d0)
        res = res + w * integrand

        ! Domain 2: [z, z+30]
        x = xm2 + xl2 * t
        w = wa(i) * xl2
        arg = x - z
        IF (arg > 500.0d0) THEN
          integrand = 0.0d0
        ELSE
          integrand = (x**n) / (EXP(arg) + 1.0d0)
        END IF
        res = res + w * integrand
      END DO

    ELSE
      ! Single domain mapping for non-degenerate regime
      xm1 = 0.5d0 * MAX(z + 100.0d0, 100.0d0)
      xl1 = 0.5d0 * MAX(z + 100.0d0, 100.0d0)
      
      DO i = 1, n_pts
        t = xa(i)
        x = xm1 + xl1 * t
        w = wa(i) * xl1
        arg = x - z
        IF (arg > MaxExponent) THEN
          integrand = 0.0d0
        ELSE
          integrand = (x**n) / (EXP(arg) + 1.0d0)
        END IF
        res = res + w * integrand
      END DO
    END IF
  END FUNCTION Num_Fn

  SUBROUTINE ABC_Scattering_Coefficients( E1, E3, mu, Q, m2, m4, A1, B1, C1, A2, B2, C2, A3, B3, C3 )

    REAL(DP), INTENT(IN)  :: E1, E3, mu, Q, m2, m4
    REAL(DP), INTENT(OUT) :: A1, B1, C1, A2, B2, C2, A3, B3, C3

    A1 = A1_f( E1, E3, mu )
    B1 = B1_f( E1, E3, mu, Q )
    C1 = C1_f( E1, E3, mu, Q, m2 )

    A2 = A1_f( -E3, -E1, mu )
    B2 = B1_f( -E3, -E1, mu, Q )
    C2 = C1_f( -E3, -E1, mu, Q, m2 )

    A3 = 0.0d0
    B3 = 0.0d0
    C3 = C3_f( E1, E3, mu, m2, m4 )

  END SUBROUTINE ABC_Scattering_Coefficients

  PURE ELEMENTAL FUNCTION A1_f( E1, E3, mu )
    REAL(DP), INTENT(IN) :: E1, E3, mu
    REAL(DP) :: A1_f

    ! --- A1 (Guo Eq. A6a) ---
    A1_f = E1*E3*(1.0d0 - mu)**2 * (E1**2 + E1*E3*(3.0d0+mu) + E3**2)   ! Guo Eq. A6a

  END FUNCTION A1_f

  PURE ELEMENTAL FUNCTION B1_f( E1, E3, mu, Q )
    REAL(DP), INTENT(IN) :: E1, E3, mu, Q
    REAL(DP) :: B1_f

    ! --- A1 (Guo Eq. A6a) ---
    B1_f = E1**2*E3*(1.0d0 - mu)**2 * ( 2.0d0*E1**2 + E1*E3*(3.0d0-mu) - E3**2*(1.0d0+3.0d0*mu) ) &
      + Q * (1.0d0 - mu) * ( E1**3 + E1**2*E3*(2.0d0+mu) - E1*E3**2*(2.0d0+mu) - E3**3 )   ! Guo Eq. A6b

  END FUNCTION B1_f

  PURE ELEMENTAL FUNCTION C1_f( E1, E3, mu, Q, m2 )
    REAL(DP), INTENT(IN) :: E1, E3, mu, Q, m2
    REAL(DP) :: C1_f

    ! Guo Eq. A6c
    C1_f = E1**3*E3*(1.0d0 - mu)**2                                                &
        * (E1**2 - 2.0d0*E1*E3*mu + E3**2 * (-0.5d0 + 1.5d0*mu**2) )               &
        + Q*E1*(1.0d0 - mu)                                                        &
        * ( E1**3 - E1**2*E3*mu + E1*E3**2*(-2.0d0 + mu**2) + E3**3*mu )           &
        + Q**2 * ( E1**2*mu - E1*E3*(1.5d0 + 0.5d0*mu**2) + E3**2*mu )             &
        + 0.5d0*E1*E3*(1.0d0 - mu**2) * (E1**2 - 2.0d0*E1*E3*mu + E3**2) * m2**2

  END FUNCTION C1_f

  PURE ELEMENTAL FUNCTION C3_f( E1, E3, mu, m2, m4 )
    REAL(DP), INTENT(IN) :: E1, E3, mu, m2, m4
    REAL(DP) :: C3_f

    ! Guo Eq. A6c
    ! Lohs has the extra m2*m4 Guo does not, typo in Guo!
    C3_f = (1.0d0 - mu) * (E1**2 - 2.0d0*E1*E3*mu + E3**2)**2 * m2*m4

  END FUNCTION C3_f

  !================================================================================
  !  gauleg: Gauss-Legendre quadrature nodes and weights (Numerical Recipes)
  !  Identical in structure to the version in wlSemiLeptonicOpacityModule2D.F90
  !================================================================================
  SUBROUTINE gauleg( x1, x2, x, w, n )

    INTEGER,  INTENT(IN)  :: n
    REAL(DP), INTENT(IN)  :: x1, x2
    REAL(DP), INTENT(OUT), DIMENSION(n) :: x, w

    INTEGER  :: i, j, m
    REAL(DP), PARAMETER :: eps = 3.0d-14
    REAL(DP) :: p1, p2, p3, pp, xl, xm, z, z1

    m  = (n+1)/2
    xm = 0.5d0*(x2 + x1)
    xl = 0.5d0*(x2 - x1)

    DO i = 1, m
      z = COS( pi * (i - 0.25d0) / (n + 0.5d0) )
      DO
        p1 = 1.0d0
        p2 = 0.0d0
        DO j = 1, n
          p3 = p2
          p2 = p1
          p1 = ((2.0d0*j - 1.0d0)*z*p2 - (j - 1.0d0)*p3) / j
        END DO
        pp = n * (z*p1 - p2) / (z*z - 1.0d0)
        z1 = z
        z  = z1 - p1/pp
        IF (ABS(z - z1) <= eps) EXIT
      END DO
      x(i)        = xm - xl*z
      x(n+1-i)    = xm + xl*z
      w(i)        = 2.0d0 * xl / ((1.0d0 - z*z) * pp*pp)
      w(n+1-i)    = w(i)
    END DO

  END SUBROUTINE gauleg

  !=======================================================================
  ! Subroutine: gaulag
  ! Purpose: Computes the abscissas (x) and weights (w) for generalized 
  !          Gauss-Laguerre quadrature: int_0^\infty x^alf e^-x f(x) dx
  !=======================================================================
  SUBROUTINE gaulag(x, w, n, alf)
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: n
    REAL(8), INTENT(IN)  :: alf
    REAL(8), INTENT(OUT) :: x(n), w(n)
    
    INTEGER :: i, its, j
    REAL(8) :: ai, p1, p2, p3, pp, z, z1
    REAL(8), PARAMETER :: EPS = 1.0d-14
    INTEGER, PARAMETER :: MAXIT = 50
    
    DO i = 1, n
      ! Initial guesses for the roots
      IF (i == 1) THEN
        z = (1.0d0 + alf) * (3.0d0 + 0.92d0 * alf) / (1.0d0 + 2.4d0 * n + 1.8d0 * alf)
      ELSE IF (i == 2) THEN
        z = z + (15.0d0 + 6.25d0 * alf) / (1.0d0 + 0.9d0 * alf + 2.5d0 * n)
      ELSE
        ai = i - 2
        z = z + ((1.0d0 + 2.55d0 * ai) / (1.9d0 * ai) + 1.26d0 * ai * alf / &
            (1.0d0 + 3.5d0 * ai)) * (z - x(i-2)) / (1.0d0 + 0.3d0 * alf)
      END IF
      
      ! Newton-Raphson loop to refine the root
      DO its = 1, MAXIT
        p1 = 1.0d0
        p2 = 0.0d0
        DO j = 1, n
          p3 = p2
          p2 = p1
          p1 = ((2.0d0 * j - 1.0d0 + alf - z) * p2 - (j - 1.0d0 + alf) * p3) / j
        END DO
        
        ! p1 is the Laguerre polynomial L_n^alf(z), pp is its derivative
        pp = (n * p1 - (n + alf) * p2) / z
        z1 = z
        z  = z1 - p1 / pp
        IF (ABS(z - z1) <= EPS) EXIT
      END DO
      
      IF (its > MAXIT) WRITE(*,*) ABS(z - z1)
      IF (its > MAXIT) PRINT *, 'WARNING: too many iterations in gaulag'
      
      ! Store root and compute weight
      x(i) = z
      w(i) = -EXP(gammln(alf + n) - gammln(REAL(n, 8))) / (pp * n * p2)
    END DO
  END SUBROUTINE gaulag

  !=======================================================================
  ! Function: gammln
  ! Purpose: Returns the value ln[Gamma(xx)] for xx > 0.
  !=======================================================================
  FUNCTION gammln(xx)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: xx
    REAL(8) :: gammln
    REAL(8) :: ser, tmp, x, y
    REAL(8), SAVE :: cof(6) = (/ 76.18009172947146d0, -86.50532032941677d0, &
                                24.01409824083091d0, -1.231739572450155d0, &
                                  0.1208650973866179d-2, -0.5395239384953d-5 /)
    INTEGER :: j
    
    x = xx
    y = x
    tmp = x + 5.5d0
    tmp = (x + 0.5d0) * LOG(tmp) - tmp
    ser = 1.000000000190015d0
    DO j = 1, 6
      y = y + 1.0d0
      ser = ser + cof(j) / y
    END DO
    gammln = tmp + LOG(2.5066282746310005d0 * ser / x)
  END FUNCTION gammln
  
  !================================================================================
  !  FD: Fermi-Dirac distribution with overflow guard
  !  Identical in structure to the version in wlSemiLeptonicOpacityModule2D.F90
  !================================================================================
  PURE ELEMENTAL FUNCTION FD( E, mu, T )
    REAL(DP), INTENT(IN) :: E, mu, T
    REAL(DP) :: FD

    FD = 1.0d0 / (EXP(MIN((E - mu)/T, MaxExponent)) + 1.0d0)

  END FUNCTION FD

  ! Notice that in Bollig each pieace in eq. 6.1 is 
  ! (Gw_MeV**2/2/pi**2)*alpha_1*R_1 = alpha_1*(2*pi/Delta**5) * (A1*I_2 + B1*I_1 + C1*I_0)
  ! where alpha_1 = (CV + CA)**2, alpha_2 = (CV - CA)**2, and alpha_3 = CA**2 - CV**2
  ! And in Guo each piece in eq. 2 is
  ! lambda_1 * (1/16/pi*Delta**5) * (A1*I_2 + B1*I_1 + C1*I_0)
  ! Hence the routine below
  PURE SUBROUTINE CVCA_to_lambda(CV, CA, lam1, lam2, lam3)
    REAL(DP), INTENT(IN)  :: CV, CA
    REAL(DP), INTENT(OUT) :: lam1, lam2, lam3

    lam1 = 16.0d0 * Gw_MeV**2 * (CV + CA)**2
    lam2 = 16.0d0 * Gw_MeV**2 * (CV - CA)**2
    lam3 = 16.0d0 * Gw_MeV**2 * (CA**2 - CV**2)

  END SUBROUTINE CVCA_to_lambda

  ! These are used to test that our implementation of the Is integrals is correct
  SUBROUTINE Integrate_Is_GL(E1,E3,T,mu2,mu4,a,b,n,I0,I1,I2)
    REAL(DP), INTENT(IN)  :: E1,E3,T,mu2,mu4,a,b
    INTEGER , INTENT(IN)  :: n
    REAL(DP), INTENT(OUT) :: I0,I1,I2

    REAL(DP), ALLOCATABLE :: x(:), w(:)
    REAL(DP) :: E2, wt, f2, f4, E4, kern
    INTEGER :: i

    I0 = 0.0d0; I1 = 0.0d0; I2 = 0.0d0
    IF (b <= a) RETURN

    ALLOCATE(x(n),w(n))
    CALL gauleg(a,b,x,w,n)

    DO i=1,n
      E2 = x(i)
      wt = w(i)

      f2 = FD(E2, mu2, T)

      E4 = E1 + E2 - E3
      IF (E4 <= 0.0d0) CYCLE
      f4 = FD(E4, mu4, T)

      kern = f2 * (1.0d0 - f4)

      I0 = I0 + wt * kern
      I1 = I1 + wt * (E2)    * kern
      I2 = I2 + wt * (E2*E2) * kern
    END DO

    DEALLOCATE(x,w)
  END SUBROUTINE Integrate_Is_GL

  SUBROUTINE Compute_Is_Integrals_explicitly( E1, E3, E_minus, T, xMu_l2, xMu_l4, m4, I0, I1, I2 )

    REAL(DP), INTENT(IN)  :: E1, E3, E_minus, T, xMu_l2, xMu_l4, m4
    REAL(DP), INTENT(OUT) :: I0, I1, I2

    INTEGER, PARAMETER :: n1 = 64, n2 = 64
    REAL(DP) :: delta, E2min, E2max, Ecut
    REAL(DP) :: I0a, I1a, I2a, I0b, I1b, I2b

    delta = E1 - E3

    ! Kinematic and on-shell lower bound:
    E2min = MAX(E_minus, m4 - delta)
    IF (E2min < 0.0d0) E2min = 0.0d0

    ! Practical upper bound; make generous for high-T / degenerate cases
    Ecut  = MAX( 80.0d0*T, 200.0d0 )
    E2max = E2min + Ecut + MAX(xMu_l2,0.0d0) + MAX(xMu_l4-delta,0.0d0)

    ! Split integration to better resolve near-threshold + tail
    CALL Integrate_Is_GL(E1,E3,T,xMu_l2,xMu_l4,E2min, 0.5d0*(E2min+E2max), n1, I0a,I1a,I2a)
    CALL Integrate_Is_GL(E1,E3,T,xMu_l2,xMu_l4,0.5d0*(E2min+E2max), E2max, n2, I0b,I1b,I2b)

    I0 = I0a + I0b
    I1 = I1a + I1b
    I2 = I2a + I2b

  END SUBROUTINE Compute_Is_Integrals_explicitly

  ! Deal with collinear case
  SUBROUTINE CalculateCollinearH0(E1, E2, H0I, H0II)

    REAL(DP), INTENT(IN)  :: E1, E2
    REAL(DP), INTENT(OUT) :: H0I, H0II

    ! Notice that E3 == E1 in the collinear case. This is Eq. 44 and 45 from 
    ! Mezzacappa & Bruenn (1993) where E1 = E0 and E3 = E0' and E2 = Ee
    IF (E1 >= E2) THEN
      H0I  = 4.0d0/15.0d0*E2**5 + 4.0d0/3.0d0*E2**4*E1 + 8.0d0/3.0d0*E2**3*E1**2
      H0II = 4.0d0/15.0d0*E2**5 - 4.0d0/3.0d0*E2**4*E1 + 8.0d0/3.0d0*E2**3*E1**2
    ELSE
      H0I  = 4.0d0/15.0d0*E1**5 + 4.0d0/3.0d0*E1**4*E2 + 8.0d0/3.0d0*E1**3*E2**2
      H0II = 4.0d0/15.0d0*E1**5 - 4.0d0/3.0d0*E1**4*E2 + 8.0d0/3.0d0*E1**3*E2**2
    END IF

  END SUBROUTINE CalculateCollinearH0

END MODULE wlGeneralLeptonScatteringModule
