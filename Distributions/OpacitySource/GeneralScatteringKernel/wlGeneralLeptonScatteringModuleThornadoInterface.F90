MODULE wlGeneralLeptonScatteringModuleThornadoInterface

  USE wlKindModule,         ONLY: dp
  USE wlEosConstantsModule, ONLY: &
    pi, Gw_MeV, hbarc, sin2W, me, mmu, cvel
  USE wlGeneralLeptonScatteringModule, ONLY: &
    SelectLambdaFromProcessIndex, &
    DetBal_MaxExponent,           &
    Compute_Is_Integrals, FD,     &
    Compute_Is_Integrals_FD,      &
    f_gamma, xi_tol,              &
    A1_f, B1_f, C1_f, C3_f,       &
    gauleg, LegendrePolynomial,   &
    gaulag, CalculateCollinearH0

  IMPLICIT NONE
  PRIVATE

  INTEGER , PARAMETER, PUBLIC :: iProcessMin_Default = 1
  INTEGER , PARAMETER, PUBLIC :: iProcessMax_Default = 36
  INTEGER , PARAMETER, PUBLIC :: nProcessTotal       = iProcessMax_Default

  ! ==========================================================================
  ! Process selection.
  !
  ! The selection is an explicit ORDERED LIST, not a range: any subset of the
  ! 36 processes can be requested, in any combination. ProcessList defines the
  ! row order of Phout / Rout, so the caller controls the output layout.
  !
  ! iProcessMin / iProcessMax are kept only as the min/max of the list, for
  ! diagnostics -- no logic may assume the selection is contiguous.
  ! ==========================================================================
  INTEGER                     :: nProcessSel
  INTEGER , ALLOCATABLE       :: ProcessList(:)             ! row -> absolute process
  INTEGER                     :: iRowOfProcess(nProcessTotal) ! absolute -> row, 0 = unselected
  INTEGER                     :: iProcessMin
  INTEGER                     :: iProcessMax

  ! --- Named process groups, for callers that want to combine whole blocks ---
  INTEGER , PARAMETER, PUBLIC :: iGroup_NES    = 1   ! 1..12  NC on e-, e+
  INTEGER , PARAMETER, PUBLIC :: iGroup_NMS    = 2   ! 13..24 NC on mu-, mu+
  INTEGER , PARAMETER, PUBLIC :: iGroup_LFE    = 3   ! 25..28 lepton flavour exchange
  INTEGER , PARAMETER, PUBLIC :: iGroup_LFC    = 4   ! 29..32 lepton flavour conversion
  INTEGER , PARAMETER, PUBLIC :: iGroup_IMD    = 5   ! 33..36 inverse muon (33,34) and anti-muon (35,36) decay
  INTEGER , PARAMETER, PUBLIC :: nGroupTotal   = 5

  INTEGER , DIMENSION(nGroupTotal), PARAMETER :: GroupFirst = (/  1, 13, 25, 29, 33 /)
  INTEGER , DIMENSION(nGroupTotal), PARAMETER :: GroupLast  = (/ 12, 24, 28, 32, 36 /)
  REAL(DP), PARAMETER         :: ElectronMass        = me
  REAL(DP), PARAMETER         :: MuonMass            = mmu
  REAL(DP), PARAMETER         :: E2_max_Default      = 100.0d0
  INTEGER , PARAMETER         :: nE2_Default         = 32
  INTEGER                     :: nE2_Lepton

  ! --- Energy grids the coefficient tables were built on. Stored so the
  ! --- batched entry points need not be handed (iE1,E1) pairs that could
  ! --- disagree with the tables.
  INTEGER                  :: nE1_Grid, nE3_Grid
  REAL(DP), ALLOCATABLE    :: E1_Grid(:), E3_Grid(:)

  ! --- Maximum number of Legendre moments. Fixes the extent of every
  ! --- thread-private accumulator, so it must stay a PARAMETER.
  INTEGER , PARAMETER      :: nL_Max = 2

  ! ==========================================================================
  ! Angular quadrature and coefficient tables, PER MASS PAIR.
  !
  ! The cos(theta) integral is NOT done with one Gauss-Legendre rule on
  ! [-1,1]. At low T / high E the 2->2 integrand is a narrow forward peak at
  ! 1-mu ~ 1e-4..1e-3 that a 12-24 point rule on [-1,1] never samples, and
  ! IMD has a hard kinematic edge (k = k0, where disc = 0) inside [-1,1].
  ! Both are T-INDEPENDENT, so the rule is built at Init, separately for every
  ! (iE1, iE3, mass pair):
  !
  !   x = 1 - mu on [xLo, xHi]
  !     xHi = 2                                      (2->2)
  !         = min(2, Q / (k0 E1 E3))                 (IMD: threshold edge)
  !     xLo = max(xFloor, x where E_minus(x) = EmCap)
  !           (E_minus grows like 1/x or 1/sqrt(x) as x -> 0 off the
  !            diagonal; above EmCap the I_s integrals are ~exp(-(EmCap-mu)/T))
  !   Gauss-Legendre in s(x) = alpha*ln(x) + beta*x on [s(xLo), s(xHi)]:
  !     Q = 0 scattering (NES, NMS) : pure log      (alpha=1, beta=0)
  !     Q /= 0 scattering (LFE, LFC): linear in mu  (alpha=0, beta=1)
  !     IMD                         : log-linear    (alpha=1, beta=1/XcIMD)
  !   (choices from wlTestAngularQuadrature-style comparisons against a
  !    converged log-panel reference).
  !
  ! Because the nodes depend on the mass pair, A, B, C and the weights are
  ! stored per mass pair. Per case this streams the same 7 coefficients as the
  ! old shared (A_Scat, B_Scat(jFl2,jFl4), C_Scat(jFl2,jFl4)) layout, plus the
  ! weight, so the bandwidth per (iE1,iE3,iTh,case) is unchanged up to that.
  ! Slots: A(1:2), B(1:2) and C(1:2) multiply lam1 / lam2 (R1 / R2); C(3)
  ! multiplies lam3 (R3) and is zero for IMD. 1/(16 pi Delta^5) is built in.
  ! ==========================================================================
  REAL(DP), ALLOCATABLE    :: A_MP  (:, :, :, :, :)   ! (nE1,nE3,nTheta,2,nMassPair)
  REAL(DP), ALLOCATABLE    :: B_MP  (:, :, :, :, :)   ! (nE1,nE3,nTheta,2,nMassPair)
  REAL(DP), ALLOCATABLE    :: C_MP  (:, :, :, :, :)   ! (nE1,nE3,nTheta,3,nMassPair)
  ! P_{l}(mu_i) * w_i (Jacobian of the x -> s map included)
  REAL(DP), ALLOCATABLE    :: PlW_MP(:, :, :, :, :)   ! (nE1,nE3,nTheta,nL_Max,nMassPair)

  REAL(DP), PARAMETER      :: XFloor_Scat = 1.0d-10   ! smallest 1-mu, 2->2 (diagonal E1=E3)
  REAL(DP), PARAMETER      :: XFloor_IMD  = 1.0d-7    ! smallest 1-mu, IMD (x*R ~ x there)
  REAL(DP), PARAMETER      :: EmCap       = 3.0d3     ! [MeV] lepton-2 energy beyond which I_s ~ 0
  REAL(DP), PARAMETER      :: XcIMD       = 0.1d0     ! log -> linear crossover of the IMD map

  ! ==========================================================================
  ! Tabulated kinematics. E_minus, E_plus and the phase-space masks depend
  ! only on (E1,E3,costh) and the case MASS PAIR -- NOT on T, Mu_e or Mu_mu --
  ! so they are invariant for the whole run and are built once at Init. This
  ! is what keeps every SQRT, divide and data-dependent branch out of the
  ! kernel.
  !
  ! Indexed by deduplicated mass pair, not by case: the 9 distinct cases share
  ! at most 5 distinct (m2,m4,isIMD) combinations, and an electron-only run
  ! collapses to a single one. Keeps the added traffic to ~7 doubles per
  ! triple instead of 27.
  ! ==========================================================================
  INTEGER , PARAMETER      :: nMassPairMax = 5
  INTEGER                  :: nMassPair
  REAL(DP), ALLOCATABLE    :: Eminus_MP(:, :, :, :)   ! (nE1,nE3,nTheta,nMassPair), at that pair's nodes
  REAL(DP), ALLOCATABLE    :: Eplus_MP (:, :, :, :)   ! ditto; IMD slot only
  LOGICAL , ALLOCATABLE    :: Active_MP(:, :, :, :)   ! .FALSE. = no phase space

  ! --- Per-case constants, all resolved at Init ---
  INTEGER , ALLOCATABLE    :: iMassPair_Case(:)  ! case -> mass-pair slot
  LOGICAL , ALLOCATABLE    :: IsIMD_Case(:)      ! inverse (anti)muon decay case
  LOGICAL , ALLOCATABLE    :: TwoSided_Case(:)   ! IMD: I_s(E-) - I_s(E+)
  REAL(DP), ALLOCATABLE    :: E3Sign_Case(:)     ! +1 scattering, -1 IMD (enters eta')
  INTEGER , ALLOCATABLE    :: iIdx2(:), iIdx4(:) ! flattened iABC_m2/iABC_m4
  INTEGER , ALLOCATABLE    :: iChemSpecies2(:)   ! 1 = electron, 2 = muon
  INTEGER , ALLOCATABLE    :: iChemSpecies4(:)
  REAL(DP), ALLOCATABLE    :: ChemSign2(:)       ! +1 for l-, -1 for l+
  REAL(DP), ALLOCATABLE    :: ChemSign4(:)

  ! IMD threshold (Guo Eq. 16): kernel vanishes for k < k0
  REAL(DP), PARAMETER      :: k0_IMD = (MuonMass + ElectronMass) / (MuonMass - ElectronMass)

  INTEGER                  :: nThetaScattering   ! nodes per (iE1,iE3,mass pair)
  REAL(DP), ALLOCATABLE    :: xa_GLeg_ref(:), wa_GLeg_ref(:)
  REAL(DP), ALLOCATABLE    :: xa_GLag(:), wa_GLag(:)

  INTEGER , PARAMETER      :: iABC_el = 1
  INTEGER , PARAMETER      :: iABC_mu = 2

  INTEGER , PARAMETER      :: nDistinctCasesMax = 10
  INTEGER                  :: nDistinctCases
  REAL(DP), DIMENSION(nDistinctCasesMax,2) :: MassPairDistinct
  INTEGER , PARAMETER      :: i_em_em = 1 ! particle 2 is e-  and particle 4 is e-
  INTEGER , PARAMETER      :: i_ep_ep = 2 ! particle 2 is e+  and particle 4 is e+
  INTEGER , PARAMETER      :: i_mm_mm = 3 ! particle 2 is mu- and particle 4 is mu-
  INTEGER , PARAMETER      :: i_mp_mp = 4 ! particle 2 is mu+ and particle 4 is mu+
  INTEGER , PARAMETER      :: i_em_mm = 5 ! particle 2 is e-  and particle 4 is mu-
  INTEGER , PARAMETER      :: i_mm_em = 6 ! particle 2 is mu- and particle 4 is e-
  INTEGER , PARAMETER      :: i_ep_mp = 7 ! particle 2 is e+  and particle 4 is mu+
  INTEGER , PARAMETER      :: i_mp_ep = 8 ! particle 2 is mu+ and particle 4 is e+
  INTEGER , PARAMETER      :: i_IMD   = 9 ! IMD: particle 2 is e-, particle 4 is mu- (applies to 33 and 34)
  INTEGER , PARAMETER      :: i_IMDbar = 10 ! inverse anti-muon decay: particle 2 is e+, particle 4 is mu+ (35 and 36)

  INTEGER , DIMENSION(8), PARAMETER :: &
    iABC_m2 = (/ iABC_el, iABC_el, iABC_mu, iABC_mu, &
                 iABC_el, iABC_mu, iABC_el, iABC_mu /)
  INTEGER , DIMENSION(8), PARAMETER :: &
    iABC_m4 = (/ iABC_el, iABC_el, iABC_mu, iABC_mu, &
                 iABC_mu, iABC_el, iABC_mu, iABC_el /)

  INTEGER                             :: iDistinctMap(iProcessMax_Default)
  ! Inverse of ChosenToTrueCaseMap: iProcess -> COMPACT case index (1..nDistinctCases).
  ! Lets the hot loops carry one small dense index instead of the sparse 1..9 one.
  INTEGER                             :: iProcessToCase(iProcessMax_Default)
  ! Row -> compact case, i.e. iProcessToCase(ProcessList(row)) collapsed into one
  ! gather. Keeps the hot loops (and the device region) from chasing two maps.
  INTEGER , ALLOCATABLE, DIMENSION(:) :: iRowToCase(:)
  INTEGER , ALLOCATABLE, DIMENSION(:) :: ChosenToTrueCaseMap(:)
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: lam1(:), lam2(:), lam3(:)

  ! --- Device declarations. CalculateAllPhout is a ROUTINE SEQ / DECLARE TARGET
  ! --- routine, and such a routine may only touch module data that carries its
  ! --- own device declaration -- a PRESENT clause at the caller does not cover
  ! --- it. The scalars are listed too: an undeclared module scalar reads as
  ! --- garbage/0 on device.
#if   defined(WEAKLIB_OMP_OL)
  !$OMP DECLARE TARGET( A_MP, B_MP, C_MP, PlW_MP,                            &
  !$OMP                 Eminus_MP, Eplus_MP, Active_MP,                       &
  !$OMP                 iMassPair_Case, IsIMD_Case, TwoSided_Case,            &
  !$OMP                 E3Sign_Case, iIdx2, iIdx4,                            &
  !$OMP                 lam1, lam2, lam3, iRowToCase,                         &
  !$OMP                 iChemSpecies2, iChemSpecies4,                         &
  !$OMP                 ChemSign2, ChemSign4, E1_Grid, E3_Grid,               &
  !$OMP                 nThetaScattering, nDistinctCases, nProcessSel )
#elif defined(WEAKLIB_OACC)
  !$ACC DECLARE CREATE( A_MP, B_MP, C_MP, PlW_MP,                            &
  !$ACC                 Eminus_MP, Eplus_MP, Active_MP,                       &
  !$ACC                 iMassPair_Case, IsIMD_Case, TwoSided_Case,            &
  !$ACC                 E3Sign_Case, iIdx2, iIdx4,                            &
  !$ACC                 lam1, lam2, lam3, iRowToCase,                         &
  !$ACC                 iChemSpecies2, iChemSpecies4,                         &
  !$ACC                 ChemSign2, ChemSign4, E1_Grid, E3_Grid,               &
  !$ACC                 nThetaScattering, nDistinctCases, nProcessSel )
#endif

  REAL(DP), PARAMETER :: conv_fac = 1.0d0 / ( (2.0d0 * pi)**3 * hbarc )
  ! NOTE: xi_tol is imported from wlGeneralLeptonScatteringModule, NOT redeclared
  ! here. The scalar path reaches it through Compute_Is_Integrals while the
  ! batched kernel evaluates the xi -> 0 test itself, so two copies could drift
  ! apart and make the two paths disagree on which branch to take.

  PUBLIC :: InitGeneralScatteringKernels
  PUBLIC :: FinalizeGeneralScatteringKernels
  PUBLIC :: CalculateAllRout
  PUBLIC :: CalculateAllPhout
  PUBLIC :: CalculateAllPhout_3D
  PUBLIC :: DetailedBalanceFactor
  PUBLIC :: CalculateAllRoutIntegrated
  PUBLIC :: ProcessRowIndex
  PUBLIC :: ProcessesInGroup
  PUBLIC :: SelectedProcessCount
  PUBLIC :: SelectedProcessList

CONTAINS

  !================================================================================
  !
  !  ProcessRowIndex
  !
  !  Row index into the Phout / Rout arrays for an ABSOLUTE process index -- the
  !  kind ProcessIndexFromReactionString returns.
  !
  !  Use this at call sites instead of any hand-rolled index arithmetic. The
  !  selection is an arbitrary subset, so "iProcess - iProcessMin + 1" is simply
  !  wrong unless the selection happens to be the contiguous block starting at
  !  iProcessMin; this is a direct lookup, and it fails loudly rather than
  !  silently reading the wrong row when the process was never selected.
  !
  !================================================================================
  FUNCTION ProcessRowIndex( iProcess ) RESULT( idx )

    INTEGER, INTENT(IN) :: iProcess
    INTEGER :: idx

    IF ( iProcess < 1 .OR. iProcess > nProcessTotal ) THEN
      WRITE(*,*) 'ProcessRowIndex: process index out of range: ', iProcess
      STOP
    ENDIF

    idx = iRowOfProcess(iProcess)

    IF ( idx == 0 ) THEN
      WRITE(*,*) 'ProcessRowIndex: process ', iProcess, ' was not selected at Init.'
      WRITE(*,*) '  selected processes: ', ProcessList(1:nProcessSel)
      STOP
    ENDIF

  END FUNCTION ProcessRowIndex

  !--- Number of rows in Phout / Rout for the current selection. --------------
  PURE FUNCTION SelectedProcessCount() RESULT( n )
    INTEGER :: n
    n = nProcessSel
  END FUNCTION SelectedProcessCount

  !--- The selected absolute process indices, in row order. ------------------
  PURE FUNCTION SelectedProcessList() RESULT( list )
    INTEGER :: list(nProcessSel)
    list = ProcessList(1:nProcessSel)
  END FUNCTION SelectedProcessList

  !================================================================================
  !  ProcessesInGroup: the absolute process indices belonging to a named group.
  !  Handy for assembling an iProcessList_Option by hand.
  !================================================================================
  SUBROUTINE ProcessesInGroup( iGroup, list, n )

    INTEGER, INTENT(IN)  :: iGroup
    INTEGER, INTENT(OUT) :: list(:)
    INTEGER, INTENT(OUT) :: n

    INTEGER :: i

    IF ( iGroup < 1 .OR. iGroup > nGroupTotal ) THEN
      WRITE(*,*) 'ProcessesInGroup: unknown group ', iGroup
      STOP
    ENDIF

    n = GroupLast(iGroup) - GroupFirst(iGroup) + 1

    IF ( SIZE(list) < n ) THEN
      WRITE(*,*) 'ProcessesInGroup: list too small, need ', n, ' got ', SIZE(list)
      STOP
    ENDIF

    DO i = 1, n
      list(i) = GroupFirst(iGroup) + i - 1
    ENDDO

  END SUBROUTINE ProcessesInGroup

  SUBROUTINE InitGeneralScatteringKernels(E1_array, E3_array, nE1, nE3, nTheta, &
                                          iProcessMin_Option, iProcessMax_Option, &
                                          E2_max_Option, nE2_Option, &
                                          iProcessList_Option, iGroupList_Option)

    REAL(DP), INTENT(IN)           :: E1_array(nE1), E3_array(nE3)
    INTEGER , INTENT(IN)           :: nE1, nE3, nTheta
    ! --- Process selection. Give AT MOST ONE of:
    !       iProcessList_Option : explicit absolute indices, in the row order wanted
    !       iGroupList_Option   : named groups (iGroup_NES, iGroup_LFE, ...)
    !       iProcessMin/Max     : a contiguous block (kept for backward compatibility)
    !     Nothing given -> the full default range.
    INTEGER , INTENT(IN), OPTIONAL :: iProcessMin_Option, iProcessMax_Option
    INTEGER , INTENT(IN), OPTIONAL :: iProcessList_Option(:)
    INTEGER , INTENT(IN), OPTIONAL :: iGroupList_Option(:)
    REAL(DP), INTENT(IN), OPTIONAL :: E2_max_Option
    INTEGER , INTENT(IN), OPTIONAL :: nE2_Option

    REAL(DP) :: E1, E3, costh, Delta, DeltaA, E2_max
    REAL(DP) :: Q, m2, m4, k_val, disc, Em, Ep, OneMinusMu, EEOneMinusMu
    REAL(DP) :: xLo, xHi, xFloor, x, wx, sAlpha, sBeta, Den
    REAL(DP), ALLOCATABLE :: sNode(:), sWght(:)
    INTEGER  :: iE1, iE3, iTh, iProcess, iCase, iTrueCase, iL
    LOGICAL  :: IsActive


    INTEGER, DIMENSION(nDistinctCasesMax) :: unique_cases
    INTEGER :: idx, current_val
    INTEGER :: iRow, iGroup, nInGroup, nForms
    INTEGER, DIMENSION(nProcessTotal) :: tmpList, groupBuf
    REAL(DP), DIMENSION(nMassPairMax) :: MPm2, MPm4
    LOGICAL , DIMENSION(nMassPairMax) :: MPisIMD

    ! A second Init without Finalize would re-ALLOCATE allocated arrays
    IF ( ALLOCATED(E1_Grid) ) CALL FinalizeGeneralScatteringKernels()

    nThetaScattering = nTheta
    nE1_Grid         = nE1
    nE3_Grid         = nE3

    ALLOCATE( E1_Grid(nE1), E3_Grid(nE3) )
    E1_Grid(:) = E1_array(:)
    E3_Grid(:) = E3_array(:)

    IF ( PRESENT(nE2_Option) ) THEN
      nE2_Lepton = nE2_Option
    ELSE
      nE2_Lepton = nE2_Default
    ENDIF

    IF ( PRESENT(E2_max_Option) ) THEN
      E2_max = E2_max_Option
    ELSE
      E2_max = E2_max_Default
    ENDIF

    ! =======================================================================
    ! Build the process selection. Exactly one input form may be used.
    ! =======================================================================
    nForms = 0
    IF ( PRESENT(iProcessList_Option) ) nForms = nForms + 1
    IF ( PRESENT(iGroupList_Option)   ) nForms = nForms + 1
    IF ( PRESENT(iProcessMin_Option) .OR. PRESENT(iProcessMax_Option) ) &
      nForms = nForms + 1

    IF ( nForms > 1 ) THEN
      WRITE(*,*) 'InitGeneralScatteringKernels: give at most ONE of ', &
                 'iProcessList_Option, iGroupList_Option, iProcessMin/Max_Option.'
      STOP
    ENDIF

    IF ( PRESENT(iProcessList_Option) ) THEN

      nProcessSel = SIZE(iProcessList_Option)
      IF ( nProcessSel < 1 ) THEN
        WRITE(*,*) 'InitGeneralScatteringKernels: empty iProcessList_Option.'
        STOP
      ENDIF
      IF ( nProcessSel > nProcessTotal ) THEN
        WRITE(*,*) 'InitGeneralScatteringKernels: iProcessList_Option has more than ', &
                   nProcessTotal, ' entries.'
        STOP
      ENDIF
      tmpList(1:nProcessSel) = iProcessList_Option(:)

    ELSE IF ( PRESENT(iGroupList_Option) ) THEN

      nProcessSel = 0
      DO idx = 1, SIZE(iGroupList_Option)
        iGroup = iGroupList_Option(idx)
        CALL ProcessesInGroup( iGroup, groupBuf, nInGroup )
        IF ( nProcessSel + nInGroup > nProcessTotal ) THEN
          WRITE(*,*) 'InitGeneralScatteringKernels: too many (or repeated) groups in ', &
                     'iGroupList_Option.'
          STOP
        ENDIF
        DO iRow = 1, nInGroup
          nProcessSel = nProcessSel + 1
          tmpList(nProcessSel) = groupBuf(iRow)
        ENDDO
      ENDDO
      IF ( nProcessSel < 1 ) THEN
        WRITE(*,*) 'InitGeneralScatteringKernels: empty iGroupList_Option.'
        STOP
      ENDIF

    ELSE

      ! Backward-compatible contiguous block.
      IF ( PRESENT(iProcessMin_Option) ) THEN
        iProcessMin = iProcessMin_Option
      ELSE
        iProcessMin = iProcessMin_Default
      ENDIF
      IF ( PRESENT(iProcessMax_Option) ) THEN
        iProcessMax = iProcessMax_Option
      ELSE
        iProcessMax = iProcessMax_Default
      ENDIF
      IF ( iProcessMax < iProcessMin ) THEN
        WRITE(*,*) 'InitGeneralScatteringKernels: iProcessMax < iProcessMin.'
        STOP
      ENDIF
      nProcessSel = iProcessMax - iProcessMin + 1
      DO idx = 1, nProcessSel
        tmpList(idx) = iProcessMin + idx - 1
      ENDDO

    ENDIF

    ! --- validate, then commit the selection and its inverse map ---
    iRowOfProcess(:) = 0
    DO idx = 1, nProcessSel
      iProcess = tmpList(idx)
      IF ( iProcess < 1 .OR. iProcess > nProcessTotal ) THEN
        WRITE(*,*) 'InitGeneralScatteringKernels: process out of range: ', iProcess
        STOP
      ENDIF
      IF ( iRowOfProcess(iProcess) /= 0 ) THEN
        WRITE(*,*) 'InitGeneralScatteringKernels: process ', iProcess, &
                   ' requested twice (rows ', iRowOfProcess(iProcess), ' and ', idx, ')'
        STOP
      ENDIF
      iRowOfProcess(iProcess) = idx
    ENDDO

    ALLOCATE( ProcessList(nProcessSel) )
    ProcessList(:) = tmpList(1:nProcessSel)

    ! Diagnostics only -- nothing below may assume contiguity.
    iProcessMin = MINVAL(ProcessList)
    iProcessMax = MAXVAL(ProcessList)

    ALLOCATE( lam1(nProcessSel) )
    ALLOCATE( lam2(nProcessSel) )
    ALLOCATE( lam3(nProcessSel) )

    ! --- Mass pairs per distinct case (m2 = particle 2, m4 = particle 4) ---
    MassPairDistinct(i_em_em,1) = ElectronMass
    MassPairDistinct(i_em_em,2) = ElectronMass

    MassPairDistinct(i_ep_ep,1) = ElectronMass
    MassPairDistinct(i_ep_ep,2) = ElectronMass

    MassPairDistinct(i_mm_mm,1) = MuonMass
    MassPairDistinct(i_mm_mm,2) = MuonMass

    MassPairDistinct(i_mp_mp,1) = MuonMass
    MassPairDistinct(i_mp_mp,2) = MuonMass

    MassPairDistinct(i_em_mm,1) = ElectronMass
    MassPairDistinct(i_em_mm,2) = MuonMass

    MassPairDistinct(i_mm_em,1) = MuonMass
    MassPairDistinct(i_mm_em,2) = ElectronMass

    MassPairDistinct(i_ep_mp,1) = ElectronMass
    MassPairDistinct(i_ep_mp,2) = MuonMass

    MassPairDistinct(i_mp_ep,1) = MuonMass
    MassPairDistinct(i_mp_ep,2) = ElectronMass

    MassPairDistinct(i_IMD,1) = ElectronMass
    MassPairDistinct(i_IMD,2) = MuonMass

    MassPairDistinct(i_IMDbar,1) = ElectronMass
    MassPairDistinct(i_IMDbar,2) = MuonMass

    ! =======================================================================
    ! Resolve the distinct-case list FIRST. The per-case coefficient and
    ! kinematic tables below are dimensioned by nDistinctCases, so this block
    ! has to precede them (it used to sit at the end of this routine).
    ! =======================================================================
    nDistinctCases  = 0
    unique_cases(:) = -1
    DO idx = 1, nProcessSel
       iProcess = ProcessList(idx)

       CALL SelectLambdaFromProcessIndex( iProcess, &
                  lam1(idx), &
                  lam2(idx), &
                  lam3(idx))

       ! Notice that iDistinctMap is always called with the "True" iProcess
       CALL SelectiDistinctFromProcessIndex( iProcess, &
                  iDistinctMap(iProcess))

       current_val = iDistinctMap(iProcess)

       ! 2. Check if current_val is NOT present in the populated portion of unique_cases
       ! This basically looks at whether we have encountered current_val already, and if
       ! we have not, then increases nDistinctCases by 1.
       IF (.NOT. ANY(unique_cases(1:nDistinctCases) == current_val)) THEN
          nDistinctCases = nDistinctCases + 1
          unique_cases(nDistinctCases) = current_val
       END IF
    ENDDO

    ALLOCATE( ChosenToTrueCaseMap(nDistinctCases) )
    ! Needed because the selected processes (a contiguous block, an explicit
    ! list or a set of groups) may use any subset of the distinct cases.
    DO iCase=1, nDistinctCases
      ChosenToTrueCaseMap(iCase) = unique_cases(iCase)
    ENDDO

    ! Invert it: every selected process gets its compact case index.
    iProcessToCase(:) = -1
    DO idx = 1, nProcessSel
      iProcess = ProcessList(idx)
      DO iCase = 1, nDistinctCases
        IF ( ChosenToTrueCaseMap(iCase) == iDistinctMap(iProcess) ) THEN
          iProcessToCase(iProcess) = iCase
          EXIT
        ENDIF
      ENDDO
      IF ( iProcessToCase(iProcess) < 0 ) THEN
        WRITE(*,*) 'InitGeneralScatteringKernels: no compact case for iProcess =', iProcess
        STOP
      ENDIF
    ENDDO

    ! Collapse the two maps into a single row -> case gather.
    ALLOCATE( iRowToCase(nProcessSel) )
    DO idx = 1, nProcessSel
      iRowToCase(idx) = iProcessToCase( ProcessList(idx) )
    ENDDO

    ! =======================================================================
    ! Per-case constants: flavour indices, chemical-potential selectors, the
    ! E3 sign that enters eta' (IMD evaluates its I_s at (E1,-E3)), and the
    ! deduplicated mass-pair slot. Replaces SetChemPotDistinctCases being
    ! called once per theta and the iABC_m2/iABC_m4 gather in the hot loop.
    ! =======================================================================
    ALLOCATE( TwoSided_Case(nDistinctCases), E3Sign_Case(nDistinctCases) )
    ALLOCATE( IsIMD_Case(nDistinctCases), iMassPair_Case(nDistinctCases) )
    ALLOCATE( iIdx2(nDistinctCases), iIdx4(nDistinctCases) )
    ALLOCATE( iChemSpecies2(nDistinctCases), iChemSpecies4(nDistinctCases) )
    ALLOCATE( ChemSign2(nDistinctCases),     ChemSign4(nDistinctCases)     )

    nMassPair = 0
    MPm2(:) = 0.0d0 ; MPm4(:) = 0.0d0 ; MPisIMD(:) = .FALSE.

    DO iCase = 1, nDistinctCases
      iTrueCase = ChosenToTrueCaseMap(iCase)

      ! Inverse muon and anti-muon decay share masses, A/B/C_IMD tables and
      ! kinematics (same mass-pair slot below); only the chemical-potential
      ! signs differ.
      IsIMD_Case   (iCase) = ( iTrueCase == i_IMD .OR. iTrueCase == i_IMDbar )
      TwoSided_Case(iCase) = IsIMD_Case(iCase)
      IF ( IsIMD_Case(iCase) ) THEN
        E3Sign_Case(iCase) = -1.0d0
        iIdx2(iCase) = 1          ! unused for IMD, kept in range
        iIdx4(iCase) = 1
      ELSE
        E3Sign_Case(iCase) = +1.0d0
        iIdx2(iCase) = iABC_m2( iTrueCase )
        iIdx4(iCase) = iABC_m4( iTrueCase )
      ENDIF

      ! --- deduplicate the (m2,m4,isIMD) kinematic key ---
      m2 = MassPairDistinct(iTrueCase,1)
      m4 = MassPairDistinct(iTrueCase,2)
      iMassPair_Case(iCase) = -1
      DO idx = 1, nMassPair
        ! NOTE the parentheses around .EQV.: it binds LOOSER than .AND., so
        ! without them this parses as (a .AND. b .AND. c) .EQV. d and every
        ! mass pair collapses onto slot 1.
        IF ( MPm2(idx) == m2 .AND. MPm4(idx) == m4 .AND. &
             ( MPisIMD(idx) .EQV. IsIMD_Case(iCase) ) ) THEN
          iMassPair_Case(iCase) = idx
          EXIT
        ENDIF
      ENDDO
      IF ( iMassPair_Case(iCase) < 0 ) THEN
        nMassPair = nMassPair + 1
        MPm2(nMassPair)    = m2
        MPm4(nMassPair)    = m4
        MPisIMD(nMassPair) = IsIMD_Case(iCase)
        iMassPair_Case(iCase) = nMassPair
      ENDIF

      SELECT CASE (iTrueCase)
      CASE (i_em_em)   ! (e-, e-)
        iChemSpecies2(iCase) = 1; ChemSign2(iCase) = +1.0d0
        iChemSpecies4(iCase) = 1; ChemSign4(iCase) = +1.0d0
      CASE (i_ep_ep)   ! (e+, e+)
        iChemSpecies2(iCase) = 1; ChemSign2(iCase) = -1.0d0
        iChemSpecies4(iCase) = 1; ChemSign4(iCase) = -1.0d0
      CASE (i_mm_mm)   ! (mu-, mu-)
        iChemSpecies2(iCase) = 2; ChemSign2(iCase) = +1.0d0
        iChemSpecies4(iCase) = 2; ChemSign4(iCase) = +1.0d0
      CASE (i_mp_mp)   ! (mu+, mu+)
        iChemSpecies2(iCase) = 2; ChemSign2(iCase) = -1.0d0
        iChemSpecies4(iCase) = 2; ChemSign4(iCase) = -1.0d0
      CASE (i_em_mm)   ! (e-, mu-)
        iChemSpecies2(iCase) = 1; ChemSign2(iCase) = +1.0d0
        iChemSpecies4(iCase) = 2; ChemSign4(iCase) = +1.0d0
      CASE (i_mm_em)   ! (mu-, e-)
        iChemSpecies2(iCase) = 2; ChemSign2(iCase) = +1.0d0
        iChemSpecies4(iCase) = 1; ChemSign4(iCase) = +1.0d0
      CASE (i_ep_mp)   ! (e+, mu+)
        iChemSpecies2(iCase) = 1; ChemSign2(iCase) = -1.0d0
        iChemSpecies4(iCase) = 2; ChemSign4(iCase) = -1.0d0
      CASE (i_mp_ep)   ! (mu+, e+)
        iChemSpecies2(iCase) = 2; ChemSign2(iCase) = -1.0d0
        iChemSpecies4(iCase) = 1; ChemSign4(iCase) = -1.0d0
      CASE (i_IMD)     ! (e-, mu-)
        iChemSpecies2(iCase) = 1; ChemSign2(iCase) = +1.0d0
        iChemSpecies4(iCase) = 2; ChemSign4(iCase) = +1.0d0
      CASE (i_IMDbar)  ! (e+, mu+)
        iChemSpecies2(iCase) = 1; ChemSign2(iCase) = -1.0d0
        iChemSpecies4(iCase) = 2; ChemSign4(iCase) = -1.0d0
      CASE DEFAULT
        WRITE(*,*) 'InitGeneralScatteringKernels: unhandled distinct case ', iTrueCase
        STOP
      END SELECT
    ENDDO

    ! =======================================================================
    ! Angular nodes, weights, coefficient and kinematic tables, built per
    ! (iE1, iE3, mass pair) on [xLo, xHi] in x = 1 - mu (see the module
    ! header). Nothing here depends on T or the chemical potentials.
    ! =======================================================================
    ALLOCATE( A_MP  (nE1, nE3, nThetaScattering, 2,      nMassPair) )
    ALLOCATE( B_MP  (nE1, nE3, nThetaScattering, 2,      nMassPair) )
    ALLOCATE( C_MP  (nE1, nE3, nThetaScattering, 3,      nMassPair) )
    ALLOCATE( PlW_MP(nE1, nE3, nThetaScattering, nL_Max, nMassPair) )

    ALLOCATE( Eminus_MP(nE1, nE3, nThetaScattering, nMassPair) )
    ALLOCATE( Eplus_MP (nE1, nE3, nThetaScattering, nMassPair) )
    ALLOCATE( Active_MP(nE1, nE3, nThetaScattering, nMassPair) )

    A_MP = 0.0d0 ; B_MP = 0.0d0 ; C_MP = 0.0d0 ; PlW_MP = 0.0d0
    Eminus_MP = 0.0d0 ; Eplus_MP = 0.0d0 ; Active_MP = .FALSE.

    ALLOCATE( sNode(nThetaScattering), sWght(nThetaScattering) )

    DO idx = 1, nMassPair

      m2 = MPm2(idx)
      m4 = MPm4(idx)
      Q  = 0.5d0 * (m4**2 - m2**2)                                      ! Guo Eq. 5

      ! --- map s(x) = sAlpha*ln(x) + sBeta*x for this mass pair
      IF ( MPisIMD(idx) ) THEN
        sAlpha = 1.0d0 ; sBeta = 1.0d0 / XcIMD ; xFloor = XFloor_IMD
      ELSE IF ( m2 == m4 ) THEN
        sAlpha = 1.0d0 ; sBeta = 0.0d0         ; xFloor = XFloor_Scat
      ELSE
        sAlpha = 0.0d0 ; sBeta = 1.0d0         ; xFloor = XFloor_Scat
      ENDIF

      DO iE3 = 1, nE3
        DO iE1 = 1, nE1

          E1 = E1_Grid(iE1)
          E3 = E3_Grid(iE3)

          ! --- angular range with phase space
          xHi = 2.0d0
          IF ( MPisIMD(idx) ) xHi = MIN( 2.0d0, Q / ( k0_IMD * E1 * E3 ) )  ! Guo Eq. 16
          IF ( xHi <= xFloor ) CYCLE                    ! no phase space: all inactive

          ! xLo is made symmetric under E1 <-> E3 and m2 <-> m4 (xHi, the map
          ! and xFloor already are). The (iE1,iE3) node set of this pair is
          ! then IDENTICAL to the (iE3,iE1) set of the same pair and of the
          ! mass-swapped pair, so the pointwise relations
          !   Phout_p(E1,E3) = Phout_p(E3,E1) * DB          (NES, NMS)
          !   Phout_r(E1,E3) = Phout_f(E3,E1) * DB          (LFE/LFC pairs)
          !   Phout_b(E1,E3) = Phout_a(E3,E1)               (IMD pairs)
          ! hold to round-off in the DISCRETE kernels (wlTestKernelHalving),
          ! when E1_Grid = E3_Grid as in thornado.
          xLo = AngularLowerLimit( E1, E3, m2, m4, MPisIMD(idx), xFloor, xHi )
          xLo = MIN( xLo, AngularLowerLimit( E3, E1, m2, m4, MPisIMD(idx), xFloor, xHi ) )
          IF ( .NOT. MPisIMD(idx) ) THEN
            xLo = MIN( xLo, AngularLowerLimit( E1, E3, m4, m2, .FALSE., xFloor, xHi ) )
            xLo = MIN( xLo, AngularLowerLimit( E3, E1, m4, m2, .FALSE., xFloor, xHi ) )
          ENDIF
          IF ( xLo >= xHi ) CYCLE

          CALL gauleg( MapS( xLo, sAlpha, sBeta ), MapS( xHi, sAlpha, sBeta ), &
                       sNode, sWght, nThetaScattering )

          DO iTh = 1, nThetaScattering

            x     = MapSInverse( sNode(iTh), sAlpha, sBeta )
            x     = MIN( MAX( x, xLo ), xHi )
            wx    = sWght(iTh) / ( sAlpha / x + sBeta )   ! dx = ds / (ds/dx)
            costh = 1.0d0 - x

            OneMinusMu   = x
            EEOneMinusMu = E1 * E3 * OneMinusMu
            k_val        = Q / EEOneMinusMu                                  ! Guo Eq. 5

            IF ( MPisIMD(idx) ) THEN

              ! --- IMD (Guo Eqs. 16-17). DeltaA from x, not mu, to keep
              ! --- precision when x is tiny.
              DeltaA   = SQRT( MAX( (E1 + E3)**2 - 2.0d0 * EEOneMinusMu, 0.0d0 ) ) ! Guo Eq. 17
              disc     = (1.0d0 - k_val)**2 - 2.0d0 * m2**2 / EEOneMinusMu   ! Guo Eq. 17
              IsActive = ( DeltaA >= 1.0d-10 ) .AND. ( k_val >= k0_IMD ) &
                         .AND. ( disc >= 0.0d0 )                             ! Guo Eq. 16
              IF ( IsActive ) THEN
                Em = 0.5d0 * ( (E3 + E1)*(k_val - 1.0d0) - DeltaA*SQRT(disc) ) ! Guo Eq. 17
                Ep = 0.5d0 * ( (E3 + E1)*(k_val - 1.0d0) + DeltaA*SQRT(disc) ) ! Guo Eq. 17
                Em = MAX(Em, m2)
                Ep = MAX(Ep, m2)
                IsActive = ( Em < Ep )
                Eminus_MP(iE1,iE3,iTh,idx) = Em
                Eplus_MP (iE1,iE3,iTh,idx) = Ep
              ENDIF

              IF ( IsActive ) THEN
                ! The "1" set enters with lam1 (processes 34/36), the "2" set
                ! (E1 <-> E3 swap) with lam2 (33/35). lam3 is always zero.
                Den = 16.0d0 * pi * DeltaA**5
                A_MP(iE1,iE3,iTh,1,idx) = A1_f( E1, -E3, costh )        / Den
                A_MP(iE1,iE3,iTh,2,idx) = A1_f( E3, -E1, costh )        / Den
                B_MP(iE1,iE3,iTh,1,idx) = B1_f( E1, -E3, costh, Q )     / Den
                B_MP(iE1,iE3,iTh,2,idx) = B1_f( E3, -E1, costh, Q )     / Den
                C_MP(iE1,iE3,iTh,1,idx) = C1_f( E1, -E3, costh, Q, m2 ) / Den
                C_MP(iE1,iE3,iTh,2,idx) = C1_f( E3, -E1, costh, Q, m2 ) / Den
              ENDIF

            ELSE

              ! --- 2 -> 2 scattering (Guo Eq. 4). disc is a square plus a
              ! --- positive term, so it can never be negative.
              Delta    = SQRT( (E1 - E3)**2 + 2.0d0 * EEOneMinusMu )           ! Guo Eq. 4
              IsActive = ( Delta >= 1.0d-10 )
              IF ( IsActive ) THEN
                disc = (1.0d0 + k_val)**2 + 2.0d0 * m2**2 / EEOneMinusMu      ! Guo Eq. 4
                Em   = 0.5d0 * ( (E3 - E1)*(1.0d0 + k_val) + Delta*SQRT(disc) ) ! Guo Eq. 4
                Eminus_MP(iE1,iE3,iTh,idx) = MAX(Em, m2)                      ! floor at m2

                Den = 16.0d0 * pi * Delta**5
                A_MP(iE1,iE3,iTh,1,idx) = A1_f(  E1,  E3, costh )         / Den
                A_MP(iE1,iE3,iTh,2,idx) = A1_f( -E3, -E1, costh )         / Den
                B_MP(iE1,iE3,iTh,1,idx) = B1_f(  E1,  E3, costh, Q )      / Den
                B_MP(iE1,iE3,iTh,2,idx) = B1_f( -E3, -E1, costh, Q )      / Den
                C_MP(iE1,iE3,iTh,1,idx) = C1_f(  E1,  E3, costh, Q, m2 )  / Den
                C_MP(iE1,iE3,iTh,2,idx) = C1_f( -E3, -E1, costh, Q, m2 )  / Den
                C_MP(iE1,iE3,iTh,3,idx) = C3_f(  E1,  E3, costh, m2, m4 ) / Den
              ENDIF

            ENDIF

            Active_MP(iE1,iE3,iTh,idx) = IsActive

            DO iL = 1, nL_Max
              PlW_MP(iE1,iE3,iTh,iL,idx) = LegendrePolynomial(iL-1, costh) * wx
            ENDDO

          ENDDO

        ENDDO
      ENDDO
    ENDDO

    DEALLOCATE( sNode, sWght )

    ! --- Preallocate collinear grids (generated ONCE, reused every call) ---
    ALLOCATE( xa_GLeg_ref(nE2_Lepton), wa_GLeg_ref(nE2_Lepton) )
    ALLOCATE( xa_GLag(nE2_Lepton),     wa_GLag(nE2_Lepton) )

    ! Canonical GL panel on [-1,1]; affine-mapped to [0,E1] at use time.
    CALL gauleg( -1.0d0, 1.0d0, xa_GLeg_ref, wa_GLeg_ref, nE2_Lepton )

    ! Raw Gauss-Laguerre (alpha=0) nodes/weights for int_0^inf e^{-x} g(x) dx.
    CALL gaulag( xa_GLag, wa_GLag, nE2_Lepton, 0.0d0 )

    ! =======================================================================
    ! Push the read-only tables to the device. Everything listed here is
    ! written once above and never again, so a single transfer at the end of
    ! Init is all that is needed.
    !
    ! OACC: the arrays are DECLARE CREATE'd, so ALLOCATE already made the
    ! device copy; ENTER DATA COPYIN would only bump the refcount and copy
    ! NOTHING. UPDATE DEVICE is what actually moves the data.
    ! =======================================================================
#if   defined(WEAKLIB_OACC)
    !$ACC UPDATE DEVICE( A_MP, B_MP, C_MP, PlW_MP,                            &
    !$ACC                Eminus_MP, Eplus_MP, Active_MP,                       &
    !$ACC                iMassPair_Case, IsIMD_Case, TwoSided_Case,            &
    !$ACC                E3Sign_Case, iIdx2, iIdx4,                            &
    !$ACC                lam1, lam2, lam3, iRowToCase,                         &
    !$ACC                iChemSpecies2, iChemSpecies4,                         &
    !$ACC                ChemSign2, ChemSign4, E1_Grid, E3_Grid,               &
    !$ACC                nThetaScattering, nDistinctCases, nProcessSel )
#elif defined(WEAKLIB_OMP_OL)
    !$OMP TARGET ENTER DATA MAP( always, to: A_MP, B_MP, C_MP, PlW_MP,         &
    !$OMP                         Eminus_MP, Eplus_MP, Active_MP,              &
    !$OMP                         iMassPair_Case, IsIMD_Case, TwoSided_Case,   &
    !$OMP                         E3Sign_Case, iIdx2, iIdx4,                   &
    !$OMP                         lam1, lam2, lam3, iRowToCase,                &
    !$OMP                         iChemSpecies2, iChemSpecies4,                &
    !$OMP                         ChemSign2, ChemSign4, E1_Grid, E3_Grid )
    !$OMP TARGET UPDATE TO( nThetaScattering, nDistinctCases, nProcessSel )
#endif

    WRITE(*,*) 'Working with', nDistinctCases, 'Distinct Scattering Cases'

  END SUBROUTINE InitGeneralScatteringKernels

  !================================================================================
  !  Angular-quadrature helpers (Init only, host only).
  !
  !  The map s(x) = alpha*ln(x) + beta*x, x = 1 - mu, is strictly increasing for
  !  alpha, beta >= 0 (not both zero).
  !================================================================================
  PURE FUNCTION MapS( x, alpha, beta ) RESULT( sv )
    REAL(DP), INTENT(IN) :: x, alpha, beta
    REAL(DP) :: sv
    sv = beta * x
    IF ( alpha > 0.0d0 ) sv = sv + alpha * LOG( x )
  END FUNCTION MapS

  ! Inverse of MapS. For the log-linear case, Newton on t = ln(x):
  ! f(t) = alpha*t + beta*e^t - s is increasing and convex, so Newton started
  ! to the right of the root (t0 = s/alpha >= root) converges monotonically.
  PURE FUNCTION MapSInverse( sv, alpha, beta ) RESULT( x )
    REAL(DP), INTENT(IN) :: sv, alpha, beta
    REAL(DP) :: x, t, dt
    INTEGER  :: it
    IF ( alpha <= 0.0d0 ) THEN
      x = sv / beta
    ELSE IF ( beta <= 0.0d0 ) THEN
      x = EXP( sv / alpha )
    ELSE
      t = sv / alpha
      DO it = 1, 200
        dt = ( alpha * t + beta * EXP(t) - sv ) / ( alpha + beta * EXP(t) )
        t  = t - dt
        IF ( ABS(dt) <= 1.0d-15 * MAX( 1.0d0, ABS(t) ) ) EXIT
      ENDDO
      x = EXP( t )
    ENDIF
  END FUNCTION MapSInverse

  ! Kinematic lower bound E_minus on the lepton-2 energy at x = 1 - mu
  ! (Guo Eq. 4 for 2->2, Eq. 17 for IMD), floored at m2. HUGE where an IMD
  ! configuration has no phase space.
  PURE FUNCTION KinematicEminus( E1, E3, x, m2, m4, IsIMD ) RESULT( Em )
    REAL(DP), INTENT(IN) :: E1, E3, x, m2, m4
    LOGICAL , INTENT(IN) :: IsIMD
    REAL(DP) :: Em, Q, EEx, k_val, disc, D
    Q     = 0.5d0 * ( m4**2 - m2**2 )
    EEx   = E1 * E3 * x
    k_val = Q / EEx
    IF ( IsIMD ) THEN
      D    = SQRT( MAX( (E1 + E3)**2 - 2.0d0 * EEx, 0.0d0 ) )
      disc = (1.0d0 - k_val)**2 - 2.0d0 * m2**2 / EEx
      IF ( k_val < k0_IMD .OR. disc < 0.0d0 ) THEN
        Em = HUGE( 1.0d0 )
        RETURN
      ENDIF
      Em = 0.5d0 * ( (E3 + E1)*(k_val - 1.0d0) - D*SQRT(disc) )
    ELSE
      D    = SQRT( (E1 - E3)**2 + 2.0d0 * EEx )
      disc = (1.0d0 + k_val)**2 + 2.0d0 * m2**2 / EEx
      Em   = 0.5d0 * ( (E3 - E1)*(1.0d0 + k_val) + D*SQRT(disc) )
    ENDIF
    Em = MAX( Em, m2 )
  END FUNCTION KinematicEminus

  ! Smallest x in [xFloor, xHi] at which E_minus(x) <= EmCap. Below it the
  ! I_s integrals start above EmCap and are negligible for any T, mu of
  ! interest. E_minus decreases with x in the small-x region (it grows like
  ! 1/x for Q /= 0 and like 1/sqrt(x) off the diagonal for Q = 0), so the
  ! crossing is found by bisection in ln(x). Returns xFloor if E_minus(xFloor)
  ! is already below the cap, and xHi if there is no admissible x at all.
  PURE FUNCTION AngularLowerLimit( E1, E3, m2, m4, IsIMD, xFloor, xHi ) RESULT( xLo )
    REAL(DP), INTENT(IN) :: E1, E3, m2, m4, xFloor, xHi
    LOGICAL , INTENT(IN) :: IsIMD
    REAL(DP) :: xLo, a, b, c
    INTEGER  :: it
    IF ( KinematicEminus( E1, E3, xFloor, m2, m4, IsIMD ) <= EmCap ) THEN
      xLo = xFloor
      RETURN
    ENDIF
    IF ( KinematicEminus( E1, E3, xHi, m2, m4, IsIMD ) > EmCap ) THEN
      xLo = xHi
      RETURN
    ENDIF
    a = LOG( xFloor )   ! E_minus > EmCap
    b = LOG( xHi )      ! E_minus <= EmCap
    DO it = 1, 100
      c = 0.5d0 * ( a + b )
      IF ( KinematicEminus( E1, E3, EXP(c), m2, m4, IsIMD ) > EmCap ) THEN
        a = c
      ELSE
        b = c
      ENDIF
      IF ( b - a < 1.0d-12 ) EXIT
    ENDDO
    xLo = EXP( a )
  END FUNCTION AngularLowerLimit

  SUBROUTINE FinalizeGeneralScatteringKernels()

    ! Release the device copies BEFORE the host arrays go away.
    ! OACC needs nothing here: DEALLOCATE of a DECLARE CREATE'd array frees
    ! the device copy as well.
#if defined(WEAKLIB_OMP_OL)
    !$OMP TARGET EXIT DATA MAP( delete: A_MP, B_MP, C_MP, PlW_MP,              &
    !$OMP                         Eminus_MP, Eplus_MP, Active_MP,              &
    !$OMP                         iMassPair_Case, IsIMD_Case, TwoSided_Case,   &
    !$OMP                         E3Sign_Case, iIdx2, iIdx4,                   &
    !$OMP                         lam1, lam2, lam3, iRowToCase,                &
    !$OMP                         iChemSpecies2, iChemSpecies4,                &
    !$OMP                         ChemSign2, ChemSign4, E1_Grid, E3_Grid )
#endif

    DEALLOCATE( lam1, lam2, lam3 )
    DEALLOCATE( A_MP, B_MP, C_MP, PlW_MP )
    DEALLOCATE( Eminus_MP, Eplus_MP, Active_MP )
    DEALLOCATE( iMassPair_Case, IsIMD_Case, TwoSided_Case, E3Sign_Case )
    DEALLOCATE( iIdx2, iIdx4 )
    DEALLOCATE( iChemSpecies2, iChemSpecies4, ChemSign2, ChemSign4 )
    DEALLOCATE( E1_Grid, E3_Grid )
    DEALLOCATE( xa_GLeg_ref, wa_GLeg_ref, xa_GLag, wa_GLag )
    DEALLOCATE( ChosenToTrueCaseMap )
    DEALLOCATE( ProcessList, iRowToCase )

  END SUBROUTINE FinalizeGeneralScatteringKernels


  SUBROUTINE SelectiDistinctFromProcessIndex(ProcessIndex, iDistinct)

    ! Same structure as SelectLambdaFromProcessIndex

    INTEGER, INTENT(IN)  :: ProcessIndex
    INTEGER, INTENT(OUT) :: iDistinct

    SELECT CASE(ProcessIndex)
    ! ============================================================
    ! ELASTIC NC SCATTERING ONLY (same in/out flavor): 1..24
    ! ============================================================

    ! 1) nu_e + e-
    CASE(1); iDistinct = i_em_em
    ! 2) nu_e + e+
    CASE(2); iDistinct = i_ep_ep
    ! 3) nu_bar_e + e-
    CASE(3); iDistinct = i_em_em
    ! 4) nu_bar_e + e+
    CASE(4); iDistinct = i_ep_ep

    ! 5) nu_mu + e-
    CASE(5); iDistinct = i_em_em
    ! 6) nu_mu + e+
    CASE(6); iDistinct = i_ep_ep
    ! 7) nu_bar_mu + e-
    CASE(7); iDistinct = i_em_em
    ! 8) nu_bar_mu + e+
    CASE(8); iDistinct = i_ep_ep

    ! 9) nu_tau + e-
    CASE(9); iDistinct = i_em_em
    ! 10) nu_tau + e+
    CASE(10); iDistinct = i_ep_ep
    ! 11) nu_bar_tau + e-
    CASE(11); iDistinct = i_em_em
    ! 12) nu_bar_tau + e+
    CASE(12); iDistinct = i_ep_ep

    ! 13) nu_e + mu-
    CASE(13); iDistinct = i_mm_mm
    ! 14) nu_e + mu+
    CASE(14); iDistinct = i_mp_mp
    ! 15) nu_bar_e + mu-
    CASE(15); iDistinct = i_mm_mm
    ! 16) nu_bar_e + mu+
    CASE(16); iDistinct = i_mp_mp

    ! 17) nu_mu + mu-
    CASE(17); iDistinct = i_mm_mm
    ! 18) nu_mu + mu+
    CASE(18); iDistinct = i_mp_mp
    ! 19) nu_bar_mu + mu-
    CASE(19); iDistinct = i_mm_mm
    ! 20) nu_bar_mu + mu+
    CASE(20); iDistinct = i_mp_mp

    ! 21) nu_tau + mu-
    CASE(21); iDistinct = i_mm_mm
    ! 22) nu_tau + mu+
    CASE(22); iDistinct = i_mp_mp
    ! 23) nu_bar_tau + mu-
    CASE(23); iDistinct = i_mm_mm
    ! 24) nu_bar_tau + mu+
    CASE(24); iDistinct = i_mp_mp

    ! ============================================================
    ! INELASTIC FLAVOR EXCHANGE (LFE) & CONVERSION (LFC): 25..32
    ! ============================================================
    
    ! LFE (t-channel W-exchange) 
    ! 25) nu_mu + e- -> nu_e + mu-
    CASE(25); iDistinct = i_em_mm
    ! 26) nu_e + mu- -> nu_mu + e-
    CASE(26); iDistinct = i_mm_em
    ! 27) nu_bar_mu + e+ -> nu_bar_e + mu+
    CASE(27); iDistinct = i_ep_mp
    ! 28) nu_bar_e + mu+ -> nu_bar_mu + e+
    CASE(28); iDistinct = i_mp_ep

    ! LFC (s-channel W-exchange)
    ! 29) nu_mu + mu+ -> nu_e + e+
    CASE(29); iDistinct = i_mp_ep
    ! 30) nu_e + e+ -> nu_mu + mu+
    CASE(30); iDistinct = i_ep_mp
    ! 31) nu_bar_mu + mu- -> nu_bar_e + e-
    CASE(31); iDistinct = i_mm_em
    ! 32) nu_bar_e + e- -> nu_bar_mu + mu-
    CASE(32); iDistinct = i_em_mm

    CASE(33); iDistinct = i_IMD
    CASE(34); iDistinct = i_IMD ! Same as 33 but need to switch nue_bar and nu_mu

    ! 35) nu_e + e+ + nu_bar_mu -> mu+ ; 36) nu_bar_mu + e+ + nu_e -> mu+
    CASE(35, 36); iDistinct = i_IMDbar

    CASE DEFAULT
      WRITE(*,*) 'Error: Unrecognized ProcessIndex: ', ProcessIndex
      STOP
    END SELECT

  END SUBROUTINE SelectiDistinctFromProcessIndex

  !================================================================================
  !
  !  CalculateAllPhout
  !
  !  Per-point (iE1,iE3,T,Mu_e,Mu_mu) kernel. Device-callable (ROUTINE SEQ /
  !  DECLARE TARGET), so it can sit inside a caller's offloaded loop over
  !  (iX, iE3, iE1) the same way Opacity_CC_2D does for Em/Ab. Caller notes:
  !    - the private Phout the caller passes should have PARAMETER extents
  !      (e.g. (nProcessTotal, nL_Max)), or it becomes a per-thread heap
  !      allocation on the device;
  !    - make iE1 the fastest-varying loop index so the table reads
  !      (A_MP(iE1,...), Eminus_MP(iE1,...)) coalesce;
  !    - nL > nL_Max is only checked on the host; on device nL is clamped.
  !  For bulk work prefer CalculateAllPhout_3D: same arithmetic (results are
  !  bit-identical), but one tuned kernel with coalesced Phout stores.
  !
  !================================================================================
  SUBROUTINE CalculateAllPhout( iE1, iE3, E1, E3, T, Mu_e, Mu_mu, Phout, nL )
#if   defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#elif defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#endif

    INTEGER , INTENT(IN)  :: iE1, iE3
    REAL(DP), INTENT(IN)  :: E1, E3, T, Mu_e, Mu_mu
    REAL(DP), INTENT(OUT) :: Phout(nProcessSel,nL)
    INTEGER , INTENT(IN)  :: nL

    REAL(DP) :: E3eff, xMu_l2, xMu_l4, xi, y, fac
    REAL(DP) :: I0, I1, I2, I0p, I1p, I2p, R1, R2, R3
    REAL(DP) :: MuSel(2)
    INTEGER  :: iTh, iL, idx, iCase, iMP, nLeff

    ! Integrated components for the distinct cases ONLY. Extents are PARAMETERs,
    ! never the nL dummy: a runtime-extent automatic array becomes a per-thread
    ! heap allocation once this is offloaded.
    REAL(DP) :: Int_R1(nDistinctCasesMax, nL_Max)
    REAL(DP) :: Int_R2(nDistinctCasesMax, nL_Max)
    REAL(DP) :: Int_R3(nDistinctCasesMax, nL_Max)
    REAL(DP) :: eta_c (nDistinctCasesMax)
    REAL(DP) :: etap_c(nDistinctCasesMax)
    REAL(DP) :: fg_c  (nDistinctCasesMax)
    LOGICAL  :: degen_c(nDistinctCasesMax)

#if !defined(WEAKLIB_OACC) && !defined(WEAKLIB_OMP_OL)
    IF ( nL > nL_Max ) THEN
      WRITE(*,*) 'CalculateAllPhout: nL =', nL, ' exceeds nL_Max =', nL_Max
      STOP
    ENDIF
#endif
    nLeff = MIN( nL, nL_Max )

    Phout  = 0.0d0
    Int_R1 = 0.0d0
    Int_R2 = 0.0d0
    Int_R3 = 0.0d0

    ! IF ( ABS(E3 - E1) < 1.0d-10) THEN
    !   IF (nL == 1)  THEN
    !     CALL CalculateAllPhi0_collinear( iE1, E1, T, Mu_e, Mu_mu, Phout(:,1) )
    !   ELSE
    !     CALL CalculateAllPhoutGeneral_collinear( iE1, E1, T, Mu_e, Mu_mu, Phout, nL )
    !   ENDIF
    !   RETURN
    ! ENDIF

    ! -------------------------------------------------------------------
    ! 1. Per-case setup, hoisted above the theta loop: eta, eta', f_gamma
    !    and the xi -> 0 test are theta-independent. Same arithmetic as
    !    Compute_Is_Integrals (method 3), so results are unchanged.
    ! -------------------------------------------------------------------
    MuSel(1) = Mu_e
    MuSel(2) = Mu_mu

    DO iCase = 1, nDistinctCases
      xMu_l2 = ChemSign2(iCase) * MuSel( iChemSpecies2(iCase) )
      xMu_l4 = ChemSign4(iCase) * MuSel( iChemSpecies4(iCase) )

      ! IMD evaluates its I_s at (E1,-E3); scattering at (E1,E3).
      E3eff = E3Sign_Case(iCase) * E3

      eta_c (iCase) = xMu_l2 / T
      etap_c(iCase) = ( E3eff - E1 + xMu_l4 ) / T

      xi             = etap_c(iCase) - eta_c(iCase)
      degen_c(iCase) = ( ABS(xi) < xi_tol )
      fg_c  (iCase)  = f_gamma( xi )
    ENDDO

    ! -------------------------------------------------------------------
    ! 2. Angular integration over the DISTINCT cases ONLY
    ! -------------------------------------------------------------------
    DO iTh = 1, nThetaScattering
      DO iCase = 1, nDistinctCases

        iMP = iMassPair_Case(iCase)
        IF ( .NOT. Active_MP(iE1,iE3,iTh,iMP) ) CYCLE

        y = Eminus_MP(iE1,iE3,iTh,iMP) / T
        CALL Compute_Is_Integrals_FD( y, eta_c(iCase), etap_c(iCase), T, &
                                      fg_c(iCase), degen_c(iCase), I0, I1, I2 )

        IF ( TwoSided_Case(iCase) ) THEN
          ! I_s = I_s(E_minus) - I_s(E_plus)   (Guo Eq. 18)
          y = Eplus_MP(iE1,iE3,iTh,iMP) / T
          CALL Compute_Is_Integrals_FD( y, eta_c(iCase), etap_c(iCase), T, &
                                        fg_c(iCase), degen_c(iCase), I0p, I1p, I2p )
          I0 = I0 - I0p
          I1 = I1 - I1p
          I2 = I2 - I2p
        ENDIF

        ! Guo Eq. 3; the 1/(16*pi*Delta^5) normalization is baked into A/B/C
        ! (C(3) is zero for IMD).
        R1 = A_MP(iE1,iE3,iTh,1,iMP)*I2 &
           + B_MP(iE1,iE3,iTh,1,iMP)*I1 &
           + C_MP(iE1,iE3,iTh,1,iMP)*I0
        R2 = A_MP(iE1,iE3,iTh,2,iMP)*I2 &
           + B_MP(iE1,iE3,iTh,2,iMP)*I1 &
           + C_MP(iE1,iE3,iTh,2,iMP)*I0
        R3 = C_MP(iE1,iE3,iTh,3,iMP)*I0

        ! P_{iL-1}(mu) * w, tabulated at Init for this pair's nodes
        DO iL = 1, nLeff
          Int_R1(iCase,iL) = Int_R1(iCase,iL) + R1 * PlW_MP(iE1,iE3,iTh,iL,iMP)
          Int_R2(iCase,iL) = Int_R2(iCase,iL) + R2 * PlW_MP(iE1,iE3,iTh,iL,iMP)
          Int_R3(iCase,iL) = Int_R3(iCase,iL) + R3 * PlW_MP(iE1,iE3,iTh,iL,iMP)
        ENDDO

      ENDDO
    ENDDO

    ! -------------------------------------------------------------------
    ! 3. Map the integrated distinct cases back to the full iProcess array
    ! -------------------------------------------------------------------
    DO iL = 1, nLeff
      ! factors of 1/2 and 3/2 to match Bruenn (i.e. ((iL-1) + 0.5d0)),
      ! folded in here so there is no second pass over Phout.
      fac = conv_fac * ((iL-1) + 0.5d0)

      DO idx = 1, nProcessSel
        iCase = iRowToCase(idx)

        ! Apply the Lambdas to the integrated values once
        Phout(idx, iL) = fac * ( lam1(idx) * Int_R1(iCase, iL) + &
                                 lam2(idx) * Int_R2(iCase, iL) + &
                                 lam3(idx) * Int_R3(iCase, iL) )
      END DO
    END DO

  END SUBROUTINE CalculateAllPhout

  !================================================================================
  !
  !  CalculateAllPhout_3D
  !
  !  Batched, offloadable form of CalculateAllPhout: runs the whole
  !  (nX, nE3, nE1) nest as ONE kernel. A single CalculateAllPhout call is far
  !  too little work to amortize a kernel launch (~nTheta*nCase*6 rational
  !  function evaluations), so the parallelism has to come from this nest.
  !
  !  Nest order is iX -> iE3 -> iE1 deliberately. Under COLLAPSE(3) the LAST
  !  loop varies fastest, so iE1 strides across adjacent threads; with nE1
  !  leading in Phout the stores are stride-1, and the existing
  !  A_Case(iE1,iE3,iTh,...) layout is already coalesced for the same reason.
  !
  !  T, Mu_e, Mu_mu and Phout are expected to be ALREADY DEVICE-RESIDENT in the
  !  caller's data region: they are declared PRESENT, not mapped here.
  !
  !  Divisions by T are kept as divisions (rather than multiplying by a hoisted
  !  reciprocal) so that this routine is bit-for-bit identical to the scalar
  !  CalculateAllPhout. That makes the GPU-vs-CPU check an exact-equality test.
  !  The divides are a small fraction of the cost next to the Fukushima
  !  evaluations; revisit only if profiling says otherwise.
  !
  !================================================================================
  SUBROUTINE CalculateAllPhout_3D( nX, nE1, nE3, nProc, nL, T, Mu_e, Mu_mu, Phout )

    INTEGER , INTENT(IN)  :: nX, nE1, nE3, nProc, nL
    REAL(DP), INTENT(IN)  :: T(nX), Mu_e(nX), Mu_mu(nX)
    REAL(DP), INTENT(OUT) :: Phout(nE1, nE3, nProc, nL, nX)

    REAL(DP) :: E1, E3, E3eff, xMu_l2, xMu_l4, xi, y, fac
    REAL(DP) :: I0, I1, I2, I0p, I1p, I2p, R1, R2, R3
    REAL(DP) :: MuSel(2)
    INTEGER  :: iX, iE1, iE3, iTh, iL, iCase, idx, iMP
    INTEGER  :: nTh, nCase

    ! --- Thread-private accumulators. Every extent is a PARAMETER: a
    ! --- runtime-extent automatic array becomes a per-thread heap allocation
    ! --- on a device. ~0.6 kB per thread at nL_Max = 2 and 9 cases.
    REAL(DP) :: Int_R1(nDistinctCasesMax, nL_Max)
    REAL(DP) :: Int_R2(nDistinctCasesMax, nL_Max)
    REAL(DP) :: Int_R3(nDistinctCasesMax, nL_Max)
    REAL(DP) :: eta_c (nDistinctCasesMax)
    REAL(DP) :: etap_c(nDistinctCasesMax)
    REAL(DP) :: fg_c  (nDistinctCasesMax)
    LOGICAL  :: degen_c(nDistinctCasesMax)

    ! ------------------------------------------------------------------
    ! All validation happens on the host: no I/O or STOP inside the region.
    ! ------------------------------------------------------------------
    IF ( nL > nL_Max ) THEN
      WRITE(*,*) 'CalculateAllPhout_3D: nL =', nL, ' exceeds nL_Max =', nL_Max
      STOP
    ENDIF
    IF ( nE1 /= nE1_Grid .OR. nE3 /= nE3_Grid ) THEN
      WRITE(*,*) 'CalculateAllPhout_3D: (nE1,nE3) = (', nE1, ',', nE3, &
                 ') does not match the grid Init was called with: (', &
                 nE1_Grid, ',', nE3_Grid, ')'
      STOP
    ENDIF
    IF ( nProc /= nProcessSel ) THEN
      WRITE(*,*) 'CalculateAllPhout_3D: nProc =', nProc, ' expected ', nProcessSel
      STOP
    ENDIF

    ! Module scalars copied into locals so they are captured by value. A module
    ! scalar that was never mapped reads as 0 on device, which would silently
    ! skip the iTh loop and return all zeros with no error.
    nTh   = nThetaScattering
    nCase = nDistinctCases

#if   defined(WEAKLIB_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3)                                &
    !$ACC PRESENT( A_MP, B_MP, C_MP, PlW_MP,                                  &
    !$ACC          Eminus_MP, Eplus_MP, Active_MP, iMassPair_Case,             &
    !$ACC          IsIMD_Case, TwoSided_Case, E3Sign_Case, iIdx2, iIdx4,       &
    !$ACC          lam1, lam2, lam3, iRowToCase, iChemSpecies2,                &
    !$ACC          iChemSpecies4, ChemSign2, ChemSign4, E1_Grid, E3_Grid,      &
    !$ACC          T, Mu_e, Mu_mu, Phout )                                     &
    !$ACC PRIVATE( Int_R1, Int_R2, Int_R3, eta_c, etap_c, fg_c, degen_c, MuSel )
#elif defined(WEAKLIB_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3)                      &
    !$OMP PRIVATE( Int_R1, Int_R2, Int_R3, eta_c, etap_c, fg_c, degen_c, MuSel,&
    !$OMP          E1, E3, E3eff, xMu_l2, xMu_l4, xi, y, fac,                  &
    !$OMP          I0, I1, I2, I0p, I1p, I2p, R1, R2, R3, iTh, iL, iCase, idx, &
    !$OMP          iMP )
#endif
    DO iX  = 1, nX
    DO iE3 = 1, nE3
    DO iE1 = 1, nE1

      E1 = E1_Grid(iE1)
      E3 = E3_Grid(iE3)

      MuSel(1) = Mu_e (iX)
      MuSel(2) = Mu_mu(iX)

      ! ----------------------------------------------------------------
      ! Per-case setup, hoisted ABOVE the theta loop: eta, eta', f_gamma and
      ! the xi -> 0 test are all theta-INDEPENDENT. This lifts an EXP and two
      ! divides out of the innermost loop.
      ! ----------------------------------------------------------------
#if defined(WEAKLIB_OACC)
      !$ACC LOOP SEQ
#endif
      DO iCase = 1, nCase
        xMu_l2 = ChemSign2(iCase) * MuSel( iChemSpecies2(iCase) )
        xMu_l4 = ChemSign4(iCase) * MuSel( iChemSpecies4(iCase) )

        ! IMD evaluates its I_s at (E1,-E3); scattering at (E1,E3).
        E3eff = E3Sign_Case(iCase) * E3

        eta_c (iCase) = xMu_l2 / T(iX)
        etap_c(iCase) = ( E3eff - E1 + xMu_l4 ) / T(iX)

        xi             = etap_c(iCase) - eta_c(iCase)
        degen_c(iCase) = ( ABS(xi) < xi_tol )
        fg_c  (iCase)  = f_gamma( xi )
      ENDDO

      Int_R1 = 0.0d0
      Int_R2 = 0.0d0
      Int_R3 = 0.0d0

      ! ----------------------------------------------------------------
      ! Angular integral. Everything kinematic is a table lookup.
      ! ----------------------------------------------------------------
#if defined(WEAKLIB_OACC)
      !$ACC LOOP SEQ
#endif
      DO iTh = 1, nTh
#if defined(WEAKLIB_OACC)
        !$ACC LOOP SEQ
#endif
        DO iCase = 1, nCase

          iMP = iMassPair_Case(iCase)
          IF ( .NOT. Active_MP(iE1,iE3,iTh,iMP) ) CYCLE

          y = Eminus_MP(iE1,iE3,iTh,iMP) / T(iX)
          CALL Compute_Is_Integrals_FD( y, eta_c(iCase), etap_c(iCase), T(iX), &
                                        fg_c(iCase), degen_c(iCase), I0, I1, I2 )

          IF ( TwoSided_Case(iCase) ) THEN
            ! I_s = I_s(E_minus) - I_s(E_plus)   (Guo Eq. 18)
            y = Eplus_MP(iE1,iE3,iTh,iMP) / T(iX)
            CALL Compute_Is_Integrals_FD( y, eta_c(iCase), etap_c(iCase), T(iX), &
                                          fg_c(iCase), degen_c(iCase), I0p, I1p, I2p )
            I0 = I0 - I0p
            I1 = I1 - I1p
            I2 = I2 - I2p
          ENDIF

          ! Guo Eq. 3; the 1/(16*pi*Delta^5) normalization is baked into A/B/C
          ! (C(3) is zero for IMD), so there is no IMD branch.
          R1 = A_MP(iE1,iE3,iTh,1,iMP)*I2 &
             + B_MP(iE1,iE3,iTh,1,iMP)*I1 &
             + C_MP(iE1,iE3,iTh,1,iMP)*I0
          R2 = A_MP(iE1,iE3,iTh,2,iMP)*I2 &
             + B_MP(iE1,iE3,iTh,2,iMP)*I1 &
             + C_MP(iE1,iE3,iTh,2,iMP)*I0
          R3 = C_MP(iE1,iE3,iTh,3,iMP)*I0

#if defined(WEAKLIB_OACC)
          !$ACC LOOP SEQ
#endif
          DO iL = 1, nL
            Int_R1(iCase,iL) = Int_R1(iCase,iL) + R1 * PlW_MP(iE1,iE3,iTh,iL,iMP)
            Int_R2(iCase,iL) = Int_R2(iCase,iL) + R2 * PlW_MP(iE1,iE3,iTh,iL,iMP)
            Int_R3(iCase,iL) = Int_R3(iCase,iL) + R3 * PlW_MP(iE1,iE3,iTh,iL,iMP)
          ENDDO

        ENDDO
      ENDDO

      ! ----------------------------------------------------------------
      ! Map the integrated distinct cases onto the full process list.
      ! ----------------------------------------------------------------
#if defined(WEAKLIB_OACC)
      !$ACC LOOP SEQ
#endif
      DO iL = 1, nL
        ! factors of 1/2 and 3/2 to match Bruenn (i.e. ((iL-1) + 0.5d0))
        fac = conv_fac * ((iL-1) + 0.5d0)

#if defined(WEAKLIB_OACC)
        !$ACC LOOP SEQ
#endif
        DO idx = 1, nProc
          iCase = iRowToCase(idx)

          Phout(iE1,iE3,idx,iL,iX) = fac * ( lam1(idx) * Int_R1(iCase,iL)   &
                                           + lam2(idx) * Int_R2(iCase,iL)   &
                                           + lam3(idx) * Int_R3(iCase,iL) )
        ENDDO
      ENDDO

    ENDDO
    ENDDO
    ENDDO
#if defined(WEAKLIB_OACC)
    !$ACC END PARALLEL LOOP
#elif defined(WEAKLIB_OMP_OL)
    !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
#endif

  END SUBROUTINE CalculateAllPhout_3D

  !================================================================================
  !
  !  DetailedBalanceFactor
  !
  !  Phin = Phout * DetailedBalanceFactor( iRow, ... ) for row iRow of Phout.
  !
  !    scattering rows (Guo Eq. 13):  exp( (E3 - E1 - Dmu) / T )
  !      Phin is the kernel for nu(E3) + l4 -> nu(E1) + l2.
  !    inverse (anti)muon decay rows (Guo Eq. 21):  exp( -(E1 + E3 + Dmu) / T )
  !      Phout is the absorption kernel R_A (nu1(E1) + l2 + nu3(E3) -> l4),
  !      Phin the decay kernel R_D (l4 -> l2 + nu1(E1) + nu3(E3)), i.e. the
  !      emissivity kernel.
  !
  !  with Dmu = mu_l2 - mu_l4, using the signed (anti)lepton chemical potentials
  !  of that row. The relation is exact, not an approximation: the two kernels
  !  share matrix element and phase space and differ only by f4(1-f2) vs
  !  f2(1-f4), whose ratio is exactly this exponential for Fermi-Dirac leptons
  !  at common T. Using it (rather than integrating R_D separately) also
  !  guarantees the discrete collision term vanishes in equilibrium.
  !
  !  Device-callable; uses only device-declared module data.
  !
  !================================================================================
  FUNCTION DetailedBalanceFactor( iRow, E1, E3, T, Mu_e, Mu_mu ) RESULT( fac )
#if   defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#elif defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#endif

    INTEGER , INTENT(IN) :: iRow
    REAL(DP), INTENT(IN) :: E1, E3, T, Mu_e, Mu_mu
    REAL(DP) :: fac

    REAL(DP) :: MuSel(2), Delta_mu, x
    INTEGER  :: iCase

    MuSel(1) = Mu_e
    MuSel(2) = Mu_mu

    iCase    = iRowToCase(iRow)
    Delta_mu = ChemSign2(iCase) * MuSel( iChemSpecies2(iCase) ) &
             - ChemSign4(iCase) * MuSel( iChemSpecies4(iCase) )

    IF ( IsIMD_Case(iCase) ) THEN
      x = -( E1 + E3 + Delta_mu ) / T
    ELSE
      x =  ( E3 - E1 - Delta_mu ) / T
    ENDIF

    fac = EXP( MIN( x, DetBal_MaxExponent ) )

  END FUNCTION DetailedBalanceFactor

  SUBROUTINE CalculateAllRout( E1, E3, costh, T, Mu_e, Mu_mu, iE1, iE3, iTh, Rout )

    REAL(DP), INTENT(IN)  :: E1, E3, costh, T, Mu_e, Mu_mu
    INTEGER , INTENT(IN)  :: iE1, iE3, iTh
    ! Indexed by selection row, matching Phout and CalculateAllRoutIntegrated:
    ! use ProcessRowIndex(iProcess) to find the row of an absolute process index.
    REAL(DP), INTENT(OUT) :: Rout(nProcessSel)

    ! Temporary arrays to hold the distinct R components for this theta
    REAL(DP) :: R1_th(nDistinctCasesMax)
    REAL(DP) :: R2_th(nDistinctCasesMax)
    REAL(DP) :: R3_th(nDistinctCasesMax)
    INTEGER  :: idx, iCase

    ! 1. Get the raw distinct R components for this specific angle/energy.
    ! costh is NOT used: the angle enters only through the iTh-indexed tables,
    ! and node iTh sits at a different cos(theta) for every (iE1,iE3,mass pair)
    ! (see the module header). costh stays in the signature for interface
    ! stability only.
    CALL CalculateAllR1R2R3( E1, E3, T, Mu_e, Mu_mu, &
                             iE1, iE3, iTh, R1_th, R2_th, R3_th )

    ! 2. Safely initialize Rout
    Rout = 0.0d0

    ! 3. Map the distinct cases to the full iProcess array using Lambdas
    DO idx = 1, nProcessSel
      ! Fetch the compact case index for this process
      iCase = iRowToCase(idx)

      ! Combine them
      Rout(idx) = lam1(idx) * R1_th(iCase) + &
                  lam2(idx) * R2_th(iCase) + &
                  lam3(idx) * R3_th(iCase)
    END DO

  END SUBROUTINE CalculateAllRout

  SUBROUTINE CalculateAllRoutIntegrated( iE1, iE3, E1, E3, T, Mu_e, Mu_mu, Rout )

    REAL(DP), INTENT(IN)  :: E1, E3, T, Mu_e, Mu_mu
    INTEGER , INTENT(IN)  :: iE1, iE3
    ! Compact indexing, as in CalculateAllRout.
    REAL(DP), INTENT(OUT) :: Rout(nProcessSel)

    ! Temporary arrays to hold the distinct R components for this theta
    REAL(DP) :: R1_th(nDistinctCasesMax)
    REAL(DP) :: R2_th(nDistinctCasesMax)
    REAL(DP) :: R3_th(nDistinctCasesMax)
    INTEGER  :: idx, iTh, iCase
    REAL(DP) :: w

    ! Arrays to hold the integrated components for the distinct cases ONLY
    REAL(DP) :: Int_R1(nDistinctCasesMax)
    REAL(DP) :: Int_R2(nDistinctCasesMax)
    REAL(DP) :: Int_R3(nDistinctCasesMax)

    Int_R1(:) = 0.0d0
    Int_R2(:) = 0.0d0
    Int_R3(:) = 0.0d0
    DO iTh=1, nThetaScattering

      CALL CalculateAllR1R2R3( E1, E3, T, Mu_e, Mu_mu, &
                               iE1, iE3, iTh, R1_th, R2_th, R3_th )

        ! Only integrate the unique cases we actually need
        ! Weight of this case's mass-pair node (P_0 = 1)
        DO iCase = 1, nDistinctCases
          w = PlW_MP(iE1,iE3,iTh,1,iMassPair_Case(iCase))
          Int_R1(iCase) = Int_R1(iCase) + R1_th(iCase) * w
          Int_R2(iCase) = Int_R2(iCase) + R2_th(iCase) * w
          Int_R3(iCase) = Int_R3(iCase) + R3_th(iCase) * w
        END DO

    ENDDO
    ! 2. Safely initialize Rout
    Rout = 0.0d0

    ! -------------------------------------------------------------------
    ! 2. Map the integrated distinct cases back to the full iProcess array
    ! -------------------------------------------------------------------
    DO idx = 1, nProcessSel
      iCase = iRowToCase(idx)

      ! Apply the Lambdas to the integrated values once
      Rout(idx) = lam1(idx) * Int_R1(iCase) + &
                  lam2(idx) * Int_R2(iCase) + &
                  lam3(idx) * Int_R3(iCase)
      
    END DO
    ! Clamp it to zero. This is the angle-INTEGRATED zeroth moment, which is a
    ! rate and must be non-negative. Note CalculateAllPhout deliberately does NOT
    ! clamp: it returns Legendre moments, and Phi_l for l >= 1 is legitimately
    ! negative over much of the (E1,E3) plane.
    Rout = MAX(Rout, 0.0d0)
    
  END SUBROUTINE CalculateAllRoutIntegrated

  SUBROUTINE CalculateAllR1R2R3( E1, E3, T, Mu_e, Mu_mu, iE1, iE3, iTh, R1, R2, R3 )

    REAL(DP), INTENT(IN)  :: E1, E3, T, Mu_e, Mu_mu
    INTEGER , INTENT(IN)  :: iE1, iE3, iTh
    ! R1/R2/R3 are indexed COMPACTLY (1..nDistinctCases), same as I0/I1/I2.
    REAL(DP), INTENT(OUT) :: R1(nDistinctCasesMax)
    REAL(DP), INTENT(OUT) :: R2(nDistinctCasesMax)
    REAL(DP), INTENT(OUT) :: R3(nDistinctCasesMax)

    REAL(DP) :: I0(nDistinctCasesMax), I1(nDistinctCasesMax), I2(nDistinctCasesMax)
    INTEGER  :: iCase, iMP

    ! Ensure unused cases are safely zeroed out
    R1 = 0.0d0; R2 = 0.0d0; R3 = 0.0d0

    CALL CalculateDistinctI0I1I2( E1, E3, T, Mu_e, Mu_mu, iE1, iE3, iTh, I0, I1, I2 )

    ! Guo Eq. 3. The 1/(16*pi*Delta^5) normalization is built into A/B/C at Init,
    ! and the IMD pair carries its own DeltaA^5 version in the same slots -- so
    ! there is no branch here. Indexing scalars
    ! rather than passing A(:)/B(:)/C(:) array sections also avoids a dummy
    ! temporary per call, which matters once this runs on a device.
    DO iCase = 1, nDistinctCases
      iMP = iMassPair_Case(iCase)
      R1(iCase) = A_MP(iE1,iE3,iTh,1,iMP)*I2(iCase) &
                + B_MP(iE1,iE3,iTh,1,iMP)*I1(iCase) &
                + C_MP(iE1,iE3,iTh,1,iMP)*I0(iCase)
      R2(iCase) = A_MP(iE1,iE3,iTh,2,iMP)*I2(iCase) &
                + B_MP(iE1,iE3,iTh,2,iMP)*I1(iCase) &
                + C_MP(iE1,iE3,iTh,2,iMP)*I0(iCase)
      R3(iCase) = C_MP(iE1,iE3,iTh,3,iMP)*I0(iCase)
    ENDDO

  END SUBROUTINE CalculateAllR1R2R3

  SUBROUTINE SetChemPotDistinctCases(Mu_e, Mu_mu, ChemPotPair)

    REAL(DP), INTENT(IN)  :: Mu_e, Mu_mu
    REAL(DP), INTENT(OUT) :: ChemPotPair(nDistinctCasesMax,2)

    ChemPotPair(i_em_em,1) = Mu_e
    ChemPotPair(i_em_em,2) = Mu_e

    ChemPotPair(i_ep_ep,1) = -Mu_e
    ChemPotPair(i_ep_ep,2) = -Mu_e

    ChemPotPair(i_mm_mm,1) = Mu_mu
    ChemPotPair(i_mm_mm,2) = Mu_mu

    ChemPotPair(i_mp_mp,1) = -Mu_mu
    ChemPotPair(i_mp_mp,2) = -Mu_mu

    ChemPotPair(i_em_mm,1) = Mu_e
    ChemPotPair(i_em_mm,2) = Mu_mu

    ChemPotPair(i_mm_em,1) = Mu_mu
    ChemPotPair(i_mm_em,2) = Mu_e

    ChemPotPair(i_ep_mp,1) = -Mu_e
    ChemPotPair(i_ep_mp,2) = -Mu_mu

    ChemPotPair(i_mp_ep,1) = -Mu_mu
    ChemPotPair(i_mp_ep,2) = -Mu_e

    ChemPotPair(i_IMD,1) = Mu_e
    ChemPotPair(i_IMD,2) = Mu_mu

    ChemPotPair(i_IMDbar,1) = -Mu_e
    ChemPotPair(i_IMDbar,2) = -Mu_mu

  END SUBROUTINE SetChemPotDistinctCases

  SUBROUTINE CalculateDistinctI0I1I2( E1, E3, T, Mu_e, Mu_mu, iE1, iE3, iTh, I0, I1, I2 )

    ! All kinematics (E_minus, E_plus, phase-space masks) were tabulated at Init:
    ! they depend only on (E1,E3,costh) and the case mass pair, never on T or the
    ! chemical potentials. What is left here is the only genuinely thermodynamic
    ! part, the incomplete Fermi-Dirac integrals.

    REAL(DP), INTENT(IN)  :: E1, E3, T, Mu_e, Mu_mu
    INTEGER , INTENT(IN)  :: iE1, iE3, iTh
    ! Indexed compactly (1..nDistinctCases). Declared at the PARAMETER bound so the
    ! dummy extent does not depend on mutable module state.
    REAL(DP), INTENT(OUT) :: I0(nDistinctCasesMax), I1(nDistinctCasesMax), I2(nDistinctCasesMax)

    REAL(DP) :: xMu_l2, xMu_l4, E3eff, MuSel(2)
    REAL(DP) :: I0p, I1p, I2p
    INTEGER  :: iCase, iMP

    I0 = 0.0d0
    I1 = 0.0d0
    I2 = 0.0d0

    MuSel(1) = Mu_e
    MuSel(2) = Mu_mu

    DO iCase = 1, nDistinctCases

      iMP = iMassPair_Case(iCase)
      IF ( .NOT. Active_MP(iE1,iE3,iTh,iMP) ) CYCLE

      xMu_l2 = ChemSign2(iCase) * MuSel( iChemSpecies2(iCase) )
      xMu_l4 = ChemSign4(iCase) * MuSel( iChemSpecies4(iCase) )

      ! IMD evaluates its I_s at (E1,-E3); scattering at (E1,E3).
      E3eff = E3Sign_Case(iCase) * E3

      CALL Compute_Is_Integrals( E1, E3eff, Eminus_MP(iE1,iE3,iTh,iMP), &
                                 T, xMu_l2, xMu_l4, I0(iCase), I1(iCase), I2(iCase) )

      IF ( TwoSided_Case(iCase) ) THEN
        ! I_s = I_s(E_minus) - I_s(E_plus)   (Guo Eq. 18)
        CALL Compute_Is_Integrals( E1, E3eff, Eplus_MP(iE1,iE3,iTh,iMP), &
                                   T, xMu_l2, xMu_l4, I0p, I1p, I2p )
        I0(iCase) = I0(iCase) - I0p
        I1(iCase) = I1(iCase) - I1p
        I2(iCase) = I2(iCase) - I2p
      ENDIF

    ENDDO

  END SUBROUTINE CalculateDistinctI0I1I2

  SUBROUTINE CalculateAllPhi0_collinear( iE1, E1, xT, Mu_e, Mu_mu, Phi0 )

    INTEGER,  INTENT(IN)  :: iE1
    REAL(DP), INTENT(IN)  :: E1, xT, Mu_e, Mu_mu
    REAL(DP), INTENT(OUT) :: Phi0(nProcessSel)

    INTEGER  :: iE2, iCase, iTrueCase, idx
    REAL(DP) :: E2, wE2, xMu_l2, xMu_l4, f2, f4, occ, xnode
    REAL(DP) :: val_H0I, val_H0II
    REAL(DP) :: ChemPotPairDistinct(nDistinctCasesMax,2)
    REAL(DP) :: integrand0I(nDistinctCasesMax), integrand0II(nDistinctCasesMax)

    integrand0I  = 0.0d0
    integrand0II = 0.0d0
    CALL SetChemPotDistinctCases(Mu_e, Mu_mu, ChemPotPairDistinct)

    DO iCase = 1, nDistinctCases
      iTrueCase = ChosenToTrueCaseMap(iCase)

      xMu_l2 = ChemPotPairDistinct(iTrueCase,1)
      xMu_l4 = ChemPotPairDistinct(iTrueCase,2)

      ! ================================================================
      ! PANEL 1: [0, E1] -- finite, use Gauss-Legendre.
      ! Affine-map the preallocated canonical [-1,1] grid to [0,E1].
      ! E2 = (E1/2)*(x_ref + 1),  weight scaled by (E1/2).
      ! ================================================================
      DO iE2 = 1, nE2_Lepton
        E2  = 0.5d0 * E1 * ( xa_GLeg_ref(iE2) + 1.0d0 )
        wE2 = 0.5d0 * E1 * wa_GLeg_ref(iE2)

        CALL CalculateCollinearH0(E1, E2, val_H0I, val_H0II)

        f2 = FD(E2, xMu_l2, xT)
        IF (ABS(xMu_l2 - xMu_l4) < 1.0d-10) THEN
          f4 = f2
        ELSE
          f4 = FD(E2, xMu_l4, xT)
        ENDIF
        occ = f2 * (1.0d0 - f4)

        integrand0I(iTrueCase)  = integrand0I(iTrueCase)  + wE2 * occ * val_H0I
        integrand0II(iTrueCase) = integrand0II(iTrueCase) + wE2 * occ * val_H0II
      ENDDO

      ! ================================================================
      ! PANEL 2: [E1, infinity) -- semi-infinite, use Gauss-Laguerre.
      ! Map raw GLag nodes (weight e^{-x}) via E2 = E1 + x*T.
      ! The built-in e^{-x} must be undone (multiply by e^{+x}) since the
      ! true integrand is f2*(1-f4)*H0, not a pure exponential; xT is the
      ! dE2 = xT dx Jacobian.
      ! ================================================================
      DO iE2 = 1, nE2_Lepton
        xnode = xa_GLag(iE2)
        E2    = E1 + xnode * xT
        wE2   = wa_GLag(iE2) * EXP(xnode) * xT

        CALL CalculateCollinearH0(E1, E2, val_H0I, val_H0II)

        f2 = FD(E2, xMu_l2, xT)
        IF (ABS(xMu_l2 - xMu_l4) < 1.0d-10) THEN
          f4 = f2
        ELSE
          f4 = FD(E2, xMu_l4, xT)
        ENDIF
        occ = f2 * (1.0d0 - f4)

        integrand0I(iTrueCase)  = integrand0I(iTrueCase)  + wE2 * occ * val_H0I
        integrand0II(iTrueCase) = integrand0II(iTrueCase) + wE2 * occ * val_H0II
      ENDDO

    ENDDO

    DO idx = 1, nProcessSel
      iTrueCase = iDistinctMap( ProcessList(idx) )

      Phi0(idx) = lam1(idx) * integrand0I(iTrueCase) + &
                  lam2(idx) * integrand0II(iTrueCase)
    END DO

    Phi0(:) = Phi0(:) / E1**4 * 4066026.8989943061_DP

  END SUBROUTINE CalculateAllPhi0_collinear

  SUBROUTINE CalculateAllPhoutGeneral_collinear( iE1, E1, xT, Mu_e, Mu_mu, Phout, nL )
    
    INTEGER , INTENT(IN)  :: iE1
    REAL(DP), INTENT(IN)  :: E1, xT, Mu_e, Mu_mu
    REAL(DP), INTENT(OUT) :: Phout(nProcessSel, nL)
    INTEGER , INTENT(IN)  :: nL

    REAL(DP) :: Phi0(nProcessSel)
    INTEGER  :: iL
    
    CALL CalculateAllPhi0_collinear( iE1, E1, xT, Mu_e, Mu_mu, Phi0 )

    Phout(:,1) = Phi0(:) * 0.5d0
    DO iL = 2, nL
      ! Here you have to reconstruct the kernel using eq. 42 from Mezzacappa & Bruenn 1993, 
      ! and then integrate over angle to get all other Phis. My guess is that from Phi0 you can
      ! directly calculate all Phis perhaps without having to do the Theta loop but maybe not.
      ! Don't have time now to implement this. Phi0 should be enough.
      Phout(:,iL) = 0.0d0
    ENDDO

  END SUBROUTINE CalculateAllPhoutGeneral_collinear

END MODULE wlGeneralLeptonScatteringModuleThornadoInterface
