PROGRAM wlTestAbEm2DFastSettings

  !=========================================================================
  !  Cheapest quadrature setting of the semi-leptonic charged-current 2D
  !  kernels (OpacitySource/GeneralEmAb) that meets a physically weighted
  !  accuracy metric. The methods are called through their explicit entry
  !  points (wlSemiLeptonicOpacityModule2D), independent of the build's
  !  default method.
  !
  !  Settings (all with the Opacity_CC_2D argument list and units):
  !    B "Momentum"  Opacity_CC_2D_Momentum_Nodes (default method):
  !                  nE2 outer nodes per outer segment, nE3 (<= 12) inner
  !                  nodes per kink-free inner segment; grid nE2_B x nE3_B
  !    A "Energy"    Opacity_CC_2D_Energy_Nodes ("standard", GSI-like
  !                  method, wlSemiLeptonicIntegration2D_Energy): linear
  !                  Gauss-Legendre, nE2 E2 nodes x nE3 E3 nodes; all pairs
  !                  of nE_A
  !    C "Tail"      Opacity_CC_2D_Tail_Nodes (wlSemiLeptonicIntegration2D_Tail):
  !                  no thermal cutoff, semi-infinite tail rule; nE_C
  !  The [0,1] Gauss-Legendre rules are built once per setting with gauleg,
  !  outside the timing, and passed to the stateless _Nodes kernels.
  !  Reference: C at nE_Ref = 16. If DoCheck, C at nE_Check = 24 is also
  !  run and scored against the reference like any setting; its max
  !  normalised error (<< 1 expected) measures the reference uncertainty.
  !
  !  Samples (as in wlTestAbEm2DConvergence):
  !    - every PointStride-th row of ThermoConditions.dat
  !    - nEnu log-spaced neutrino energies in [EnuMin, EnuMax]
  !    - electrons (Mue, me) and muons (Mum, mmu)
  !    - reactions 1-4 and corrections 0-3
  !
  !  Metric. Per (point, lepton l, energy, correction), from the reference:
  !    kappa_tot(nu_l)    = kappa_1
  !    kappa_tot(nubar_l) = kappa_2 + kappa_3
  !  Reaction 1 belongs to nu_l, reactions 2, 3, 4 to nubar_l. A sample is
  !  ignored if kappa_tot,ref of its species < KappaFloor (cm^-1).
  !  Weighted error:
  !    reactions 1-3: |kappa - kappa_ref| / kappa_tot,ref
  !    reaction 4   : |j - j_ref|/|j_ref| * kappa_3,ref / kappa_tot,ref
  !                   (emissivity of the channel of reaction 3, detailed
  !                   balance); j_ref = 0: error 0 if j = 0, else 1 * weight
  !                   (counted only if the weight >= ZeroSigWeight)
  !  NaN, or a zero where the reference is nonzero and its weight
  !  (kappa_ref/kappa_tot,ref, resp. kappa_3,ref/kappa_tot,ref) is
  !  >= ZeroSigWeight, is a failure (normalised error set to BadErr).
  !  Tolerance: TolLow for Enu <= E_split, TolHigh above, with
  !    E_split = max(mu_species, 0) + nT_Split * T,
  !    mu_nu_e = Mue + Mup - Mun, mu_nu_mu = Mum + Mup - Mun,
  !    mu_nubar = -mu_nu.
  !  Normalised error = weighted error / tolerance. A setting PASSES if the
  !  max normalised error over all non-ignored samples (all reactions,
  !  both leptons, corrections 0-3) is <= 1.
  !  For context the strict metric of wlTestAbEm2DConvergence is also
  !  given: max plain relative error |v - ref|/|ref| on corrections 0-2,
  !  on its significance mask (ref > 0, >= SigFrac * max over energies at
  !  the same (point, lepton, reaction, correction), >= ChanFrac * max of
  !  the channel; zero or NaN count as 1).
  !
  !  Reading the output:
  !    1. Per setting: us/eval (= wall * threads / samples), max
  !       normalised error at low and high E, % of scored samples failing,
  !       NaN/zero failures, max normalised error on corrections 0-2,
  !       strict max (corrections 0-2), worst sample (lepton, reaction,
  !       correction, Enu, T, Enu/E_split), PASS/FAIL on corrections 0-3
  !       (the metric) and on corrections 0-2 only (for context:
  !       correction 3, the q^2-dependent form factors clipped where
  !       |M|^2 < 0, is discontinuous and converges slowly in every scheme,
  !       also in the reference).
  !    2. All runs sorted by cost and the cheapest passing setting (and
  !       the cheapest passing on corrections 0-2 only).
  !    3. Sensitivity: cheapest passing setting for TolLow in TolLow_S x
  !       nT_Split in nTS_S (TolHigh fixed), from the stored values; for
  !       corrections 0-3 and 0-2.
  !    4. Worst samples of the cheapest passing setting (if none passes on
  !       corrections 0-3: of the cheapest passing on corrections 0-2, on
  !       corrections 0-2).
  !
  !  All values are stored (8 bytes x samples x runs, ~180 MB for all
  !  points). Output: stdout and, if WriteSamples, one line per scored
  !  sample with the normalised error of every run in
  !  AbEm2DFastSettings.dat (large).
  !=========================================================================

  USE wlKindModule, ONLY: dp
  USE wlSemiLeptonicOpacityModule2D, ONLY: &
    Opacity_CC_2D_Energy_Nodes, Opacity_CC_2D_Momentum_Nodes, Opacity_CC_2D_Tail_Nodes, gauleg
  USE wlEosConstantsModule, ONLY: &
    me, mmu
#if defined(WEAKLIB_OMP)
  USE omp_lib, ONLY: &
    omp_get_max_threads
#endif

  IMPLICIT NONE

  INTEGER , PARAMETER :: PointStride = 1        ! use every PointStride-th row
  INTEGER , PARAMETER :: nEnu  = 12             ! neutrino energies
  REAL(DP), PARAMETER :: EnuMin = 1.0d0, EnuMax = 300.0d0
  INTEGER , PARAMETER :: nLep = 2, nReac = 4, nCorr = 4

  ! ---- settings grid
  INTEGER , PARAMETER :: nB2 = 7, nB3 = 4
  INTEGER , PARAMETER :: nE2_B(nB2) = (/ 6, 8, 10, 12, 16, 20, 24 /)
  INTEGER , PARAMETER :: nE3_B(nB3) = (/ 4, 6, 8, 12 /)
  INTEGER , PARAMETER :: nA = 5                 ! A: all pairs (nE2, nE3) of nE_A
  INTEGER , PARAMETER :: nE_A(nA) = (/ 8, 12, 16, 20, 24 /)
  INTEGER , PARAMETER :: nC = 3
  INTEGER , PARAMETER :: nE_C(nC) = (/ 4, 6, 8 /)
  INTEGER , PARAMETER :: nE_Ref = 16            ! reference: C at this nE
  LOGICAL , PARAMETER :: DoCheck = .TRUE.       ! also run C at nE_Check
  INTEGER , PARAMETER :: nE_Check = 24

  ! ---- metric
  REAL(DP), PARAMETER :: KappaFloor = 1.0d-10   ! cm^-1, ignore species below
  REAL(DP), PARAMETER :: TolLow  = 1.0d-2       ! Enu <= E_split
  REAL(DP), PARAMETER :: TolHigh = 0.2d0        ! Enu >  E_split
  REAL(DP), PARAMETER :: nT_Split = 8.0d0       ! E_split = max(mu,0) + nT_Split*T
  REAL(DP), PARAMETER :: ZeroSigWeight = 1.0d-6 ! zero fails if its weight >= this
  REAL(DP), PARAMETER :: BadErr = 1.0d30        ! weighted error of a failure
  INTEGER , PARAMETER :: nTolS = 3, nTSS = 3    ! sensitivity grid
  REAL(DP), PARAMETER :: TolLow_S(nTolS) = (/ 3.0d-3, 1.0d-2, 3.0d-2 /)
  REAL(DP), PARAMETER :: nTS_S(nTSS) = (/ 5.0d0, 8.0d0, 12.0d0 /)

  ! ---- strict (old) metric, as in wlTestAbEm2DConvergence
  REAL(DP), PARAMETER :: SigFrac  = 1.0d-6
  REAL(DP), PARAMETER :: ChanFrac = 1.0d-15

  INTEGER , PARAMETER :: nWorst = 12
  LOGICAL , PARAMETER :: WriteSamples = .FALSE.
  CHARACTER(LEN=3), PARAMETER :: LepName(nLep) = (/ 'e  ', 'mu ' /)
  CHARACTER(LEN=1), PARAMETER :: SchTag(3) = (/ 'A', 'B', 'C' /)
  INTEGER , PARAMETER :: nRunMax = nB2*nB3 + nA*nA + nC + 2

  INTEGER  :: nThermoPoints, nPoints, nS, nThreads, nRun, nUsed
  INTEGER  :: i, iRow, k, k2, l, n, q, iRef, iChk, iBest, iBest02
  INTEGER  :: RunSch(nRunMax), RunN2(nRunMax), RunN3(nRunMax), Order(nRunMax)
  LOGICAL  :: RunIsSet(nRunMax)                 ! .FALSE. for reference / check
  REAL(DP) :: EnuA(nEnu), Row(12)
  REAL(DP) :: Wall(nRunMax), UsEval(nRunMax)
  REAL(DP) :: MaxLow(nRunMax), MaxHigh(nRunMax), FracFail(nRunMax), StrictMax(nRunMax)
  REAL(DP) :: MaxC02(nRunMax), xL, xH, fF
  INTEGER  :: nBad(nRunMax), iWorst(nRunMax), nB, jW
  LOGICAL  :: Pass(nRunMax), PassC02(nRunMax)

  REAL(DP), ALLOCATABLE :: Thermo(:,:)              ! (12, nPoints)
  INTEGER , ALLOCATABLE :: PointRow(:)
  REAL(DP), ALLOCATABLE :: Val(:,:)                 ! (nS, nRun)
  REAL(DP), ALLOCATABLE :: Ktot(:), WRef(:)         ! species total, weight
  LOGICAL , ALLOCATABLE :: Used(:), Sig(:)

  ! Log-spaced neutrino energies
  DO l = 1, nEnu
    EnuA(l) = EnuMin * (EnuMax/EnuMin)**( DBLE(l-1) / DBLE(nEnu-1) )
  END DO

  ! Run list: B grid, A grid, C list, reference, check
  n = 0
  DO k = 1, nB2
    DO k2 = 1, nB3
      n = n + 1; RunSch(n) = 2; RunN2(n) = nE2_B(k); RunN3(n) = nE3_B(k2)
    END DO
  END DO
  DO k = 1, nA
    DO k2 = 1, nA
      n = n + 1; RunSch(n) = 1; RunN2(n) = nE_A(k); RunN3(n) = nE_A(k2)
    END DO
  END DO
  DO k = 1, nC
    n = n + 1; RunSch(n) = 3; RunN2(n) = nE_C(k); RunN3(n) = 0
  END DO
  RunIsSet(1:n) = .TRUE.
  n = n + 1; RunSch(n) = 3; RunN2(n) = nE_Ref; RunN3(n) = 0
  RunIsSet(n) = .FALSE.
  iRef = n
  iChk = 0
  IF( DoCheck ) THEN
    n = n + 1; RunSch(n) = 3; RunN2(n) = nE_Check; RunN3(n) = 0
    RunIsSet(n) = .FALSE.
    iChk = n
  END IF
  nRun = n

  ! Read thermodynamic data
  ! Columns: T Rho Ye Ym Mue Mum Mun Mup Un Up M*n M*p
  ! (relativistic chemical potentials, MeV)
  OPEN(UNIT=123, FILE='ThermoConditions.dat', STATUS='OLD', ACTION='READ')
  READ(123,*) nThermoPoints
  READ(123,*)
  nPoints = ( nThermoPoints - 1 ) / PointStride + 1
  ALLOCATE( Thermo(12, nPoints), PointRow(nPoints) )
  i = 0
  DO iRow = 1, nThermoPoints
    READ(123,*) Row
    IF( MOD(iRow-1, PointStride) /= 0 ) CYCLE
    i = i + 1
    Thermo(:,i) = Row
    PointRow(i) = iRow
  END DO
  CLOSE(123)

  ! Sample layout (fastest first): energy, correction, reaction, lepton, point
  nS = nEnu*nCorr*nReac*nLep*nPoints
  ALLOCATE( Val(nS, nRun), Ktot(nS), WRef(nS), Used(nS), Sig(nS) )

  nThreads = 1
#if defined(WEAKLIB_OMP)
  nThreads = omp_get_max_threads()
#endif

  WRITE(*,'(A)') ' wlTestAbEm2DFastSettings: cheapest quadrature meeting a weighted accuracy metric'
  WRITE(*,'(A)') '   A = Energy (nE2 x nE3 linear GL), B = Momentum (nE2 outer x nE3 inner per segment),' // &
                 ' C = Tail (nE)'
  WRITE(*,'(A,I6,A,I3,A,I3,A,F6.1,A,F6.1,A)') ' Thermo points:', nPoints, '  (stride', PointStride, &
      ')   energies:', nEnu, ' in [', EnuMin, ',', EnuMax, '] MeV'
  WRITE(*,'(A,I9,A,I3,A,I3,A,F7.1,A)') ' Samples:', nS, '   runs:', nRun, '   threads:', nThreads, &
      '   value storage', 8.0d0*DBLE(nS)*DBLE(nRun)/1.0d6, ' MB'
  WRITE(*,'(A,A,A,ES8.1,A)') ' Reference: ', TRIM(RunLabel(iRef)), ';  KappaFloor', KappaFloor, ' cm^-1'
  WRITE(*,'(A,ES8.1,A,ES8.1,A,F5.1,A)') ' Tolerance: ', TolLow, ' for Enu <= E_split, ', TolHigh, &
      ' above;  E_split = max(mu,0) +', nT_Split, ' T'

  DO n = 1, nRun
    CALL ComputeRun( n )
  END DO
  UsEval(1:nRun) = Wall(1:nRun) * nThreads / DBLE(nS) * 1.0d6

  CALL SetReference
  CALL SetSignificance
  nUsed = COUNT( Used )
  WRITE(*,'(A,I9,A,I9,A)') ' Scored samples:', nUsed, '  (', nS - nUsed, ' ignored, kappa_tot,ref < KappaFloor)'

  DO n = 1, nRun
    CALL RunStats( n, TolLow, nT_Split, nCorr-1, MaxLow(n), MaxHigh(n), FracFail(n), nBad(n), iWorst(n) )
    Pass(n)      = MAX( MaxLow(n), MaxHigh(n) ) <= 1.0d0
    CALL RunStats( n, TolLow, nT_Split, 2, xL, xH, fF, nB, jW )
    MaxC02(n)    = MAX( xL, xH )
    PassC02(n)   = MaxC02(n) <= 1.0d0
    StrictMax(n) = StrictMaxOf( n )
  END DO

  ! ---- 1. per setting
  WRITE(*,'(/,A)') ' Per setting (normalised error = weighted error / tolerance; PASS if max <= 1):'
  CALL TableHeader
  DO n = 1, nRun
    CALL TableLine( n )
  END DO
  IF( iChk > 0 ) THEN
    WRITE(*,'(/,A,A,A,ES10.2,A,ES10.2)') ' Reference check: ', TRIM(RunLabel(iChk)), &
        ' vs reference, max normalised error', MAX( MaxLow(iChk), MaxHigh(iChk) ), &
        '  (should be << 1);  corrections 0-2:', MaxC02(iChk)
  END IF

  ! ---- 2. sorted by cost
  CALL SortByCost
  WRITE(*,'(/,A)') ' All runs sorted by cost:'
  CALL TableHeader
  DO q = 1, nRun
    CALL TableLine( Order(q) )
  END DO
  iBest   = Cheapest( Pass )
  iBest02 = Cheapest( PassC02 )
  IF( iBest > 0 ) THEN
    WRITE(*,'(/,A,A,A,F8.2,A,2ES10.2)') ' Cheapest passing setting: ', TRIM(RunLabel(iBest)), '   us/eval', &
        UsEval(iBest), '   max norm. err low/high E', MaxLow(iBest), MaxHigh(iBest)
  ELSE
    WRITE(*,'(/,A)') ' Cheapest passing setting: none'
  END IF
  IF( iBest02 > 0 ) THEN
    WRITE(*,'(A,A,A,F8.2,A,ES10.2,A)') ' (corrections 0-2 only:   ', TRIM(RunLabel(iBest02)), '   us/eval', &
        UsEval(iBest02), '   max norm. err', MaxC02(iBest02), ')'
  ELSE
    WRITE(*,'(A)') ' (corrections 0-2 only: none)'
  END IF

  ! ---- 3. sensitivity
  CALL Sensitivity

  ! ---- 4. worst samples of the cheapest passing setting
  IF( iBest > 0 ) THEN
    CALL ListWorst( iBest, nCorr-1 )
  ELSE IF( iBest02 > 0 ) THEN
    CALL ListWorst( iBest02, 2 )
  END IF

  IF( WriteSamples ) CALL WriteSampleFile

  DEALLOCATE( Thermo, PointRow, Val, Ktot, WRef, Used, Sig )

CONTAINS

  ! ---- sample index <-> (energy, correction, reaction, lepton, point)
  INTEGER FUNCTION EnuOf( j )
    INTEGER, INTENT(IN) :: j
    EnuOf = MOD( j-1, nEnu ) + 1
  END FUNCTION EnuOf

  INTEGER FUNCTION CorrOf( j )
    INTEGER, INTENT(IN) :: j
    CorrOf = MOD( (j-1)/nEnu, nCorr )
  END FUNCTION CorrOf

  INTEGER FUNCTION ReacOf( j )
    INTEGER, INTENT(IN) :: j
    ReacOf = MOD( (j-1)/(nEnu*nCorr), nReac ) + 1
  END FUNCTION ReacOf

  INTEGER FUNCTION LepOf( j )
    INTEGER, INTENT(IN) :: j
    LepOf = MOD( (j-1)/(nEnu*nCorr*nReac), nLep ) + 1
  END FUNCTION LepOf

  INTEGER FUNCTION PointOf( j )
    INTEGER, INTENT(IN) :: j
    PointOf = (j-1)/(nEnu*nCorr*nReac*nLep) + 1
  END FUNCTION PointOf

  INTEGER FUNCTION WithReac( j, ir )
    ! sample with the same (energy, correction, lepton, point), reaction ir
    INTEGER, INTENT(IN) :: j, ir
    WithReac = j + nEnu*nCorr*( ir - ReacOf(j) )
  END FUNCTION WithReac


  CHARACTER(LEN=8) FUNCTION RunLabel( nn )
    INTEGER, INTENT(IN) :: nn
    IF( RunSch(nn) == 3 ) THEN
      WRITE(RunLabel,'(A1,I3)') SchTag(3), RunN2(nn)
    ELSE
      WRITE(RunLabel,'(A1,I3,A1,I2)') SchTag(RunSch(nn)), RunN2(nn), 'x', RunN3(nn)
    END IF
  END FUNCTION RunLabel


  SUBROUTINE ComputeRun( nn )

    INTEGER, INTENT(IN) :: nn

    REAL(DP) :: x2(RunN2(nn)), w2(RunN2(nn)), x3(MAX(RunN3(nn),1)), w3(MAX(RunN3(nn),1))
    REAL(DP) :: Mul, ml
    INTEGER  :: j, ip, il, iSc, n2, n3, count_start, count_end, count_rate

    iSc = RunSch(nn)
    n2  = RunN2(nn)
    n3  = MAX( RunN3(nn), 1 )
    CALL gauleg( 0.0d0, 1.0d0, x2, w2, n2 )
    IF( iSc /= 3 ) CALL gauleg( 0.0d0, 1.0d0, x3, w3, n3 )

    CALL SYSTEM_CLOCK( count_start, count_rate )
#if defined(WEAKLIB_OMP)
    !$OMP PARALLEL DO SCHEDULE(DYNAMIC,16) PRIVATE( j, ip, il, Mul, ml )
#endif
    DO j = 1, nS
      ip = PointOf(j)
      il = LepOf(j)
      IF( il == 1 ) THEN
        Mul = Thermo(5,ip)
        ml  = me
      ELSE
        Mul = Thermo(6,ip)
        ml  = mmu
      END IF
      SELECT CASE( iSc )
      CASE( 1 )
        CALL Opacity_CC_2D_Energy_Nodes( CorrOf(j), ReacOf(j), EnuA(EnuOf(j)), Val(j,nn), &
            Thermo(1,ip), Mul, Thermo(7,ip), Thermo(8,ip), ml, Thermo(11,ip), Thermo(12,ip), &
            Thermo(9,ip), Thermo(10,ip), n2, x2, w2, n3, x3, w3 )
      CASE( 2 )
        CALL Opacity_CC_2D_Momentum_Nodes( CorrOf(j), ReacOf(j), EnuA(EnuOf(j)), Val(j,nn), &
            Thermo(1,ip), Mul, Thermo(7,ip), Thermo(8,ip), ml, Thermo(11,ip), Thermo(12,ip), &
            Thermo(9,ip), Thermo(10,ip), n2, x2, w2, n3, x3, w3 )
      CASE DEFAULT
        CALL Opacity_CC_2D_Tail_Nodes( CorrOf(j), ReacOf(j), EnuA(EnuOf(j)), Val(j,nn), &
            Thermo(1,ip), Mul, Thermo(7,ip), Thermo(8,ip), ml, Thermo(11,ip), Thermo(12,ip), &
            Thermo(9,ip), Thermo(10,ip), n2, x2, w2 )
      END SELECT
    END DO
#if defined(WEAKLIB_OMP)
    !$OMP END PARALLEL DO
#endif
    CALL SYSTEM_CLOCK( count_end )
    Wall(nn) = DBLE( count_end - count_start ) / DBLE( count_rate )

    WRITE(*,'(A,A8,A,F9.2,A)') '   computed ', RunLabel(nn), '  in', Wall(nn), ' s'

  END SUBROUTINE ComputeRun


  SUBROUTINE SetReference

    ! species totals and weights of the reference; mask of scored samples
    INTEGER  :: j, ir
    REAL(DP) :: k1, k2, k3, kt

    DO j = 1, nS
      ir = ReacOf(j)
      k1 = Val(WithReac(j,1), iRef)
      k2 = Val(WithReac(j,2), iRef)
      k3 = Val(WithReac(j,3), iRef)
      IF( ir == 1 ) THEN
        kt = k1
      ELSE
        kt = k2 + k3
      END IF
      Ktot(j) = kt
      Used(j) = kt >= KappaFloor               ! .FALSE. also for NaN
      IF( .NOT. Used(j) ) THEN
        WRef(j) = 0.0d0
      ELSE IF( ir == 4 ) THEN
        WRef(j) = ABS( k3 ) / kt
      ELSE
        WRef(j) = ABS( Val(j,iRef) ) / kt
      END IF
    END DO

  END SUBROUTINE SetReference


  SUBROUTINE SetSignificance

    ! strict metric mask, as in wlTestAbEm2DConvergence
    REAL(DP) :: MaxE(nS/nEnu), MaxChan(nCorr, nReac, nLep), r
    INTEGER  :: j, b

    MaxE    = 0.0d0
    MaxChan = 0.0d0
    DO j = 1, nS
      b = (j-1)/nEnu + 1
      r = Val(j,iRef)
      IF( r > 0.0d0 ) THEN
        MaxE(b) = MAX( MaxE(b), r )
        MaxChan(CorrOf(j)+1, ReacOf(j), LepOf(j)) = &
          MAX( MaxChan(CorrOf(j)+1, ReacOf(j), LepOf(j)), r )
      END IF
    END DO
    DO j = 1, nS
      b = (j-1)/nEnu + 1
      r = Val(j,iRef)
      Sig(j) = r > 0.0d0 .AND. r >= SigFrac*MaxE(b) .AND. &
               r >= ChanFrac*MaxChan(CorrOf(j)+1, ReacOf(j), LepOf(j))
    END DO

  END SUBROUTINE SetSignificance


  REAL(DP) FUNCTION WErr( j, nn )
    ! weighted error of the scored sample j in run nn; BadErr for a failure
    INTEGER, INTENT(IN) :: j, nn
    REAL(DP) :: v, r

    v = Val(j,nn)
    r = Val(j,iRef)
    IF( v /= v ) THEN
      WErr = BadErr
    ELSE IF( v == 0.0d0 .AND. r /= 0.0d0 .AND. WRef(j) >= ZeroSigWeight ) THEN
      WErr = BadErr
    ELSE IF( ReacOf(j) == 4 ) THEN
      IF( r == 0.0d0 ) THEN
        IF( v == 0.0d0 .OR. WRef(j) < ZeroSigWeight ) THEN
          WErr = 0.0d0
        ELSE
          WErr = WRef(j)
        END IF
      ELSE
        WErr = ABS( v - r ) / ABS( r ) * WRef(j)
      END IF
    ELSE
      WErr = ABS( v - r ) / Ktot(j)
    END IF

  END FUNCTION WErr


  REAL(DP) FUNCTION ESplit( j, nTS )
    ! max(mu_species, 0) + nTS*T, mu_nu_l = Mul + Mup - Mun, mu_nubar = -mu_nu
    INTEGER , INTENT(IN) :: j
    REAL(DP), INTENT(IN) :: nTS
    INTEGER  :: ip
    REAL(DP) :: mu
    ip = PointOf(j)
    mu = Thermo(4+LepOf(j),ip) + Thermo(8,ip) - Thermo(7,ip)
    IF( ReacOf(j) /= 1 ) mu = -mu
    ESplit = MAX( mu, 0.0d0 ) + nTS*Thermo(1,ip)
  END FUNCTION ESplit


  REAL(DP) FUNCTION NormErr( j, nn, TL, nTS )
    INTEGER , INTENT(IN) :: j, nn
    REAL(DP), INTENT(IN) :: TL, nTS
    IF( EnuA(EnuOf(j)) <= ESplit( j, nTS ) ) THEN
      NormErr = WErr( j, nn ) / TL
    ELSE
      NormErr = WErr( j, nn ) / TolHigh
    END IF
  END FUNCTION NormErr


  SUBROUTINE RunStats( nn, TL, nTS, cMax, xLow, xHigh, fFail, nB, jW )

    ! max normalised error at low/high E, fraction failing, number of
    ! NaN/zero failures and worst sample of run nn for (TolLow, nT_Split)
    ! = (TL, nTS), over corrections 0..cMax
    INTEGER , INTENT(IN)  :: nn, cMax
    REAL(DP), INTENT(IN)  :: TL, nTS
    REAL(DP), INTENT(OUT) :: xLow, xHigh, fFail
    INTEGER , INTENT(OUT) :: nB, jW
    INTEGER  :: j, nF
    REAL(DP) :: e, we, xW

    xLow  = 0.0d0
    xHigh = 0.0d0
    xW    = -1.0d0
    nF = 0
    nB = 0
    jW = 0
    DO j = 1, nS
      IF( .NOT. Used(j) .OR. CorrOf(j) > cMax ) CYCLE
      we = WErr( j, nn )
      IF( we >= BadErr ) nB = nB + 1
      IF( EnuA(EnuOf(j)) <= ESplit( j, nTS ) ) THEN
        e = we / TL
        xLow = MAX( xLow, e )
      ELSE
        e = we / TolHigh
        xHigh = MAX( xHigh, e )
      END IF
      IF( e > 1.0d0 ) nF = nF + 1
      IF( e > xW ) THEN
        xW = e
        jW = j
      END IF
    END DO
    fFail = DBLE( nF ) / DBLE( MAX( nUsed, 1 ) )

  END SUBROUTINE RunStats


  REAL(DP) FUNCTION StrictMaxOf( nn )
    ! max plain relative error, corrections 0-2, significant samples
    INTEGER, INTENT(IN) :: nn
    INTEGER  :: j
    REAL(DP) :: v, e
    StrictMaxOf = 0.0d0
    DO j = 1, nS
      IF( .NOT. Sig(j) .OR. CorrOf(j) > 2 ) CYCLE
      v = Val(j,nn)
      IF( v /= v .OR. v == 0.0d0 ) THEN
        e = 1.0d0
      ELSE
        e = ABS( v - Val(j,iRef) ) / Val(j,iRef)
      END IF
      StrictMaxOf = MAX( StrictMaxOf, e )
    END DO
  END FUNCTION StrictMaxOf


  SUBROUTINE TableHeader
    WRITE(*,'(A)') '  run       us/eval | max norm.err: low E    high E  fail[%] NaN/0 |   C0-2   | strict C0-2 |' // &
                   ' worst: lep reac corr Enu[MeV]  T[MeV] E/Esplit | all  C0-2'
  END SUBROUTINE TableHeader


  SUBROUTINE TableLine( nn )
    INTEGER, INTENT(IN) :: nn
    INTEGER :: j
    CHARACTER(LEN=4) :: Verdict, V02
    j = iWorst(nn)
    IF( nn == iRef ) THEN
      Verdict = 'ref'
    ELSE IF( nn == iChk ) THEN
      Verdict = 'chk'
    ELSE IF( Pass(nn) ) THEN
      Verdict = 'PASS'
    ELSE
      Verdict = 'FAIL'
    END IF
    V02 = Verdict
    IF( RunIsSet(nn) ) THEN
      V02 = 'FAIL'
      IF( PassC02(nn) ) V02 = 'PASS'
    END IF
    IF( j > 0 .AND. nn /= iRef ) THEN
      WRITE(*,'(2X,A8,F10.2,A,12X,2ES10.2,F8.3,I6,A,ES9.2,A,ES11.2,A,A9,2I5,F9.2,F8.2,F9.2,A,A5,1X,A4)') &
          RunLabel(nn), UsEval(nn), ' |', MaxLow(nn), MaxHigh(nn), 100.0d0*FracFail(nn), nBad(nn), ' |', &
          MaxC02(nn), ' |', StrictMax(nn), ' |', TRIM(LepName(LepOf(j))), ReacOf(j), CorrOf(j), EnuA(EnuOf(j)), &
          Thermo(1,PointOf(j)), EnuA(EnuOf(j))/ESplit( j, nT_Split ), ' |', Verdict, V02
    ELSE
      WRITE(*,'(2X,A8,F10.2,A,12X,2ES10.2,F8.3,I6,A,ES9.2,A,ES11.2,A,A45,A,A5,1X,A4)') &
          RunLabel(nn), UsEval(nn), ' |', MaxLow(nn), MaxHigh(nn), 100.0d0*FracFail(nn), nBad(nn), ' |', &
          MaxC02(nn), ' |', StrictMax(nn), ' |', ' ', ' |', Verdict, V02
    END IF
  END SUBROUTINE TableLine


  SUBROUTINE SortByCost
    ! insertion sort of the run indices by us/eval
    INTEGER :: a, b, t
    DO a = 1, nRun
      Order(a) = a
    END DO
    DO a = 2, nRun
      t = Order(a)
      b = a - 1
      DO WHILE( b >= 1 )
        IF( UsEval(Order(b)) <= UsEval(t) ) EXIT
        Order(b+1) = Order(b)
        b = b - 1
      END DO
      Order(b+1) = t
    END DO
  END SUBROUTINE SortByCost


  INTEGER FUNCTION Cheapest( Ok )
    ! cheapest setting with Ok (reference and check excluded)
    LOGICAL, INTENT(IN) :: Ok(nRunMax)
    INTEGER :: nn
    Cheapest = 0
    DO nn = 1, nRun
      IF( .NOT. RunIsSet(nn) .OR. .NOT. Ok(nn) ) CYCLE
      IF( Cheapest == 0 ) THEN
        Cheapest = nn
      ELSE IF( UsEval(nn) < UsEval(Cheapest) ) THEN
        Cheapest = nn
      END IF
    END DO
  END FUNCTION Cheapest


  SUBROUTINE Sensitivity

    LOGICAL  :: Ok(nRunMax)
    REAL(DP) :: sL, sH, sF
    INTEGER  :: nn, sB, sW, iT, iTS, iB(nTSS), cM

    DO cM = nCorr-1, 2, -1
      IF( cM == nCorr-1 ) THEN
        WRITE(*,'(/,A,ES8.1,A)') ' Sensitivity: cheapest passing setting (us/eval); TolHigh =', TolHigh, &
            ' fixed; corrections 0-3'
      ELSE
        WRITE(*,'(/,A)') ' Same, corrections 0-2 only:'
      END IF
      WRITE(*,'(A,3(F14.0,6X))') '  TolLow \ nT_Split', ( nTS_S(iTS), iTS = 1, nTSS )
      DO iT = 1, nTolS
        DO iTS = 1, nTSS
          Ok = .FALSE.
          DO nn = 1, nRun
            IF( .NOT. RunIsSet(nn) ) CYCLE
            CALL RunStats( nn, TolLow_S(iT), nTS_S(iTS), cM, sL, sH, sF, sB, sW )
            Ok(nn) = MAX( sL, sH ) <= 1.0d0
          END DO
          iB(iTS) = Cheapest( Ok )
        END DO
        WRITE(*,'(ES12.1,8X,3A20)') TolLow_S(iT), ( SensCell( iB(iTS) ), iTS = 1, nTSS )
      END DO
    END DO

  END SUBROUTINE Sensitivity


  CHARACTER(LEN=20) FUNCTION SensCell( nn )
    INTEGER, INTENT(IN) :: nn
    IF( nn == 0 ) THEN
      SensCell = '      none'
    ELSE
      WRITE(SensCell,'(A10,A2,F6.2,A1)') RunLabel(nn), ' (', UsEval(nn), ')'
    END IF
  END FUNCTION SensCell


  SUBROUTINE ListWorst( nn, cMax )

    INTEGER, INTENT(IN) :: nn, cMax
    INTEGER  :: j, q, r, Idx(nWorst)
    REAL(DP) :: Err(nWorst), e

    Err = -1.0d0
    Idx = 0
    DO j = 1, nS
      IF( .NOT. Used(j) .OR. CorrOf(j) > cMax ) CYCLE
      e = NormErr( j, nn, TolLow, nT_Split )
      IF( e <= Err(nWorst) ) CYCLE
      q = nWorst
      DO WHILE( q > 1 )
        IF( Err(q-1) >= e ) EXIT
        q = q - 1
      END DO
      DO r = nWorst, q+1, -1
        Err(r) = Err(r-1)
        Idx(r) = Idx(r-1)
      END DO
      Err(q) = e
      Idx(q) = j
    END DO

    WRITE(*,'(/,A,I1,A,A,A)') ' Worst samples (corrections 0-', cMax, ') of the cheapest passing setting ', &
        TRIM(RunLabel(nn)), ':'
    WRITE(*,'(A)') '   row lep reac corr  Enu[MeV]    T[MeV]  E/Esplit       ref   kappa_tot   rel.err' // &
                   '  norm.err'
    DO q = 1, nWorst
      j = Idx(q)
      IF( j == 0 ) EXIT
      WRITE(*,'(I6,A4,I5,I5,3F10.3,2ES12.3,2ES10.2)') PointRow(PointOf(j)), &
          TRIM(LepName(LepOf(j))), ReacOf(j), CorrOf(j), EnuA(EnuOf(j)), Thermo(1,PointOf(j)), &
          EnuA(EnuOf(j))/ESplit( j, nT_Split ), Val(j,iRef), Ktot(j), PlainRel( j, nn ), Err(q)
    END DO
    WRITE(*,'(A)') '   (ref in cm^-1, reaction 4 in (s cm^3 MeV)^-1; kappa_tot of the species in cm^-1)'

  END SUBROUTINE ListWorst


  REAL(DP) FUNCTION PlainRel( j, nn )
    INTEGER, INTENT(IN) :: j, nn
    IF( Val(j,iRef) == 0.0d0 ) THEN
      PlainRel = 0.0d0
      IF( Val(j,nn) /= 0.0d0 ) PlainRel = 1.0d0
    ELSE
      PlainRel = ABS( Val(j,nn) - Val(j,iRef) ) / ABS( Val(j,iRef) )
    END IF
  END FUNCTION PlainRel


  SUBROUTINE WriteSampleFile

    INTEGER :: j, nn
    CHARACTER(LEN=40) :: Fmt

    OPEN(UNIT=124, FILE='AbEm2DFastSettings.dat', STATUS='REPLACE', ACTION='WRITE')
    WRITE(124,'(A)') '# scored samples only; normalised error vs reference ' // TRIM(RunLabel(iRef))
    WRITE(124,'(A)') '# units of ref: cm^-1 (reactions 1-3), (s cm^3 MeV)^-1 (reaction 4); kappa_tot cm^-1'
    WRITE(124,'(A,*(A11))') '# row lep reac corr Enu[MeV] T[MeV] E_split ref kappa_tot   norm.err for', &
        ( RunLabel(nn), nn = 1, nRun )
    WRITE(Fmt,'(A,I0,A)') '(I6,3I4,3ES13.5,2ES14.6,', nRun, 'ES11.3)'
    DO j = 1, nS
      IF( .NOT. Used(j) ) CYCLE
      WRITE(124,Fmt) PointRow(PointOf(j)), LepOf(j), ReacOf(j), CorrOf(j), &
          EnuA(EnuOf(j)), Thermo(1,PointOf(j)), ESplit( j, nT_Split ), Val(j,iRef), Ktot(j), &
          ( NormErr( j, nn, TolLow, nT_Split ), nn = 1, nRun )
    END DO
    CLOSE(124)
    WRITE(*,'(/,A)') ' Per-sample normalised errors written to AbEm2DFastSettings.dat'

  END SUBROUTINE WriteSampleFile

END PROGRAM wlTestAbEm2DFastSettings
