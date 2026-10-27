PROGRAM wlTestEmAb2DCompareSchemes

  !=========================================================================
  !  Error vs cost of the three 2D integration methods of the
  !  semi-leptonic charged-current kernels (OpacitySource/GeneralEmAb/
  !  IntegrationMethods), each selected with InitOpacity_CC_2D(nE, method)
  !  and evaluated with Opacity_CC_2D (independent of the build's default
  !  method); nE is the resolution knob of each method:
  !
  !    A "Energy"    method iEmAb2D_Energy, the "standard" (GSI-like)
  !                  method: linear Gauss-Legendre in E2 and E3, nE x nE
  !    B "Momentum"  method iEmAb2D_Momentum, the default: nE outer
  !                  nodes per segment in p2, kink-split inner rule in p3
  !                  with MAX(6,MIN(12,nE/4)) nodes per segment
  !    C "Tail"      method iEmAb2D_Tail: no thermal cutoff, Gauss-Laguerre
  !                  tail
  !
  !  Samples (as in wlTestEmAb2DConvergence):
  !    - every PointStride-th row of ThermoConditions.dat
  !    - nEnu log-spaced neutrino energies in [EnuMin, EnuMax]
  !    - electrons (Mue, me) and muons (Mum, mmu)
  !    - reactions 1-4 and corrections 0-3
  !  Every scheme is run for its resolution list (nE_A, nE_B, nE_C) and
  !  at a high resolution nE_Hi (A and B clamp nE to nQuadMax = 256; the
  !  cost of A and C grows like nE**2, so their nE_Hi are lower).
  !  Timings exclude the node generation (done once by InitOpacity_CC_2D).
  !
  !  Reference: the high-resolution run of scheme iRefScheme (default C at
  !  nE = 24). C has no thermal cutoff and agrees with an independent
  !  1024x1024 linear Gauss-Legendre reference (Tfac = 150) to ~1e-6 for
  !  corrections 0-2 from nE = 12 on. B and A keep the Tfac = 30 cutoff and
  !  fixed/linear inner rules, which leave floors of a few 1e-5 (B) and
  !  more (A). Correction 3 does not converge in any scheme (the
  !  form-factor term of |M|^2 is discontinuous), so its errors are not
  !  resolved below ~1e-2. The "reference agreement" table gives the pairwise
  !  differences between the three high-resolution runs, normalised by
  !  the reference, i.e. an estimate of the reference uncertainty: errors
  !  below those numbers are not meaningful.
  !
  !  Relative errors |v - ref|/|ref| are evaluated on "significant"
  !  samples only:
  !    ref > 0,
  !    ref >= SigFrac  * max over energies at the same (point, lepton,
  !                      reaction, correction), and
  !    ref >= ChanFrac * max of the (lepton, reaction, correction)
  !                      channel over all points and energies.
  !  v = 0 or NaN on a significant sample counts as relative error 1 and
  !  is counted in the zero/NaN columns (all corrections).
  !  Correction 3 (q^2-dependent form factors, clipped where |M|^2 < 0)
  !  converges more slowly and is reported separately.
  !
  !  Reading the output:
  !    1. Reference agreement: floor below which errors are not resolved.
  !    2. Per-scheme tables: corrections 0-2 max/p99/median, max per
  !       reaction and per lepton, correction 3 max/median, zeros/NaNs,
  !       wall time and us/eval = wall * threads / samples.
  !    3. All runs sorted by cost (us/eval), with speedup vs B at nE_Base.
  !    4. For each tolerance, the cheapest (scheme, nE) whose max (resp.
  !       p99) error for corrections 0-2 is below it (reference excluded).
  !
  !  Units: cm^-1 (reactions 1-3) and (s cm^3 MeV)^-1 (reaction 4); only
  !  relative errors are reported.
  !
  !  Output: summary on stdout and, if WriteSamples, one line per
  !  significant sample in EmAb2DCompareSchemes.dat (large: ~60 MB for
  !  all points).
  !=========================================================================

  USE wlKindModule, ONLY: dp
  USE wlSemiLeptonicOpacityModule2D, ONLY: &
    Opacity_CC_2D, InitOpacity_CC_2D, iEmAb2D_Energy, iEmAb2D_Momentum, iEmAb2D_Tail
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
  INTEGER , PARAMETER :: nSch = 3
  CHARACTER(LEN=9), PARAMETER :: SchName(nSch) = (/ 'Energy   ', 'Momentum ', 'Tail     ' /)
  CHARACTER(LEN=1), PARAMETER :: SchTag(nSch)  = (/ 'A', 'B', 'C' /)
  INTEGER , PARAMETER :: nA = 6, nB = 9, nC = 5
  INTEGER , PARAMETER :: nE_A(nA) = (/ 16, 24, 32, 48, 64, 96 /)
  INTEGER , PARAMETER :: nE_B(nB) = (/ 8, 12, 16, 20, 24, 32, 48, 50, 64 /)
  INTEGER , PARAMETER :: nE_C(nC) = (/ 6, 8, 10, 12, 16 /)
  INTEGER , PARAMETER :: nE_Hi(nSch) = (/ 128, 256, 24 /)    ! high-resolution run per scheme
  INTEGER , PARAMETER :: iRefScheme = 3         ! reference = high-res run of this scheme
  INTEGER , PARAMETER :: nE_Base = 50           ! speedups relative to scheme B at this nE
  REAL(DP), PARAMETER :: SigFrac  = 1.0d-6
  REAL(DP), PARAMETER :: ChanFrac = 1.0d-15
  INTEGER , PARAMETER :: nTol = 5
  REAL(DP), PARAMETER :: Tol(nTol) = (/ 1.0d-2, 3.0d-3, 1.0d-3, 3.0d-4, 1.0d-4 /)
  LOGICAL , PARAMETER :: WriteSamples = .FALSE.
  CHARACTER(LEN=3), PARAMETER :: LepName(nLep) = (/ 'e  ', 'mu ' /)
  INTEGER , PARAMETER :: nRun = nA + nB + nC + nSch

  INTEGER  :: nThermoPoints, nPoints, nS, nThreads
  INTEGER  :: i, iRow, iLep, k, l, iS, m, n, n2, iSch, iTol, iRef, iBase
  INTEGER  :: RunSch(nRun), RunNE(nRun), iHi(nSch), Order(nRun)
  LOGICAL  :: RunIsHi(nRun)
  REAL(DP) :: EnuA(nEnu), Row(12)
  REAL(DP) :: xMax, xP99, xMed, yMax, yP99, yMed

  REAL(DP), ALLOCATABLE :: Thermo(:,:)              ! (12, nPoints)
  INTEGER , ALLOCATABLE :: PointRow(:)
  REAL(DP), ALLOCATABLE :: Val(:,:)                 ! (nS, nRun)
  LOGICAL , ALLOCATABLE :: Sig(:)
  REAL(DP), ALLOCATABLE :: Work(:)
  REAL(DP) :: Wall(nRun), UsEval(nRun)
  REAL(DP) :: MaxC02(nRun), P99C02(nRun), MedC02(nRun), MaxC3(nRun)
  INTEGER  :: nZero(nRun), nNaN(nRun), nSig, nSigC02

  ! Log-spaced neutrino energies
  DO l = 1, nEnu
    EnuA(l) = EnuMin * (EnuMax/EnuMin)**( DBLE(l-1) / DBLE(nEnu-1) )
  END DO

  ! Run list: per scheme its resolution list, then its high-resolution run
  n = 0
  DO iSch = 1, nSch
    SELECT CASE( iSch )
    CASE( 1 )
      DO k = 1, nA
        n = n + 1; RunSch(n) = 1; RunNE(n) = nE_A(k)
      END DO
    CASE( 2 )
      DO k = 1, nB
        n = n + 1; RunSch(n) = 2; RunNE(n) = nE_B(k)
      END DO
    CASE( 3 )
      DO k = 1, nC
        n = n + 1; RunSch(n) = 3; RunNE(n) = nE_C(k)
      END DO
    END SELECT
    n = n + 1; RunSch(n) = iSch; RunNE(n) = nE_Hi(iSch)
    iHi(iSch) = n
  END DO
  RunIsHi = .FALSE.
  RunIsHi(iHi) = .TRUE.
  iRef  = iHi(iRefScheme)
  iBase = 0
  DO n = 1, nRun
    IF( RunSch(n) == 2 .AND. RunNE(n) == nE_Base .AND. .NOT. RunIsHi(n) ) iBase = n
  END DO

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
  ALLOCATE( Val(nS, nRun), Sig(nS), Work(nS) )

  nThreads = 1
#if defined(WEAKLIB_OMP)
  nThreads = omp_get_max_threads()
#endif

  WRITE(*,'(A)') ' wlTestEmAb2DCompareSchemes: error vs cost of the 2D quadrature schemes'
  WRITE(*,'(A)') '   A = Energy (nE x nE linear GL), B = Momentum (default), C = Tail (no Tfac)'
  WRITE(*,'(A,I6,A,I3,A,I3,A,F6.1,A,F6.1,A)') ' Thermo points:', nPoints, '  (stride', PointStride, &
      ')   energies:', nEnu, ' in [', EnuMin, ',', EnuMax, '] MeV'
  WRITE(*,'(A,I9,A,I3,A,I3)') ' Samples:', nS, '   runs:', nRun, '   threads:', nThreads
  WRITE(*,'(A,A,A,F7.1,A)') ' Reference: ', TRIM(RunLabel(iRef)), ';  value storage', &
      8.0d0*DBLE(nS)*DBLE(nRun)/1.0d6, ' MB'

  DO n = 1, nRun
    CALL ComputeRun( n )
  END DO
  UsEval = Wall * nThreads / DBLE(nS) * 1.0d6

  CALL SetSignificance

  nSig    = COUNT( Sig )
  nSigC02 = 0
  DO iS = 1, nS
    IF( Sig(iS) .AND. CorrOf(iS) <= 2 ) nSigC02 = nSigC02 + 1
  END DO
  WRITE(*,'(A,I9,A,I9,A)') ' Significant samples:', nSig, '  (', nSigC02, ' with corrections 0-2)'
  WRITE(*,'(A,ES8.1,A,ES8.1,A)') ' sig: ref >', SigFrac, ' * max_E(point) and >', ChanFrac, &
      ' * max(channel)'

  ! Zeros and NaNs on significant samples, summary stats per run
  DO n = 1, nRun
    nZero(n) = 0
    nNaN(n)  = 0
    DO iS = 1, nS
      IF( .NOT. Sig(iS) ) CYCLE
      IF( Val(iS,n) /= Val(iS,n) ) THEN
        nNaN(n) = nNaN(n) + 1
      ELSE IF( Val(iS,n) == 0.0d0 ) THEN
        nZero(n) = nZero(n) + 1
      END IF
    END DO
    CALL GroupStats( n, iRef, -1, -1, 0, 2, m, MaxC02(n), P99C02(n), MedC02(n) )
    MaxC3(n) = GroupMax( n, iRef, -1, -1, 3, 3 )
  END DO

  ! ---- reference agreement: pairwise differences of the high-resolution runs
  WRITE(*,'(/,A)') ' Reference agreement: |X - Y| / ref between the high-resolution runs'
  WRITE(*,'(A)')   '   (errors below these numbers are not resolved by the reference)'
  WRITE(*,'(A)')   '                 | corrections 0-2                | correction 3'
  WRITE(*,'(A)')   '   pair          |  max       p99       median   |  max       p99       median'
  DO iSch = 1, nSch
    DO k = iSch+1, nSch
      CALL GroupStats( iHi(iSch), iHi(k), -1, -1, 0, 2, m, xMax, xP99, xMed )
      CALL GroupStats( iHi(iSch), iHi(k), -1, -1, 3, 3, m, yMax, yP99, yMed )
      WRITE(*,'(3X,A5,A3,A5,A,3ES10.2,A,3ES10.2)') RunLabelShort(iHi(iSch)), ' - ', &
          RunLabelShort(iHi(k)), ' |', xMax, xP99, xMed, ' |', yMax, yP99, yMed
    END DO
  END DO

  ! ---- per-scheme tables
  DO iSch = 1, nSch
    WRITE(*,'(/,A,A,A,A,A)') ' Scheme ', SchTag(iSch), ' (', TRIM(SchName(iSch)), &
        '): relative error vs reference'
    WRITE(*,'(A)') '        | corrections 0-2                ' // &
        '| max per reaction (corr 0-2)             | max corr 0-2        | corr 3              ' // &
        '| sig. samples |'
    WRITE(*,'(A)') '     nE |  max       p99       median   ' // &
        '|  reac 1    reac 2    reac 3    reac 4   |  e         mu       |  max       median   ' // &
        '|  zero   NaN  | wall[s]   us/eval'
    DO n = 1, nRun
      IF( RunSch(n) /= iSch ) CYCLE
      CALL GroupStats( n, iRef, -1, -1, 3, 3, m, yMax, yP99, yMed )
      WRITE(*,'(I6,A1,A,3ES10.2,A,4ES10.2,A,2ES10.2,A,2ES10.2,A,2I7,A,F8.2,F10.2)') RunNE(n), &
          HiMark(n), '|', MaxC02(n), P99C02(n), MedC02(n), ' |', &
          ( GroupMax( n, iRef, -1, k, 0, 2 ), k = 1, nReac ), ' |', &
          ( GroupMax( n, iRef, iLep, -1, 0, 2 ), iLep = 1, nLep ), ' |', &
          yMax, yMed, ' |', nZero(n), nNaN(n), ' |', Wall(n), UsEval(n)
    END DO
  END DO
  WRITE(*,'(A)') ' (* = high-resolution run; R = reference; us/eval = wall x threads / samples)'

  ! ---- all runs sorted by cost
  DO n = 1, nRun
    Order(n) = n
  END DO
  DO n = 2, nRun
    k = Order(n)
    n2 = n - 1
    DO WHILE( n2 >= 1 )
      IF( UsEval(Order(n2)) <= UsEval(k) ) EXIT
      Order(n2+1) = Order(n2)
      n2 = n2 - 1
    END DO
    Order(n2+1) = k
  END DO

  WRITE(*,'(/,A,I3,A)') ' All runs sorted by cost (speedup relative to B at nE =', nE_Base, ')'
  WRITE(*,'(A)') '                                       | corrections 0-2               | corr 3'
  WRITE(*,'(A)') '   run                us/eval  speedup |  max       p99       median   |  max'
  DO l = 1, nRun
    n = Order(l)
    WRITE(*,'(3X,A15,A1,F10.2,F9.2,A,3ES10.2,A,ES10.2)') RunLabel(n), HiMark(n), UsEval(n), &
        SpeedupOf(n), ' |', MaxC02(n), P99C02(n), MedC02(n), ' |', MaxC3(n)
  END DO

  ! ---- cheapest run meeting each tolerance
  WRITE(*,'(/,A)') ' Cheapest run meeting the tolerance (corrections 0-2, all reactions, both leptons;'
  WRITE(*,'(A,I3,A)') ' reference run excluded; speedup relative to B at nE =', nE_Base, '):'
  WRITE(*,'(A)')   '           | by max                             | by p99'
  WRITE(*,'(A)')   '       tol | run              us/eval  speedup | run              us/eval  speedup'
  DO iTol = 1, nTol
    WRITE(*,'(ES10.1,A,A,A,A)') Tol(iTol), ' |', TolEntry( Cheapest( MaxC02, Tol(iTol) ) ), ' |', &
        TolEntry( Cheapest( P99C02, Tol(iTol) ) )
  END DO

  IF( WriteSamples ) CALL WriteSampleFile

  DEALLOCATE( Thermo, PointRow, Val, Sig, Work )

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


  CHARACTER(LEN=15) FUNCTION RunLabel( nn )
    ! e.g. 'B Momentum 50'
    INTEGER, INTENT(IN) :: nn
    WRITE(RunLabel,'(A1,1X,A9,I4)') SchTag(RunSch(nn)), SchName(RunSch(nn)), RunNE(nn)
  END FUNCTION RunLabel


  CHARACTER(LEN=10) FUNCTION RunLabelShort( nn )
    ! e.g. 'B256'
    INTEGER, INTENT(IN) :: nn
    WRITE(RunLabelShort,'(A1,I0)') SchTag(RunSch(nn)), RunNE(nn)
  END FUNCTION RunLabelShort


  CHARACTER(LEN=1) FUNCTION HiMark( nn )
    INTEGER, INTENT(IN) :: nn
    HiMark = ' '
    IF( RunIsHi(nn) ) HiMark = '*'
    IF( nn == iRef  ) HiMark = 'R'
  END FUNCTION HiMark


  SUBROUTINE ComputeRun( nn )

    INTEGER, INTENT(IN) :: nn

    REAL(DP) :: Mul, ml
    INTEGER  :: j, ip, il, iSc, nQ, count_start, count_end, count_rate

    iSc = RunSch(nn)
    nQ  = RunNE(nn)

    ! Quadrature rules of the scheme at resolution nQ, computed once (not
    ! timed); scheme index = method index (1 Energy, 2 Momentum, 3 Tail)
    SELECT CASE( iSc )
    CASE( 1 )
      CALL InitOpacity_CC_2D( nQ, iEmAb2D_Energy )
    CASE( 3 )
      CALL InitOpacity_CC_2D( nQ, iEmAb2D_Tail )
    CASE DEFAULT
      CALL InitOpacity_CC_2D( nQ, iEmAb2D_Momentum )
    END SELECT

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
      CALL Opacity_CC_2D( CorrOf(j), ReacOf(j), EnuA(EnuOf(j)), Val(j,nn), &
          Thermo(1,ip), Mul, Thermo(7,ip), Thermo(8,ip), ml, Thermo(11,ip), Thermo(12,ip), &
          Thermo(9,ip), Thermo(10,ip) )
    END DO
#if defined(WEAKLIB_OMP)
    !$OMP END PARALLEL DO
#endif
    CALL SYSTEM_CLOCK( count_end )
    Wall(nn) = DBLE( count_end - count_start ) / DBLE( count_rate )

    WRITE(*,'(A,A15,A,F9.2,A)') '   computed ', RunLabel(nn), '  in', Wall(nn), ' s'

  END SUBROUTINE ComputeRun


  SUBROUTINE SetSignificance

    ! max over energies per (point, lepton, reaction, correction) and max
    ! per (lepton, reaction, correction) channel, of the reference
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


  REAL(DP) FUNCTION RelErr( j, na, nb )
    ! |Val(na) - Val(nb)| / ref; zero or NaN in either counts as 1
    INTEGER, INTENT(IN) :: j, na, nb
    REAL(DP) :: va, vb
    va = Val(j,na)
    vb = Val(j,nb)
    IF( va /= va .OR. va == 0.0d0 .OR. vb /= vb .OR. vb == 0.0d0 ) THEN
      RelErr = 1.0d0
    ELSE
      RelErr = ABS( va - vb ) / Val(j,iRef)
    END IF
  END FUNCTION RelErr


  LOGICAL FUNCTION InGroup( j, iLepSel, iReacSel, cMin, cMax )
    ! iLepSel, iReacSel < 0: all
    INTEGER, INTENT(IN) :: j, iLepSel, iReacSel, cMin, cMax
    InGroup = Sig(j) .AND. CorrOf(j) >= cMin .AND. CorrOf(j) <= cMax
    IF( iLepSel  > 0 ) InGroup = InGroup .AND. LepOf(j)  == iLepSel
    IF( iReacSel > 0 ) InGroup = InGroup .AND. ReacOf(j) == iReacSel
  END FUNCTION InGroup


  SUBROUTINE GroupStats( na, nb, iLepSel, iReacSel, cMin, cMax, nG, gMax, gP99, gMed )

    ! stats of RelErr( ., na, nb ) over the significant samples of a group
    INTEGER , INTENT(IN)  :: na, nb, iLepSel, iReacSel, cMin, cMax
    INTEGER , INTENT(OUT) :: nG
    REAL(DP), INTENT(OUT) :: gMax, gP99, gMed
    INTEGER :: j

    nG = 0
    DO j = 1, nS
      IF( .NOT. InGroup( j, iLepSel, iReacSel, cMin, cMax ) ) CYCLE
      nG = nG + 1
      Work(nG) = RelErr( j, na, nb )
    END DO
    IF( nG == 0 ) THEN
      gMax = 0.0d0
      gP99 = 0.0d0
      gMed = 0.0d0
      RETURN
    END IF
    CALL HeapSort( Work, nG )
    gMax = Work(nG)
    gP99 = Work( MAX( 1, CEILING( 0.99d0*nG ) ) )
    IF( MOD(nG,2) == 1 ) THEN
      gMed = Work((nG+1)/2)
    ELSE
      gMed = 0.5d0*( Work(nG/2) + Work(nG/2+1) )
    END IF

  END SUBROUTINE GroupStats


  REAL(DP) FUNCTION GroupMax( na, nb, iLepSel, iReacSel, cMin, cMax )
    INTEGER, INTENT(IN) :: na, nb, iLepSel, iReacSel, cMin, cMax
    INTEGER :: j
    GroupMax = 0.0d0
    DO j = 1, nS
      IF( .NOT. InGroup( j, iLepSel, iReacSel, cMin, cMax ) ) CYCLE
      GroupMax = MAX( GroupMax, RelErr( j, na, nb ) )
    END DO
  END FUNCTION GroupMax


  REAL(DP) FUNCTION SpeedupOf( nn )
    INTEGER, INTENT(IN) :: nn
    IF( iBase > 0 .AND. UsEval(nn) > 0.0d0 ) THEN
      SpeedupOf = UsEval(iBase) / UsEval(nn)
    ELSE
      SpeedupOf = 0.0d0
    END IF
  END FUNCTION SpeedupOf


  INTEGER FUNCTION Cheapest( err, tl )
    ! cheapest run (excluding the reference) with err <= tl; 0 if none
    REAL(DP), INTENT(IN) :: err(nRun), tl
    INTEGER :: nn
    Cheapest = 0
    DO nn = 1, nRun
      IF( nn == iRef .OR. err(nn) > tl ) CYCLE
      IF( Cheapest == 0 ) THEN
        Cheapest = nn
      ELSE IF( UsEval(nn) < UsEval(Cheapest) ) THEN
        Cheapest = nn
      END IF
    END DO
  END FUNCTION Cheapest


  CHARACTER(LEN=34) FUNCTION TolEntry( nn )
    INTEGER, INTENT(IN) :: nn
    IF( nn == 0 ) THEN
      WRITE(TolEntry,'(A)') ' none'
    ELSE
      WRITE(TolEntry,'(1X,A15,F9.2,F9.2)') RunLabel(nn), UsEval(nn), SpeedupOf(nn)
    END IF
  END FUNCTION TolEntry


  SUBROUTINE WriteSampleFile

    INTEGER :: j, nn
    CHARACTER(LEN=32) :: Fmt

    OPEN(UNIT=124, FILE='EmAb2DCompareSchemes.dat', STATUS='REPLACE', ACTION='WRITE')
    WRITE(124,'(A)') '# significant samples only; rel. error vs reference ' // TRIM(RunLabel(iRef))
    WRITE(124,'(A)') '# units of ref: cm^-1 (reactions 1-3), (s cm^3 MeV)^-1 (reaction 4)'
    WRITE(124,'(A,*(1X,A10))') '# row lep reac corr Enu[MeV] T[MeV] ref   err for runs:', &
        ( RunLabelShort(nn), nn = 1, nRun )
    WRITE(Fmt,'(A,I0,A)') '(I6,3I4,2ES13.5,ES14.6,', nRun, 'ES11.3)'
    DO j = 1, nS
      IF( .NOT. Sig(j) ) CYCLE
      WRITE(124,Fmt) PointRow(PointOf(j)), LepOf(j), ReacOf(j), CorrOf(j), &
          EnuA(EnuOf(j)), Thermo(1,PointOf(j)), Val(j,iRef), ( RelErr( j, nn, iRef ), nn = 1, nRun )
    END DO
    CLOSE(124)
    WRITE(*,'(/,A)') ' Per-sample errors written to EmAb2DCompareSchemes.dat'

  END SUBROUTINE WriteSampleFile


  SUBROUTINE HeapSort( a, na )

    REAL(DP), INTENT(INOUT) :: a(:)
    INTEGER , INTENT(IN)    :: na
    INTEGER  :: ii, iEnd
    REAL(DP) :: tmp

    DO ii = na/2, 1, -1
      CALL SiftDown( a, ii, na )
    END DO
    DO iEnd = na, 2, -1
      tmp     = a(1)
      a(1)    = a(iEnd)
      a(iEnd) = tmp
      CALL SiftDown( a, 1, iEnd-1 )
    END DO

  END SUBROUTINE HeapSort


  SUBROUTINE SiftDown( a, iStart, iEnd )

    REAL(DP), INTENT(INOUT) :: a(:)
    INTEGER , INTENT(IN)    :: iStart, iEnd
    INTEGER  :: root, child
    REAL(DP) :: tmp

    root = iStart
    DO
      child = 2*root
      IF( child > iEnd ) EXIT
      IF( child < iEnd ) THEN
        IF( a(child) < a(child+1) ) child = child + 1
      END IF
      IF( a(root) >= a(child) ) EXIT
      tmp      = a(root)
      a(root)  = a(child)
      a(child) = tmp
      root = child
    END DO

  END SUBROUTINE SiftDown

END PROGRAM wlTestEmAb2DCompareSchemes
