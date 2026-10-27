PROGRAM wlTestEmAb2DConvergence

  !=========================================================================
  !  Convergence of the semi-leptonic charged-current 2D kernels
  !  (OpacitySource/GeneralEmAb) with the resolution parameter nE of one
  !  integration method, iMethod (default: iEmAb2D_Default, i.e. the method
  !  Opacity_CC_2D uses in this build; set it to iEmAb2D_Energy /
  !  _Momentum / _Tail to test another one). The meaning of nE depends on
  !  the method:
  !    ENERGY   ("standard"): nE x nE linear Gauss-Legendre nodes in E2, E3
  !    MOMENTUM (default)   : nE outer p2 nodes per outer segment,
  !                           MAX(6,MIN(12,nE/4)) inner p3 nodes per segment
  !    TAIL                 : nE nodes per piece, Gauss-Laguerre tail
  !  The nE lists and the benchmark nE_Ref depend on the method (the cost of
  !  ENERGY and TAIL grows like nE**2). For every nE in the list this test
  !  calls InitOpacity_CC_2D(nE, iMethod) (nodes built once per nE, not
  !  timed) and computes, with Opacity_CC_2D,
  !
  !    - every PointStride-th row of ThermoConditions.dat
  !    - nEnu log-spaced neutrino energies in [EnuMin, EnuMax]
  !    - electrons (Mue, me) and muons (Mum, mmu)
  !    - reactions 1-4 and corrections 0-3
  !
  !  and compares with the same method at the benchmark nE = nE_Ref.
  !
  !  This is self-convergence: errors common to all nE of the method (for
  !  ENERGY and MOMENTUM the thermal cutoff Tfac = 30, a floor of a few
  !  1e-5 max) are not visible here; wlTestEmAb2DCompareSchemes compares
  !  the methods against the Tfac-free TAIL reference.
  !
  !  Relative errors |v - ref|/|ref| are evaluated on "significant"
  !  samples only:
  !    ref > 0,
  !    ref >= SigFrac  * max over energies at the same (point, lepton,
  !                      reaction, correction), and
  !    ref >= ChanFrac * max of the (lepton, reaction, correction)
  !                      channel over all points and energies
  !  which excludes the deep Boltzmann tails. Samples where v = 0 or v is
  !  NaN while ref is significant are counted and reported (a zero counts
  !  as relative error 1).
  !
  !  Correction 3 (q^2 dependence of the form factors) is reported
  !  separately: the form-factor term is dropped when it would make |M|^2
  !  negative, which makes the integrand discontinuous and limits its
  !  convergence at high energies.
  !
  !  Units: cm^-1 (reactions 1-3) and (s cm^3 MeV)^-1 (reaction 4); only
  !  relative errors are reported.
  !
  !  Output: summary on stdout and, if WriteSamples, one line per
  !  significant sample in EmAb2DConvergence.dat (indices, Enu, T, ref and
  !  the relative error for every nE in nE_List).
  !=========================================================================

  USE wlKindModule, ONLY: dp
  USE wlSemiLeptonicOpacityModule2D, ONLY: &
    Opacity_CC_2D, InitOpacity_CC_2D, &
    iEmAb2D_Default, iEmAb2D_Energy, iEmAb2D_Momentum, iEmAb2D_Tail
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
  INTEGER , PARAMETER :: iMethod = iEmAb2D_Default   ! integration method tested
  ! nE lists, benchmark and speedup base per method
  INTEGER , PARAMETER :: nE_List_Momentum(13) = &
    (/ 6, 8, 10, 12, 16, 20, 24, 32, 40, 48, 50, 64, 96 /)
  INTEGER , PARAMETER :: nE_List_Energy(13) = &
    (/ 8, 12, 16, 20, 24, 32, 40, 48, 50, 64, 80, 96, 128 /)
  INTEGER , PARAMETER :: nE_List_Tail(6) = &
    (/ 4, 6, 8, 10, 12, 16 /)
  INTEGER , PARAMETER :: nE_Ref_Momentum = 256, nE_Base_Momentum = 50
  INTEGER , PARAMETER :: nE_Ref_Energy   = 256, nE_Base_Energy   = 50
  INTEGER , PARAMETER :: nE_Ref_Tail     =  24, nE_Base_Tail     =  8
  REAL(DP), PARAMETER :: SigFrac  = 1.0d-6
  REAL(DP), PARAMETER :: ChanFrac = 1.0d-15
  INTEGER , PARAMETER :: nTol = 5
  REAL(DP), PARAMETER :: Tol(nTol) = (/ 1.0d-2, 3.0d-3, 1.0d-3, 3.0d-4, 1.0d-4 /)
  REAL(DP), PARAMETER :: TolWorst = 1.0d-3      ! worst samples listed at the smallest nE meeting this
  INTEGER , PARAMETER :: nWorst = 8
  LOGICAL , PARAMETER :: WriteSamples = .TRUE.
  CHARACTER(LEN=3), PARAMETER :: LepName(nLep) = (/ 'e  ', 'mu ' /)

  INTEGER  :: nThermoPoints, nPoints, nS, nThreads, nTest, nE_Ref, nE_Base
  INTEGER , ALLOCATABLE :: nE_List(:)
  CHARACTER(LEN=8) :: MethodName
  INTEGER  :: i, iRow, iLep, k, l, iS, m, n, iTol, iBase, iW
  REAL(DP) :: EnuA(nEnu), Row(12)
  REAL(DP) :: xMax, xP99, xMed

  REAL(DP), ALLOCATABLE :: Thermo(:,:)              ! (12, nPoints)
  INTEGER , ALLOCATABLE :: PointRow(:)
  REAL(DP), ALLOCATABLE :: Ref(:), Val(:,:)         ! (nS), (nS, nTest)
  LOGICAL , ALLOCATABLE :: Sig(:)
  REAL(DP), ALLOCATABLE :: Work(:)
  REAL(DP) :: t_Ref
  REAL(DP), ALLOCATABLE :: t_Test(:), MaxC02(:), P99C02(:)
  INTEGER , ALLOCATABLE :: nZero(:), nNaN(:)
  INTEGER  :: nSig, nSigC02

  SELECT CASE( iMethod )
  CASE( iEmAb2D_Energy )
    MethodName = 'ENERGY'
    nE_Ref  = nE_Ref_Energy
    nE_Base = nE_Base_Energy
    ALLOCATE( nE_List(SIZE(nE_List_Energy)) )
    nE_List = nE_List_Energy
  CASE( iEmAb2D_Tail )
    MethodName = 'TAIL'
    nE_Ref  = nE_Ref_Tail
    nE_Base = nE_Base_Tail
    ALLOCATE( nE_List(SIZE(nE_List_Tail)) )
    nE_List = nE_List_Tail
  CASE DEFAULT
    MethodName = 'MOMENTUM'
    nE_Ref  = nE_Ref_Momentum
    nE_Base = nE_Base_Momentum
    ALLOCATE( nE_List(SIZE(nE_List_Momentum)) )
    nE_List = nE_List_Momentum
  END SELECT
  nTest = SIZE( nE_List )
  ALLOCATE( t_Test(nTest), MaxC02(nTest), P99C02(nTest), nZero(nTest), nNaN(nTest) )

  ! Log-spaced neutrino energies
  DO l = 1, nEnu
    EnuA(l) = EnuMin * (EnuMax/EnuMin)**( DBLE(l-1) / DBLE(nEnu-1) )
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
  ALLOCATE( Ref(nS), Val(nS, nTest), Sig(nS), Work(nS) )

  nThreads = 1
#if defined(WEAKLIB_OMP)
  nThreads = omp_get_max_threads()
#endif

  WRITE(*,'(A)') ' wlTestEmAb2DConvergence: relative error vs resolution nE'
  WRITE(*,'(A,A,A,I2,A)') ' Integration method: ', TRIM(MethodName), &
      '   (default of Opacity_CC_2D in this build: method', iEmAb2D_Default, ')'
  WRITE(*,'(A,I6,A,I3,A,I3,A,F6.1,A,F6.1,A)') ' Thermo points:', nPoints, '  (stride', PointStride, &
      ')   energies:', nEnu, ' in [', EnuMin, ',', EnuMax, '] MeV'
  WRITE(*,'(A,I9,A,I4,A,I3)') ' Samples:', nS, '   benchmark nE:', nE_Ref, '   threads:', nThreads

  CALL ComputeAll( nE_Ref, Ref, t_Ref )
  DO n = 1, nTest
    CALL ComputeAll( nE_List(n), Val(:,n), t_Test(n) )
  END DO

  CALL SetSignificance

  nSig    = COUNT( Sig )
  nSigC02 = 0
  DO iS = 1, nS
    IF( Sig(iS) .AND. CorrOf(iS) <= 2 ) nSigC02 = nSigC02 + 1
  END DO
  WRITE(*,'(A,I9,A,I9,A)') ' Significant samples:', nSig, '  (', nSigC02, ' with corrections 0-2)'
  WRITE(*,'(A,ES8.1,A,ES8.1,A)') ' sig: ref >', SigFrac, ' * max_E(point) and >', ChanFrac, &
      ' * max(channel)'

  ! Zeros and NaNs on significant samples
  DO n = 1, nTest
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
  END DO

  iBase = 0
  DO n = 1, nTest
    IF( nE_List(n) == nE_Base ) iBase = n
  END DO

  ! ---- main table: corrections 0-2
  WRITE(*,'(/,A)') ' Corrections 0-2 (all leptons): relative error vs benchmark'
  WRITE(*,'(A)')   '   nE |  max       p99       median   |  max reac 1  reac 2    reac 3    reac 4  ' // &
                   '|  zero  NaN | wall[s]  us/eval  speedup'
  DO n = 1, nTest
    CALL GroupStats( n, -1, -1, 0, 2, m, xMax, xP99, xMed )
    MaxC02(n) = xMax
    P99C02(n) = xP99
    WRITE(*,'(I6,A,3ES10.2,A,4ES10.2,A,2I6,A,F8.2,F9.2,F9.2)') nE_List(n), ' |', xMax, xP99, xMed, ' |', &
        ( GroupMax( n, -1, k, 0, 2 ), k = 1, nReac ), ' |', nZero(n), nNaN(n), ' |', &
        t_Test(n), t_Test(n)*nThreads/DBLE(nS)*1.0d6, SpeedupOf(n)
  END DO
  WRITE(*,'(A,I4,A,F8.2,A,F9.2,A)') ' benchmark nE =', nE_Ref, ':  wall', t_Ref, ' s,', &
      t_Ref*nThreads/DBLE(nS)*1.0d6, ' us/eval'
  WRITE(*,'(A,I3,A)') ' (us/eval = wall x threads / samples; speedup relative to nE =', nE_Base, ')'

  ! ---- per lepton, corrections 0-2
  DO iLep = 1, nLep
    WRITE(*,'(/,A,A,A)') ' Corrections 0-2, lepton ', TRIM(LepName(iLep)), ':'
    WRITE(*,'(A)') '   nE |  max       p99       median   |  max reac 1  reac 2    reac 3    reac 4'
    DO n = 1, nTest
      CALL GroupStats( n, iLep, -1, 0, 2, m, xMax, xP99, xMed )
      WRITE(*,'(I6,A,3ES10.2,A,4ES10.2)') nE_List(n), ' |', xMax, xP99, xMed, ' |', &
          ( GroupMax( n, iLep, k, 0, 2 ), k = 1, nReac )
    END DO
  END DO

  ! ---- correction 3
  WRITE(*,'(/,A)') ' Correction 3 (form factors; FF clip limits convergence at high E):'
  WRITE(*,'(A)') '   nE |  max       p99       median   |  max reac 1  reac 2    reac 3    reac 4'
  DO n = 1, nTest
    CALL GroupStats( n, -1, -1, 3, 3, m, xMax, xP99, xMed )
    WRITE(*,'(I6,A,3ES10.2,A,4ES10.2)') nE_List(n), ' |', xMax, xP99, xMed, ' |', &
        ( GroupMax( n, -1, k, 3, 3 ), k = 1, nReac )
  END DO

  ! ---- smallest nE meeting each tolerance (corrections 0-2, all reactions and leptons)
  WRITE(*,'(/,A)') ' Smallest nE meeting the tolerance (corrections 0-2, all reactions, both leptons):'
  WRITE(*,'(A)')   '      tol |  by max  by p99'
  DO iTol = 1, nTol
    WRITE(*,'(ES10.1,A,2A8)') Tol(iTol), ' |', SmallestNE( MaxC02, Tol(iTol) ), &
        SmallestNE( P99C02, Tol(iTol) )
  END DO

  ! ---- worst samples at the smallest nE meeting TolWorst by p99
  iW = 0
  DO n = nTest, 1, -1
    IF( P99C02(n) <= TolWorst ) iW = n
  END DO
  IF( iW == 0 ) iW = nTest
  CALL ListWorst( iW )

  IF( WriteSamples ) CALL WriteSampleFile

  DEALLOCATE( Thermo, PointRow, Ref, Val, Sig, Work )
  DEALLOCATE( nE_List, t_Test, MaxC02, P99C02, nZero, nNaN )

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


  SUBROUTINE ComputeAll( nQ, Out, Wall )

    INTEGER , INTENT(IN)  :: nQ
    REAL(DP), INTENT(OUT) :: Out(nS), Wall

    REAL(DP) :: Mul, ml
    INTEGER  :: j, ip, il, count_start, count_end, count_rate

    ! Quadrature rules of iMethod at resolution nQ, computed once (not timed)
    CALL InitOpacity_CC_2D( nQ, iMethod )

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
      CALL Opacity_CC_2D( CorrOf(j), ReacOf(j), EnuA(EnuOf(j)), Out(j), &
          Thermo(1,ip), Mul, Thermo(7,ip), Thermo(8,ip), ml, Thermo(11,ip), Thermo(12,ip), &
          Thermo(9,ip), Thermo(10,ip) )
    END DO
#if defined(WEAKLIB_OMP)
    !$OMP END PARALLEL DO
#endif
    CALL SYSTEM_CLOCK( count_end )
    Wall = DBLE( count_end - count_start ) / DBLE( count_rate )

    WRITE(*,'(A,I4,A,F9.2,A)') '   computed nE =', nQ, '  in', Wall, ' s'

  END SUBROUTINE ComputeAll


  SUBROUTINE SetSignificance

    ! max over energies per (point, lepton, reaction, correction) and max
    ! per (lepton, reaction, correction) channel
    REAL(DP) :: MaxE(nS/nEnu), MaxChan(nCorr, nReac, nLep)
    INTEGER  :: j, b

    MaxE    = 0.0d0
    MaxChan = 0.0d0
    DO j = 1, nS
      b = (j-1)/nEnu + 1
      IF( Ref(j) > 0.0d0 ) THEN
        MaxE(b) = MAX( MaxE(b), Ref(j) )
        MaxChan(CorrOf(j)+1, ReacOf(j), LepOf(j)) = &
          MAX( MaxChan(CorrOf(j)+1, ReacOf(j), LepOf(j)), Ref(j) )
      END IF
    END DO
    DO j = 1, nS
      b = (j-1)/nEnu + 1
      Sig(j) = Ref(j) > 0.0d0 .AND. Ref(j) >= SigFrac*MaxE(b) .AND. &
               Ref(j) >= ChanFrac*MaxChan(CorrOf(j)+1, ReacOf(j), LepOf(j))
    END DO

  END SUBROUTINE SetSignificance


  REAL(DP) FUNCTION RelErr( j, nn )
    ! zero or NaN on a significant sample count as relative error 1
    INTEGER, INTENT(IN) :: j, nn
    IF( Val(j,nn) /= Val(j,nn) .OR. Val(j,nn) == 0.0d0 ) THEN
      RelErr = 1.0d0
    ELSE
      RelErr = ABS( Val(j,nn) - Ref(j) ) / Ref(j)
    END IF
  END FUNCTION RelErr


  LOGICAL FUNCTION InGroup( j, iLepSel, iReacSel, cMin, cMax )
    ! iLepSel, iReacSel < 0: all
    INTEGER, INTENT(IN) :: j, iLepSel, iReacSel, cMin, cMax
    InGroup = Sig(j) .AND. CorrOf(j) >= cMin .AND. CorrOf(j) <= cMax
    IF( iLepSel  > 0 ) InGroup = InGroup .AND. LepOf(j)  == iLepSel
    IF( iReacSel > 0 ) InGroup = InGroup .AND. ReacOf(j) == iReacSel
  END FUNCTION InGroup


  SUBROUTINE GroupStats( nn, iLepSel, iReacSel, cMin, cMax, nG, gMax, gP99, gMed )

    INTEGER , INTENT(IN)  :: nn, iLepSel, iReacSel, cMin, cMax
    INTEGER , INTENT(OUT) :: nG
    REAL(DP), INTENT(OUT) :: gMax, gP99, gMed
    INTEGER :: j

    nG = 0
    DO j = 1, nS
      IF( .NOT. InGroup( j, iLepSel, iReacSel, cMin, cMax ) ) CYCLE
      nG = nG + 1
      Work(nG) = RelErr( j, nn )
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


  REAL(DP) FUNCTION GroupMax( nn, iLepSel, iReacSel, cMin, cMax )
    INTEGER, INTENT(IN) :: nn, iLepSel, iReacSel, cMin, cMax
    INTEGER :: j
    GroupMax = 0.0d0
    DO j = 1, nS
      IF( .NOT. InGroup( j, iLepSel, iReacSel, cMin, cMax ) ) CYCLE
      GroupMax = MAX( GroupMax, RelErr( j, nn ) )
    END DO
  END FUNCTION GroupMax


  REAL(DP) FUNCTION SpeedupOf( nn )
    INTEGER, INTENT(IN) :: nn
    IF( iBase > 0 .AND. t_Test(nn) > 0.0d0 ) THEN
      SpeedupOf = t_Test(iBase) / t_Test(nn)
    ELSE
      SpeedupOf = 0.0d0
    END IF
  END FUNCTION SpeedupOf


  CHARACTER(LEN=8) FUNCTION SmallestNE( err, tl )
    ! smallest nE from which on (all larger nE included) err <= tl
    REAL(DP), INTENT(IN) :: err(nTest), tl
    INTEGER :: nn, nFirst
    nFirst = 0
    DO nn = nTest, 1, -1
      IF( err(nn) > tl ) EXIT
      nFirst = nn
    END DO
    IF( nFirst == 0 ) THEN
      WRITE(SmallestNE,'(A8)') '>' // TRIM(IntStr(nE_List(nTest)))
    ELSE
      WRITE(SmallestNE,'(I8)') nE_List(nFirst)
    END IF
  END FUNCTION SmallestNE


  CHARACTER(LEN=8) FUNCTION IntStr( iv )
    INTEGER, INTENT(IN) :: iv
    WRITE(IntStr,'(I0)') iv
  END FUNCTION IntStr


  SUBROUTINE ListWorst( nn )

    INTEGER, INTENT(IN) :: nn
    INTEGER  :: j, q, r, Idx(nWorst)
    REAL(DP) :: Err(nWorst), e

    Err = -1.0d0
    Idx = 0
    DO j = 1, nS
      IF( .NOT. InGroup( j, -1, -1, 0, 2 ) ) CYCLE
      e = RelErr( j, nn )
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

    WRITE(*,'(/,A,I4,A)') ' Worst samples (corrections 0-2) at nE =', nE_List(nn), ':'
    WRITE(*,'(A)') '   row lep reac corr   Enu[MeV]    T[MeV]  Rho[g/cm3]     ref           rel.err'
    DO q = 1, nWorst
      j = Idx(q)
      IF( j == 0 ) EXIT
      WRITE(*,'(I6,A4,I5,I5,2F11.3,ES12.3,ES14.4,ES12.3)') PointRow(PointOf(j)), &
          TRIM(LepName(LepOf(j))), ReacOf(j), CorrOf(j), EnuA(EnuOf(j)), &
          Thermo(1,PointOf(j)), Thermo(2,PointOf(j)), Ref(j), Err(q)
    END DO

  END SUBROUTINE ListWorst


  SUBROUTINE WriteSampleFile

    INTEGER :: j, nn
    CHARACTER(LEN=32) :: Fmt

    OPEN(UNIT=124, FILE='EmAb2DConvergence.dat', STATUS='REPLACE', ACTION='WRITE')
    WRITE(124,'(A)') '# significant samples only; rel. error vs benchmark nE = ' // TRIM(IntStr(nE_Ref))
    WRITE(124,'(A)') '# units of ref: cm^-1 (reactions 1-3), (s cm^3 MeV)^-1 (reaction 4)'
    WRITE(124,'(A,*(I11))') '# row lep reac corr Enu[MeV] T[MeV] ref   err(nE) for nE =', nE_List
    WRITE(Fmt,'(A,I0,A)') '(I6,3I4,2ES13.5,ES14.6,', nTest, 'ES11.3)'
    DO j = 1, nS
      IF( .NOT. Sig(j) ) CYCLE
      WRITE(124,Fmt) PointRow(PointOf(j)), LepOf(j), ReacOf(j), CorrOf(j), &
          EnuA(EnuOf(j)), Thermo(1,PointOf(j)), Ref(j), ( RelErr( j, nn ), nn = 1, nTest )
    END DO
    CLOSE(124)
    WRITE(*,'(/,A)') ' Per-sample errors written to EmAb2DConvergence.dat'

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

END PROGRAM wlTestEmAb2DConvergence
