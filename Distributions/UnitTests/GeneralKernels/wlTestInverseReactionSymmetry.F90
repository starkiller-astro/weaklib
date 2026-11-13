PROGRAM wlTestInverseReactionSymmetry

  ! Checks the relations between each leptonic reaction and its inverse
  ! (Bollig 2018 thesis), for processes 1..36:
  !   NES/NMS  1..24 : Phin_p(E1,E3)   = Phout_p(E3,E1)          (Sec. 9.1)
  !   LFE/LFC 25..32 : Phin_fwd(E1,E3) = Phout_rev(E3,E1)        (Eqs. 6.40-6.43)
  !                    Phout_fwd(E1,E3)= Phin_rev(E3,E1)
  !                    pairs (25,26) (27,28) (29,30) (31,32)
  !   IMD     33..36 : Phout/in_34(E1,E3) = Phout/in_33(E3,E1)   (Eq. 6.59)
  !                    and 36 vs 35
  ! for l = 0, 1, on both the interface (CalculateAllPhout_3D) and the
  ! reference (CalculatePhoutPhin) path. Then, with arbitrary non-equilibrium
  ! neutrino spectra, the discrete number-conservation error these relations
  ! control: dN_a + dN_b = 0 for 2->2 (dN = 0 elastic), dN_1 = dN_3 for IMD
  ! (Bollig Eq. 9.45).
  ! Finally, the thornado rate assembly, where in-rates use a detailed balance
  ! factor built from the (DG-limited) equilibrium distributions J0: storing
  ! both directions of a pair conserves number only if J0 is the exact FD,
  ! while building the partner from the transposed stored row always does.

  USE wlKindModule, ONLY: dp
  USE wlEosConstantsModule, ONLY: ElectronMass => me, MuonMass => mmu
  USE wlGeneralLeptonScatteringModule, ONLY: &
    CalculatePhoutPhin
  USE wlGeneralLeptonScatteringModuleThornadoInterface, ONLY: &
    InitGeneralScatteringKernels, &
    FinalizeGeneralScatteringKernels, &
    CalculateAllPhout_3D, &
    DetailedBalanceFactor, &
    ProcessRowIndex

  IMPLICIT NONE

  INTEGER , PARAMETER :: nE = 16, nTh = 24, nL = 2, nX = 5, nP = 36
  INTEGER , PARAMETER :: nPair = 4, nIMD = 2
  INTEGER , PARAMETER :: PairFwd(nPair) = (/ 25, 27, 29, 31 /)
  INTEGER , PARAMETER :: PairRev(nPair) = (/ 26, 28, 30, 32 /)
  INTEGER , PARAMETER :: IMD1(nIMD)     = (/ 33, 35 /)   ! nu_1 = nu_bar_e / nu_e
  INTEGER , PARAMETER :: IMD3(nIMD)     = (/ 34, 36 /)   ! the same with nu_1 <-> nu_3
  REAL(DP), PARAMETER :: dlnE     = LOG( 1.3d0 )
  REAL(DP), PARAMETER :: RelFloor = 1.0d-10  ! ignore entries below this * max|ref|
  REAL(DP), PARAMETER :: AbsFloor = 1.0d-10  ! [cm^-1] ignore entries with Phi*E3^3 below this
  REAL(DP), PARAMETER :: L1Reg    = 1.0d-6   ! l=1: error relative to max(|ref|, L1Reg*max|ref|)
  REAL(DP), PARAMETER :: RelTol   = 1.0d-6
  REAL(DP), PARAMETER :: ConsTol  = 1.0d-6

  REAL(DP) :: E(nE), w(nE), fA(nE), fB(nE), T(nX), Me(nX), Mm(nX)
  REAL(DP) :: PhOut(nE,nE,nP,nL,nX), PhIn(nE,nE,nP,nL,nX)
  REAL(DP) :: RfOut(nE,nE,nP,nL,nX), RfIn(nE,nE,nP,nL,nX)
  REAL(DP) :: Phout_r(nL), Phin_r(nL), xMu2, xMu4, m2, m4
  INTEGER  :: i, j, iX, iP, iRow
  LOGICAL  :: AllOK

  DO i = 1, nE
    E(i) = 2.0d0 * 1.3d0**(i-1)
    w(i) = E(i)**3 * dlnE            ! E^2 dE on the log grid
    ! Arbitrary, non-thermal spectra for the two neutrinos of a pair
    fA(i) = 0.7d0 / ( EXP( (E(i) - 25.0d0) / 6.0d0 ) + 1.0d0 )
    fB(i) = 0.2d0 * EXP( -E(i) / 30.0d0 ) &
          + 0.3d0 / ( EXP( (E(i) - 60.0d0) / 10.0d0 ) + 1.0d0 )
  END DO
  ! (T, Mu_e, Mu_mu) [MeV]; point 4 is hot with mu_mu well below m_mu,
  ! point 5 cold and degenerate.
  T  = (/  5.0d0, 15.2d0,  30.0d0, 40.0d0,   8.0d0 /)
  Me = (/ 20.0d0, 44.4d0, 120.0d0, 10.0d0, 180.0d0 /)
  Mm = (/105.0d0, 37.3d0, 150.0d0, 20.0d0, 175.0d0 /)

  CALL InitGeneralScatteringKernels( E, E, nE, nE, nTh, &
         iProcessList_Option = (/ (iP, iP = 1, nP) /) )

  AllOK = .TRUE.

  ! ------------------------------------------- interface and reference kernels
  CALL CalculateAllPhout_3D( nX, nE, nE, nP, nL, T, Me, Mm, PhOut )

  DO iX = 1, nX
    DO iP = 1, nP
      iRow = ProcessRowIndex( iP )
      CALL ProcessLeptons( iP, Me(iX), Mm(iX), xMu2, xMu4, m2, m4 )
      DO j = 1, nE
        DO i = 1, nE
          PhIn(i,j,iP,:,iX) = PhOut(i,j,iRow,:,iX) &
            * DetailedBalanceFactor( iRow, E(i), E(j), T(iX), Me(iX), Mm(iX) )
          CALL CalculatePhoutPhin( E(i), E(j), T(iX), xMu2, xMu4, m2, m4, &
                                   iP, Phout_r, Phin_r, nL, nTh )
          RfOut(i,j,iP,:,iX) = Phout_r
          RfIn (i,j,iP,:,iX) = Phin_r
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  ! Rows are in process order (ProcList = 1..nP); reorder defensively anyway
  PhOut = PhOut(:,:,(/ (ProcessRowIndex(iP), iP = 1, nP) /),:,:)

  CALL RunChecks( 'interface', PhOut, PhIn  )
  CALL RunChecks( 'reference', RfOut, RfIn  )

  CALL ThornadoAssembly( PhOut )

  CALL FinalizeGeneralScatteringKernels()

  IF ( AllOK ) THEN
    WRITE(*,*) 'wlTestInverseReactionSymmetry: PASSED'
  ELSE
    WRITE(*,*) 'wlTestInverseReactionSymmetry: FAILED'
  ENDIF

CONTAINS

  ! Signed chemical potentials and masses of l2, l4 for an absolute process
  ! index, as assigned in the interface (SelectiDistinctFromProcessIndex).
  SUBROUTINE ProcessLeptons( iProc, Mu_e, Mu_mu, xMu2, xMu4, m2, m4 )
    INTEGER , INTENT(IN)  :: iProc
    REAL(DP), INTENT(IN)  :: Mu_e, Mu_mu
    REAL(DP), INTENT(OUT) :: xMu2, xMu4, m2, m4
    REAL(DP) :: s

    s = +1.0d0
    SELECT CASE (iProc)
    CASE (1:12)                       ! odd: e-, even: e+
      IF ( MOD(iProc,2) == 0 ) s = -1.0d0
      xMu2 = s*Mu_e;  m2 = ElectronMass; xMu4 = s*Mu_e;  m4 = ElectronMass
    CASE (13:24)                      ! odd: mu-, even: mu+
      IF ( MOD(iProc,2) == 0 ) s = -1.0d0
      xMu2 = s*Mu_mu; m2 = MuonMass;     xMu4 = s*Mu_mu; m4 = MuonMass
    CASE (25, 32, 33, 34)             ! e- -> mu-
      xMu2 = Mu_e;    m2 = ElectronMass; xMu4 = Mu_mu;   m4 = MuonMass
    CASE (26, 31)                     ! mu- -> e-
      xMu2 = Mu_mu;   m2 = MuonMass;     xMu4 = Mu_e;    m4 = ElectronMass
    CASE (27, 30, 35, 36)             ! e+ -> mu+
      xMu2 = -Mu_e;   m2 = ElectronMass; xMu4 = -Mu_mu;  m4 = MuonMass
    CASE (28, 29)                     ! mu+ -> e+
      xMu2 = -Mu_mu;  m2 = MuonMass;     xMu4 = -Mu_e;   m4 = ElectronMass
    CASE DEFAULT
      WRITE(*,*) 'ProcessLeptons: unknown process ', iProc
      STOP
    END SELECT
  END SUBROUTINE ProcessLeptons

  SUBROUTINE RunChecks( Path, KOut, KIn )
    CHARACTER(*), INTENT(IN) :: Path
    REAL(DP)    , INTENT(IN) :: KOut(nE,nE,nP,nL,nX), KIn(nE,nE,nP,nL,nX)
    REAL(DP) :: ConsErr, ConsMax
    INTEGER  :: iP, k, iX

    WRITE(*,*)
    WRITE(*,'(A,A)') ' ==== ', Path
    WRITE(*,'(A)') '   pointwise: A(E1,E3) vs B(E3,E1)          max rel err  l=0       l=1'

    DO iP = 1, 24
      CALL CheckTransposed( 'Phin  vs Phout', iP, iP, KIn(:,:,iP,:,:), KOut(:,:,iP,:,:) )
    ENDDO
    DO k = 1, nPair
      CALL CheckTransposed( 'Phin  vs Phout', PairFwd(k), PairRev(k), &
                            KIn (:,:,PairFwd(k),:,:), KOut(:,:,PairRev(k),:,:) )
      CALL CheckTransposed( 'Phout vs Phin ', PairFwd(k), PairRev(k), &
                            KOut(:,:,PairFwd(k),:,:), KIn (:,:,PairRev(k),:,:) )
    ENDDO
    DO k = 1, nIMD
      CALL CheckTransposed( 'Phout vs Phout', IMD3(k), IMD1(k), &
                            KOut(:,:,IMD3(k),:,:), KOut(:,:,IMD1(k),:,:) )
      CALL CheckTransposed( 'Phin  vs Phin ', IMD3(k), IMD1(k), &
                            KIn (:,:,IMD3(k),:,:), KIn (:,:,IMD1(k),:,:) )
    ENDDO

    WRITE(*,'(A)') '   number conservation, non-equilibrium spectra   max |net| / gross'
    DO iP = 1, 24
      ConsMax = 0.0d0
      DO iX = 1, nX
        ConsMax = MAX( ConsMax, Cons2to2( KOut(:,:,iP,1,iX), KIn(:,:,iP,1,iX), &
                                          KOut(:,:,iP,1,iX), KIn(:,:,iP,1,iX), fA, fA ) )
      ENDDO
      CALL ReportCons( 'dN          ', iP, iP, ConsMax )
    ENDDO
    DO k = 1, nPair
      ConsMax = 0.0d0
      DO iX = 1, nX
        ConsErr = Cons2to2( KOut(:,:,PairFwd(k),1,iX), KIn(:,:,PairFwd(k),1,iX), &
                            KOut(:,:,PairRev(k),1,iX), KIn(:,:,PairRev(k),1,iX), fA, fB )
        ConsMax = MAX( ConsMax, ConsErr )
      ENDDO
      CALL ReportCons( 'dN_a + dN_b ', PairFwd(k), PairRev(k), ConsMax )
    ENDDO
    DO k = 1, nIMD
      ConsMax = 0.0d0
      DO iX = 1, nX
        ConsErr = ConsIMD( KOut(:,:,IMD1(k),1,iX), KIn(:,:,IMD1(k),1,iX), &
                           KOut(:,:,IMD3(k),1,iX), KIn(:,:,IMD3(k),1,iX), fA, fB )
        ConsMax = MAX( ConsMax, ConsErr )
      ENDDO
      CALL ReportCons( 'dN_1 - dN_3 ', IMD1(k), IMD3(k), ConsMax )
    ENDDO

  END SUBROUTINE RunChecks

  ! Max error of A(i,j) against B(j,i) over all points. l=0 is pointwise
  ! relative; l=1 changes sign, so its error is taken relative to
  ! max(|B|, L1Reg*max|B|) to stay meaningful at the zero crossings.
  SUBROUTINE CheckTransposed( Label, iProcA, iProcB, A, B )
    CHARACTER(*), INTENT(IN) :: Label
    INTEGER     , INTENT(IN) :: iProcA, iProcB
    REAL(DP)    , INTENT(IN) :: A(nE,nE,nL,nX), B(nE,nE,nL,nX)
    REAL(DP) :: Err(nL), Scale, Ref, Den
    INTEGER  :: iL, iX, i, j

    Err = 0.0d0
    DO iX = 1, nX
      DO iL = 1, nL
        Scale = MAXVAL( ABS( B(:,:,iL,iX) ) )
        IF ( Scale <= 0.0d0 ) CYCLE
        DO j = 1, nE
          DO i = 1, nE
            Ref = B(j,i,iL,iX)
            IF ( ABS(Ref) <= RelFloor * Scale )    CYCLE
            IF ( ABS(Ref) * E(j)**3 <= AbsFloor )  CYCLE
            IF ( iL == 1 ) THEN
              Den = ABS(Ref)
            ELSE
              Den = MAX( ABS(Ref), L1Reg * Scale )
            ENDIF
            Err(iL) = MAX( Err(iL), ABS( A(i,j,iL,iX) - Ref ) / Den )
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    WRITE(*,'(A,A,A,I3,A,I3,A,2ES11.3)') '   ', Label, '  ', iProcA, ' vs', iProcB, &
      '                ', Err
    IF ( ANY( Err > RelTol ) ) THEN
      AllOK = .FALSE.
      WRITE(*,'(A)') '     ^^^ FAIL'
    ENDIF
  END SUBROUTINE CheckTransposed

  ! nu_a(E1) + l2 <-> nu_b(E3) + l4 with the forward row F (rows indexed
  ! (E_a,E_b)) and the reverse row R (rows indexed (E_b,E_a)): net nu_a plus
  ! net nu_b number rate, over the gross rate. Elastic: F = R, fa = fb.
  FUNCTION Cons2to2( FOut, FIn, ROut, RIn, fa, fb ) RESULT( Err )
    REAL(DP), INTENT(IN) :: FOut(nE,nE), FIn(nE,nE), ROut(nE,nE), RIn(nE,nE)
    REAL(DP), INTENT(IN) :: fa(nE), fb(nE)
    REAL(DP) :: Err, dNa, dNb, Gross, t1, t2, t3, t4
    INTEGER  :: i, j

    dNa = 0.0d0; dNb = 0.0d0; Gross = 0.0d0
    DO j = 1, nE
      DO i = 1, nE
        t1 = w(i)*w(j) * fa(i) * (1.0d0 - fb(j)) * FOut(i,j)
        t2 = w(i)*w(j) * (1.0d0 - fa(i)) * fb(j) * FIn (i,j)
        t3 = w(i)*w(j) * fb(j) * (1.0d0 - fa(i)) * ROut(j,i)
        t4 = w(i)*w(j) * (1.0d0 - fb(j)) * fa(i) * RIn (j,i)
        dNa   = dNa - t1 + t2
        dNb   = dNb - t3 + t4
        Gross = Gross + t1 + t2 + t3 + t4
      ENDDO
    ENDDO
    Err = 0.0d0
    IF ( Gross > 0.0d0 ) Err = ABS( dNa + dNb ) / Gross
  END FUNCTION Cons2to2

  ! nu_1(E1) + l2 + nu_3(E3) <-> l4, absorption Phout / production Phin
  ! (Bollig Eq. 9.39, l=0). Row P1 has nu_1 as the first neutrino, row P3
  ! the partner as first; both neutrinos must change by the same number.
  FUNCTION ConsIMD( P1Out, P1In, P3Out, P3In, f1, f3 ) RESULT( Err )
    REAL(DP), INTENT(IN) :: P1Out(nE,nE), P1In(nE,nE), P3Out(nE,nE), P3In(nE,nE)
    REAL(DP), INTENT(IN) :: f1(nE), f3(nE)
    REAL(DP) :: Err, dN1, dN3, Gross, t1, t2, t3, t4
    INTEGER  :: i, j

    dN1 = 0.0d0; dN3 = 0.0d0; Gross = 0.0d0
    DO j = 1, nE
      DO i = 1, nE
        t1 = w(i)*w(j) * (1.0d0 - f1(i)) * (1.0d0 - f3(j)) * P1In (i,j)
        t2 = w(i)*w(j) * f1(i) * f3(j)                     * P1Out(i,j)
        t3 = w(i)*w(j) * (1.0d0 - f3(j)) * (1.0d0 - f1(i)) * P3In (j,i)
        t4 = w(i)*w(j) * f3(j) * f1(i)                     * P3Out(j,i)
        dN1   = dN1 + t1 - t2
        dN3   = dN3 + t3 - t4
        Gross = Gross + t1 + t2 + t3 + t4
      ENDDO
    ENDDO
    Err = 0.0d0
    IF ( Gross > 0.0d0 ) Err = ABS( dN1 - dN3 ) / Gross
  END FUNCTION ConsIMD

  ! Emulates thornado's ComputeNeutrinoOpacityRates_MuonLFEC/_MuonIMD (l=0).
  ! J0 is either the FD consistent with the kernel chemistry, or a perturbed
  ! and capped copy standing in for LimitEquilibriumDistributions_DG.
  !   (a) both rows of a pair stored, each with its own J0 detailed balance
  !   (b) only the thornado-stored row; the partner uses its transpose
  SUBROUTINE ThornadoAssembly( KOut )
    REAL(DP), INTENT(IN) :: KOut(nE,nE,nP,nL,nX)
    INTEGER , PARAMETER :: Stored(nPair)  = (/ 25, 27, 30, 32 /)
    INTEGER , PARAMETER :: Partner(nPair) = (/ 26, 28, 29, 31 /)
    REAL(DP), PARAMETER :: J0Cap = 1.0d0 - EPSILON( 1.0d0 )
    REAL(DP) :: J0A(nE), J0B(nE), DB(nE,nE), F(nE,nE), FIn(nE,nE)
    REAL(DP) :: R(nE,nE), RIn(nE,nE), Err(3)
    REAL(DP) :: xMu2, xMu4, m2, m4, MuA, MuB
    INTEGER  :: k, iX, iJ0, i, j

    WRITE(*,*)
    WRITE(*,'(A)') ' ==== thornado-style assembly (interface kernels), max |net| / gross'
    WRITE(*,'(A)') '   stored/partner         (a) exact J0  (a) limited J0  (b) limited J0'

    DO k = 1, nPair + nIMD
      Err = 0.0d0
      DO iX = 1, nX
        IF ( k <= nPair ) THEN
          CALL ProcessLeptons( Stored(k), Me(iX), Mm(iX), xMu2, xMu4, m2, m4 )
        ELSE
          CALL ProcessLeptons( IMD1(k-nPair), Me(iX), Mm(iX), xMu2, xMu4, m2, m4 )
        ENDIF
        DO iJ0 = 1, 2
          ! Chemical equilibrium of the reaction fixes Mu_A -/+ Mu_B
          MuB = 0.5d0 * xMu2
          IF ( k <= nPair ) THEN
            MuA = MuB + xMu4 - xMu2     ! nu_A + l2 <-> nu_B + l4
          ELSE
            MuA = xMu4 - xMu2 - MuB     ! nu_1 + l2 + nu_3 <-> l4
          ENDIF
          DO i = 1, nE
            IF ( iJ0 == 1 ) THEN
              J0A(i) = FDist( E(i), MuA, T(iX) )
              J0B(i) = FDist( E(i), MuB, T(iX) )
            ELSE
              J0A(i) = MIN( FDist( E(i), MuA + 0.05d0*T(iX)*SIN(1.7d0*i), T(iX) ), J0Cap )
              J0B(i) = MIN( FDist( E(i), MuB + 0.05d0*T(iX)*COS(1.3d0*i), T(iX) ), J0Cap )
            ENDIF
          ENDDO

          IF ( k <= nPair ) THEN
            F = KOut(:,:,Stored(k),1,iX)
            DO j = 1, nE
              DO i = 1, nE
                DB(i,j) = ( J0A(i) * ( 1.0d0 - J0B(j) ) ) / ( J0B(j) * ( 1.0d0 - J0A(i) ) )
              ENDDO
            ENDDO
            FIn = F * DB
            ! (a): partner row with its own factor, DB_B(j,i) = 1 / DB_A(i,j)
            R   = KOut(:,:,Partner(k),1,iX)
            RIn = R / TRANSPOSE( DB )
            Err(iJ0) = MAX( Err(iJ0), Cons2to2( F, FIn, R, RIn, fA, fB ) )
            IF ( iJ0 == 2 ) &
              Err(3) = MAX( Err(3), Cons2to2( F, FIn, TRANSPOSE(FIn), TRANSPOSE(F), fA, fB ) )
          ELSE
            ! IMD: Phout = absorption, Phin = production = Phout * DB
            F = KOut(:,:,IMD1(k-nPair),1,iX)
            DO j = 1, nE
              DO i = 1, nE
                DB(i,j) = ( J0A(i) * J0B(j) ) / ( ( 1.0d0 - J0A(i) ) * ( 1.0d0 - J0B(j) ) )
              ENDDO
            ENDDO
            FIn = F * DB
            R   = KOut(:,:,IMD3(k-nPair),1,iX)
            RIn = R * TRANSPOSE( DB )
            Err(iJ0) = MAX( Err(iJ0), ConsIMD( F, FIn, R, RIn, fA, fB ) )
            IF ( iJ0 == 2 ) &
              Err(3) = MAX( Err(3), ConsIMD( F, FIn, TRANSPOSE(F), TRANSPOSE(FIn), fA, fB ) )
          ENDIF
        ENDDO
      ENDDO

      IF ( k <= nPair ) THEN
        WRITE(*,'(A,I3,A,I3,A,3ES15.3)') '   ', Stored(k), ' /', Partner(k), '          ', Err
      ELSE
        WRITE(*,'(A,I3,A,I3,A,3ES15.3)') '   ', IMD1(k-nPair), ' /', IMD3(k-nPair), '          ', Err
      ENDIF
      ! (a) with the consistent J0 checks the emulation; (b) is the thornado scheme
      IF ( Err(1) > ConsTol .OR. Err(3) > ConsTol ) THEN
        AllOK = .FALSE.
        WRITE(*,'(A)') '     ^^^ FAIL'
      ENDIF
    ENDDO

  END SUBROUTINE ThornadoAssembly

  PURE REAL(DP) FUNCTION FDist( En, Mu, Tk )
    REAL(DP), INTENT(IN) :: En, Mu, Tk
    FDist = 1.0d0 / ( EXP( (En - Mu) / Tk ) + 1.0d0 )
  END FUNCTION FDist

  SUBROUTINE ReportCons( Label, iProcA, iProcB, Err )
    CHARACTER(*), INTENT(IN) :: Label
    INTEGER     , INTENT(IN) :: iProcA, iProcB
    REAL(DP)    , INTENT(IN) :: Err

    WRITE(*,'(A,A,I3,A,I3,A,ES11.3)') '   ', Label, iProcA, ' /', iProcB, &
      '                         ', Err
    IF ( Err > ConsTol ) THEN
      AllOK = .FALSE.
      WRITE(*,'(A)') '     ^^^ FAIL'
    ENDIF
  END SUBROUTINE ReportCons

END PROGRAM wlTestInverseReactionSymmetry
