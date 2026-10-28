!=======================================================================
!  wlSemiLeptonicIntegration2D_Tail             (TAIL integration method)
!
!  The semi-leptonic CC rates (Guo 2020 formalism) integrated WITHOUT the
!  thermal cutoff Tfac. Same units and particle mapping as Opacity_CC_2D
!  (stateless kernel Opacity_CC_2D_Tail_Nodes); all physics (|M|^2,
!  kinematic windows, occupation factors) comes from
!  wlSemiLeptonicPhysicsModule2D, the inner kinks from E3_Segments of
!  wlSemiLeptonicIntegration2D_Momentum.
!
!  nE (resolution knob): number of Gauss-Legendre nodes per piece, used
!  for BOTH the outer pieces (in p2) and the inner segments (in p3). The
!  number of pieces/segments is set by the kinematics and Fermi edges
!  (see below), so the cost per call is roughly
!    ( nPiece2*nE + nLag ) x ( nSeg3*nE [+ nLag] )  |M|^2 evaluations.
!
!  Design
!  ------
!  Outer integral over the nucleon momentum p2 on [p2min, inf):
!   * p2min is kinematic only (Range_pn; Range_pn_D with the reaction-3
!     branch, which has the kinematics of reaction 4 without its thermal
!     lower cut).
!   * Finite part [E2(p2min), E2s] split at
!       - the kinematic events of the outer integrand (Outer_Events:
!         P2 = Enu, and the P2 at which a pair of collinear roots appears
!         (SQRT-type point), a collinear root reaches E4 = M4 or E3 = M3,
!         or the decay window opens), in closed form;
!       - graded points e +/- (2, 6, 18)*T around the "Fermi edges"
!         e = Mu2 (f2 or 1-f2), the effective edge X of the inner
!         occupation product (captures: X = Mu3 + Mu4 - Enu; decays:
!         X = Mu4 - Enu - MAX(Mu3, M3)) and the threshold E2(p2min).
!     Pieces ending at a kinematic event get a quadratic end clustering
!     (x = u**2 at that end, smoothstep if both ends), which makes the
!     SQRT-type behaviour there analytic; pieces ending at Fermi points
!     are plain Gauss-Legendre.
!   * Tail [E2s, inf) with Gauss-Laguerre, x = (E2 - E2s)/T, nLag = 12.
!     Split point E2s = MAX( E2(p2min) + cThr*T, MAX(Mu2, X) + cEdge*T,
!     last outer event below E2s + LEv*T, plus cThr*T ). It lies past
!     every Fermi edge of the occupation product, so what is left decays
!     like poly(E2)*EXP(-E2/T): Pauli blocking can no longer cancel the
!     exponential (the failure of a split at mu +/- 3T), and the scale is
!     exactly T.
!  Inner integral over the lepton at fixed E2 (Inner_E3):
!   * window from Q0_Window_Envelope (captures: as Ebounds; decays: with
!     q0_cap = HUGE, i.e. without the Q0_FALLBACK_TFAC cap);
!   * split at the kinks from E3_Segments and at graded points around
!     e = Mu3, the E3 at which E4 = Mu4, and the window start; window
!     ends and kinks get the quadratic end clustering (the angle-
!     integrated |M|^2 is SQRT-type there);
!   * decays: if the window is open (s3 = +1 and |W| < S) or extends more
!     than LMask*T past E3s, [E3s, inf) is done by the same Gauss-Laguerre
!     rule (nodes beyond a finite window end are dropped; their weight is
!     < EXP(-LMask)), E3s = MAX( E3a + cThr*T, edges + cEdge*T,
!     last kink + cThr*T ).
!
!  Measured (stride-10 rows of ThermoConditions.dat = 100 states, 10 Enu
!  in 2-300 MeV, e/mu, reactions 1-4; c0-2 = corrections 0-2, relative
!  errors over the significant samples)
!  -------------------------------------------------------------------
!   vs RL1024 (1024x1024 linear GL, Tfac = 150; own accuracy ~1e-6):
!     scheme C nE =  8: c0-2 max 2.4e-5, p99 3.5e-6
!     scheme C nE = 12: c0-2 max 1.4e-6, p99 3.4e-7   (RL1024 floor)
!     scheme C nE = 16: c0-2 max 1.0e-6, p99 3.2e-7   (RL1024 floor)
!     scheme B Opacity_CC_2D nE = 64/128/256: c0-2 max 1.6e-4/4.2e-5/4.2e-5
!   time per evaluation (gfortran -O2, 1 thread, pinned, same sample):
!     scheme C nE = 8/12/16: 2.5e-4/4.8e-4/7.6e-4 s
!     scheme B (momentum) nE = 64/128/256: 1.5e-4/2.8e-4/6.0e-4 s
!   vs a converged self-reference (nE = 48, 4x subdivided inner, nLag = 32):
!     nE = 8: 2.5e-5 | 10: 6.6e-6 | 12: 1.5e-6 | 14: 8.9e-7 | 16: 5.4e-7
!     nE = 20: 1.5e-7
!   The 2-4e-5 floor of scheme B comes from the inner rule (unresolved
!   SQRT-type window ends/kinks, missing outer kinematic events), not
!   from the Tfac truncation itself; both are removed here.
!  Tail rule (tail error only, c0-2 max, same split points):
!     Gauss-Laguerre, tau = T    : n = 12 4.8e-10, n = 16 1.5e-11
!     Gauss-Laguerre, tau from the log-derivative of the occupation
!       product at the split point (clamped to [T, 2T]):
!                                  n = 8 7.5e-6, 12 4.0e-7, 16 4.1e-8
!       (tau = 0.8T is as good as T, 1.25T is 20x worse: the remaining
!       poly(E) growth wants tau <= T, never larger)
!     exp-sinh, x = EXP(pi/2*SINH(t)): n = 24 5.6e-7, 32 7.6e-9,
!       40 2.7e-9, 48 4.3e-10
!     DE for exponential decay, x = EXP(t - EXP(-t)): n = 12 8.4e-7,
!       16 8.2e-9, 24 1.4e-10 (insensitive to tau)
!   All rules are tabulated, so the cost per node is the same and the
!   error per node decides: Laguerre with tau = T reaches ~5e-10 with
!   12 nodes, exp-sinh needs ~48.
!
!  Known limitations
!   * Correction 3 (form-factor block of Calc_Ampsq) has jump
!     discontinuities inside the kink-free segments (the pion-pole branch
!     switches); no rule converges there faster than O(h): c3 max errors
!     are ~1e-2 for this scheme, scheme B and RL1024 alike.
!   * Outer events further than LEv*T past E2s are not resolved (the
!     integrand there is < EXP(-LEv) of its value at E2s).
!=======================================================================
MODULE wlSemiLeptonicIntegration2D_Tail

  USE wlKindModule, ONLY: dp
  USE wlEosConstantsModule, ONLY: pi, Gw_MeV, Vud, hbarc, hbar
  USE wlSemiLeptonicPhysicsModule2D, ONLY: &
    AmpCoeffsType, Calc_AmpCoeffs, Calc_Ampsq, Range_pn, Range_pn_D, FD, &
    Q0_Window_Envelope, opacity_unit_factor
  USE wlSemiLeptonicIntegration2D_Momentum, ONLY: &
    E3_Segments, nSegMax

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Opacity_CC_2D_Tail_Nodes

  ! Gauss-Laguerre rule for INT_0^inf g(x) dx: nodes xLag, weights
  ! wLag = w_i*EXP(x_i) (the EXP(-x) weight is divided out)
  INTEGER,  PARAMETER :: nLag = 12
  REAL(DP), PARAMETER :: xLag(nLag) = (/ &
    0.11572211735802067527d0, 0.61175748451513066539d0, 1.5126102697764187868d0,  &
    2.8337513377435072286d0,  4.5992276394183484846d0,  6.8445254531151773478d0,  &
    9.6213168424568670439d0,  13.00605499330634772d0,   17.116855187462255728d0,  &
    22.15109037939700567d0,   28.487967250984000313d0,  37.099121044466920337d0 /)
  REAL(DP), PARAMETER :: wLag(nLag) = (/ &
    0.29720963604441079779d0, 0.69646298043059723068d0, 1.1077813946157521533d0,  &
    1.5384642390428291248d0,  1.9983276062742427652d0,  2.5007457691008668748d0,  &
    3.0653215182823862362d0,  3.7232891107827716d0,     4.5298140299817350621d0,  &
    5.5972584618353210927d0,  7.2129954609258770082d0,  10.543837461910081112d0 /)

  ! Split points (in units of T)
  REAL(DP), PARAMETER :: cThr  = 6.0d0    ! tail start past thresholds / kinematic events
  REAL(DP), PARAMETER :: cEdge = 3.0d0    ! tail start past the Fermi edges
  REAL(DP), PARAMETER :: LMask = 36.0d0   ! closed decay window: tail rule if longer than this
  REAL(DP), PARAMETER :: LEv   = 20.0d0   ! outer events this far past E2s still move E2s
  REAL(DP), PARAMETER :: KinkSearch = 200.0d0  ! kinks searched in [E3a, E3a + KinkSearch*T]

  ! Graded breakpoints e +/- Grade*T around each Fermi edge
  INTEGER,  PARAMETER :: nGrade = 3
  REAL(DP), PARAMETER :: Grade(nGrade) = (/ 2.0d0, 6.0d0, 18.0d0 /)
  INTEGER,  PARAMETER :: nEdge  = 3

  INTEGER,  PARAMETER :: nEvMax = 21               ! 1 + 5 conditions x 2 branches x 2 roots
  INTEGER,  PARAMETER :: nBpMax = 2 + nEvMax + 2*nEdge*nGrade

CONTAINS

SUBROUTINE Opacity_CC_2D_Tail_Nodes(WhichCorrection, ReactionIndex, Enu, OpaA, xT, xMul, xMun, xMup, &
     xml, xmn, xmp, xUn, xUp, nE, xa, wxa)

  ! Tail method with explicit Gauss-Legendre nodes and weights on [0,1]
  ! (CALL gauleg(0.d0, 1.d0, xa, wxa, nE)); Opacity_CC_2D passes the rule
  ! precomputed by InitOpacity_CC_2D.
  INTEGER , INTENT(IN)  :: ReactionIndex, WhichCorrection, nE
  REAL(DP), INTENT(IN)  :: Enu,xT,xMul,xMun,xMup,xml,xmn,xmp,xUn,xUp
  REAL(DP), INTENT(IN)  :: xa(nE), wxa(nE)
  REAL(DP), INTENT(OUT) :: OpaA
#if defined(WEAKLIB_OMP_OL)
  !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
  !$ACC ROUTINE SEQ
#endif

  INTEGER  :: IncludeCorrections(3)

  SELECT CASE(WhichCorrection)
  CASE (0)
     IncludeCorrections = (/ 0, 0, 0 /)
  CASE (1)
     IncludeCorrections = (/ 1, 0, 0 /)
  CASE (2)
     IncludeCorrections = (/ 1, 1, 0 /)
  CASE (3)
     IncludeCorrections = (/ 1, 1, 1 /)
  CASE DEFAULT
     ! Invalid input: no I/O or STOP in device-callable code; return zero
     OpaA = 0.0d0
     RETURN
  END SELECT

  SELECT CASE(ReactionIndex)
  CASE (1)
    CALL Integral_2D_Tail(Enu, xT, xmn, xml, xmp, xUn, xUp, xMun,  xMul, xMup,  1.0d0, 1, &
                          IncludeCorrections, nE, xa, wxa, OpaA)
  CASE (2)
    CALL Integral_2D_Tail(Enu, xT, xmp, xml, xmn, xUp, xUn, xMup, -xMul, xMun, -1.0d0, 2, &
                          IncludeCorrections, nE, xa, wxa, OpaA)
  CASE (3)
    CALL Integral_2D_Tail(Enu, xT, xmp, xml, xmn, xUp, xUn, xMup,  xMul, xMun, -1.0d0, 3, &
                          IncludeCorrections, nE, xa, wxa, OpaA)
  CASE (4)
    CALL Integral_2D_Tail(Enu, xT, xmp, xml, xmn, xUp, xUn, xMup,  xMul, xMun, -1.0d0, 4, &
                          IncludeCorrections, nE, xa, wxa, OpaA)
  CASE DEFAULT
    OpaA = 0.0d0
  END SELECT

END SUBROUTINE Opacity_CC_2D_Tail_Nodes

!=======================================================================
!  Integral_2D_Tail: outer (nucleon) integral for all four reactions.
!  s3 = -1: captures (E4 = Enu + E2 - E3), s3 = +1: (inverse) decay
!  (E4 = Enu + E2 + E3). Occupations as in Integral_2D / Integral_2D_D.
!=======================================================================
SUBROUTINE Integral_2D_Tail(Enu, T, Mass2, Mass3, Mass4, U2, U4, Mu2, Mu3, Mu4, anti, ReactionIndex, &
     IncludeCorrections, nE, xa, wxa, res)

  REAL(DP), INTENT(IN)  :: Enu, T, Mass2, Mass3, Mass4, U2, U4, Mu2, Mu3, Mu4, anti
  INTEGER , INTENT(IN)  :: ReactionIndex, IncludeCorrections(3), nE
  REAL(DP), INTENT(IN)  :: xa(nE), wxa(nE)
  REAL(DP), INTENT(OUT) :: res
#if defined(WEAKLIB_OMP_OL)
  !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
  !$ACC ROUTINE SEQ
#endif

  REAL(DP) :: s3, p2min, p2max, p4max, E2lo, E2s, X
  REAL(DP) :: edge(nEdge), ev(nEvMax), bp(nBpMax), pbp(nBpMax)
  LOGICAL  :: bk(nBpMax)
  REAL(DP) :: pa, pb, xm, wm, P2, E2f, E2, wE2, occ2, sum3
  INTEGER  :: nEv, nb, nFin, i, iSeg, iq
  TYPE(AmpCoeffsType) :: Coeffs

  res = 0.0d0
  IF( ReactionIndex <= 2 ) THEN
    s3 = -1.0d0
    CALL Range_pn( Enu, T, Mass2, Mass3, Mass4, U2, U4, Mu2, p2min, p2max )
    X  = Mu3 + Mu4 - Enu
  ELSE
    s3 = +1.0d0
    ! kinematic p2min: the reaction-3 branch (same kinematics as reaction
    ! 4, without its thermal lower cut on p2); p2max/p4max are not used
    CALL Range_pn_D( Enu, T, Mass2, Mass3, Mass4, U2, U4, Mu2, Mu4, 3, p2min, p2max, p4max )
    X  = Mu4 - Enu - MAX( Mu3, Mass3 )
  END IF
  IF( p2min > p2max ) RETURN   ! kinematically forbidden

  CALL Calc_AmpCoeffs( Mass2, Mass4, U2, U4, anti, IncludeCorrections, Coeffs )

  ! Tail split point: past the threshold, the Fermi edges and the
  ! kinematic events of the outer integrand
  E2lo = SQRT( p2min**2 + Mass2**2 ) + U2
  E2s  = MAX( E2lo + cThr*T, MAX( Mu2, X ) + cEdge*T )
  CALL Outer_Events( s3, Enu, Mass2, Mass3, Mass4, U2, U4, nEv, ev )
  DO i = 1, nEv   ! ascending
    ev(i) = SQRT( ev(i)**2 + Mass2**2 ) + U2
    IF( ev(i) > E2s - cThr*T .AND. ev(i) < E2s + LEv*T ) E2s = MAX( E2s, ev(i) + cThr*T )
  END DO

  edge(1) = Mu2
  edge(2) = X
  edge(3) = E2lo
  CALL Breakpoints( E2lo, E2s, .FALSE., edge, T, nEv, ev, nb, bp, bk )
  DO i = 1, nb
    pbp(i) = SQRT( MAX( (bp(i)-U2)**2 - Mass2**2, 0.0d0 ) )
  END DO
  pbp(1) = p2min

  ! Finite pieces (nE nodes each), then the Gauss-Laguerre tail
  nFin = ( nb - 1 )*nE
  DO i = 1, nFin + nLag
    IF( i <= nFin ) THEN
      iSeg = ( i - 1 )/nE + 1
      iq   = MOD( i - 1, nE ) + 1
      pa = pbp(iSeg)
      pb = pbp(iSeg+1)
      IF( pb <= pa ) CYCLE
      CALL End_Map( xa(iq), bk(iSeg), bk(iSeg+1), xm, wm )
      P2  = pa + xm*( pb - pa )
      E2f = SQRT( P2**2 + Mass2**2 )
      E2  = E2f + U2
      wE2 = wm*wxa(iq)*( pb - pa )*P2/E2f
    ELSE
      E2  = E2s + T*xLag(i-nFin)
      P2  = SQRT( MAX( (E2-U2)**2 - Mass2**2, 0.0d0 ) )
      wE2 = T*wLag(i-nFin)
    END IF
    IF( ReactionIndex == 4 ) THEN
      occ2 = FD(Mu2, E2, T)   ! 1 - f2
    ELSE
      occ2 = FD(E2, Mu2, T)   ! f2
    END IF
    IF( occ2*wE2 <= 0.0d0 ) CYCLE
    CALL Inner_E3( s3, Enu, E2, P2, T, Mass2, Mass3, Mass4, U2, U4, Mu3, Mu4, &
                   ReactionIndex, IncludeCorrections, Coeffs, nE, xa, wxa, sum3 )
    res = res + occ2*sum3*wE2
  END DO

  IF( ReactionIndex <= 3 ) THEN
    ! This is in 1/cm (hbarc in MeV cm)
    res = res*(Gw_MeV*Vud)**2/16.0d0/(pi**5)/hbarc
    ! opacity_unit_factor = 1 unless output_in_km_inverse is set
    res = res*opacity_unit_factor
  ELSE
    res = res*(Gw_MeV*Vud)**2/32.0d0/(pi**7)*Enu**2/hbarc**3/hbar
  END IF

END SUBROUTINE Integral_2D_Tail

!=======================================================================
!  Inner_E3: lepton integral at fixed (E2, P2),
!    sum3 = INT dE3 occ3*occ4*|M|^2/Enu**2
!=======================================================================
SUBROUTINE Inner_E3( s3, Enu, E2, P2, T, Mass2, Mass3, Mass4, U2, U4, Mu3, Mu4, &
                     ReactionIndex, IncludeCorrections, Coeffs, nE, xa, wxa, sum3 )

  REAL(DP), INTENT(IN)  :: s3, Enu, E2, P2, T, Mass2, Mass3, Mass4, U2, U4, Mu3, Mu4
  INTEGER , INTENT(IN)  :: ReactionIndex, IncludeCorrections(3), nE
  TYPE(AmpCoeffsType), INTENT(IN) :: Coeffs
  REAL(DP), INTENT(IN)  :: xa(nE), wxa(nE)
  REAL(DP), INTENT(OUT) :: sum3
#if defined(WEAKLIB_OMP_OL)
  !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
  !$ACC ROUTINE SEQ
#endif

  REAL(DP) :: cap, q0a, q0b, E3a, E3b, E3s, hi
  REAL(DP) :: p3a(nSegMax), p3b(nSegMax), kink(nSegMax), edge(nEdge)
  REAL(DP) :: bp(nBpMax), pbp(nBpMax)
  LOGICAL  :: bk(nBpMax), found, useTail
  REAL(DP) :: pa, pb, xm, wm, P3, E3, E4, wE3, occ, xamp
  INTEGER  :: nSeg, nKink, nb, nFin, nTail, i, j, iSeg, iq

  sum3 = 0.0d0

  ! Kinematic window; for decays without the thermal q0 cap
  IF( s3 < 0.0d0 ) THEN
    cap = Enu - Mass3
  ELSE
    cap = HUGE(1.0d0)
  END IF
  CALL Q0_Window_Envelope( s3, Enu, P2, Mass3, Mass4, U4-E2, cap, q0a, q0b, found )
  IF( .NOT. found ) RETURN
  IF( s3 < 0.0d0 ) THEN
    E3a = Enu - q0b
    E3b = Enu - q0a
  ELSE
    E3a = q0a - Enu
    E3b = q0b - Enu
    IF( q0b >= 0.5d0*HUGE(1.0d0) ) E3b = HUGE(1.0d0)   ! open window
  END IF
  IF( E3b - E3a <= 0.0d0 ) RETURN

  ! Kinks of the angle-integrated |M|^2 (E3_Segments needs a finite end)
  CALL E3_Segments( s3, Enu, E2, P2, Mass3, Mass4, U4, E3a, MIN( E3b, E3a + KinkSearch*T ), &
                    nSeg, p3a, p3b )
  nKink = 0
  kink  = 0.0d0
  DO i = 2, nSeg
    nKink = nKink + 1
    kink(nKink) = SQRT( p3a(i)**2 + Mass3**2 )
  END DO

  ! Fermi edges in E3: lepton (Mu3) and nucleon 4 (E4 = Mu4); window start
  edge(1) = Mu3
  edge(2) = s3*( Mu4 - Enu - E2 )
  edge(3) = E3a

  ! Decays: semi-infinite (or long) windows get the Gauss-Laguerre tail
  useTail = .FALSE.
  hi  = E3b
  E3s = E3b
  IF( s3 > 0.0d0 ) THEN
    E3s = MAX( E3a + cThr*T, MAX( edge(1), edge(2) ) + cEdge*T )
    IF( nKink > 0 ) E3s = MAX( E3s, kink(nKink) + cThr*T )
    IF( E3b - E3s > LMask*T ) THEN
      useTail = .TRUE.
      hi = E3s
    END IF
  END IF

  CALL Breakpoints( E3a, hi, .NOT. useTail, edge, T, nKink, kink, nb, bp, bk )
  DO i = 1, nb
    pbp(i) = SQRT( MAX( bp(i)**2 - Mass3**2, 0.0d0 ) )
  END DO

  ! Finite segments (nE nodes each, in p3), then the tail if any
  nFin  = ( nb - 1 )*nE
  nTail = MERGE( nLag, 0, useTail )
  DO j = 1, nFin + nTail
    IF( j <= nFin ) THEN
      iSeg = ( j - 1 )/nE + 1
      iq   = MOD( j - 1, nE ) + 1
      pa = pbp(iSeg)
      pb = pbp(iSeg+1)
      IF( pb <= pa ) CYCLE
      CALL End_Map( xa(iq), bk(iSeg), bk(iSeg+1), xm, wm )
      P3  = pa + xm*( pb - pa )
      E3  = SQRT( P3**2 + Mass3**2 )
      wE3 = wm*wxa(iq)*( pb - pa )*P3/E3
    ELSE
      E3  = E3s + T*xLag(j-nFin)
      wE3 = T*wLag(j-nFin)
      IF( E3 > E3b ) CYCLE   ! closed but long window: drop nodes past its end
    END IF
    E4 = Enu + E2 + s3*E3
    SELECT CASE( ReactionIndex )
    CASE( 1, 2 )
      occ = FD(Mu3, E3, T)*FD(Mu4, E4, T)   ! (1-f3)(1-f4)
    CASE( 3 )
      occ = FD(E3, Mu3, T)*FD(Mu4, E4, T)   ! f3 (1-f4)
    CASE DEFAULT
      occ = FD(Mu3, E3, T)*FD(E4, Mu4, T)   ! (1-f3) f4, no neutrino blocking
    END SELECT
    IF( occ <= 0.0d0 ) CYCLE
    CALL Calc_Ampsq( Enu, E2, -s3*E3, Mass2, Mass3, Mass4, U2, U4, ReactionIndex, &
                     IncludeCorrections, Coeffs, xamp )
    sum3 = sum3 + occ*xamp*wE3/Enu**2
  END DO

END SUBROUTINE Inner_E3

!=======================================================================
!  Outer_Events: nucleon momenta P2 at which the outer integrand is not
!  analytic. With W = Enu + E2 - U4, P**2 = (Enu +/- P2)**2 (collinear
!  momenta |Enu -/+ P2|) and s = W**2 - P**2:
!    s = (M4+M3)**2, (M4-M3)**2 : a pair of collinear roots (inner kinks)
!                                 appears/disappears (SQRT-type point)
!    s = 0                      : decay window opens/closes at q0 -> inf
!    s - 2*W*M4 + M4**2 - M3**2 = 0    : collinear root at E4 = M4
!    s + 2*s3*W*M3 + M3**2 - M4**2 = 0 : collinear root at E3 = M3
!    P2 = Enu                   : |Enu - P2| changes branch
!  Each condition reads c1*SQRT(P2**2+M2**2) = c0 + c2*P2 and is solved
!  exactly (roots of the squared equation with the wrong sign dropped).
!  Returned sorted ascending, as momenta.
!=======================================================================
SUBROUTINE Outer_Events( s3, Enu, Mass2, Mass3, Mass4, U2, U4, nEv, ev )

  REAL(DP), INTENT(IN)  :: s3, Enu, Mass2, Mass3, Mass4, U2, U4
  INTEGER , INTENT(OUT) :: nEv
  REAL(DP), INTENT(OUT) :: ev(nEvMax)
#if defined(WEAKLIB_OMP_OL)
  !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
  !$ACC ROUTINE SEQ
#endif

  REAL(DP) :: a, al, be, c0, c1, c2, qa, qb, qc, disc, r(2), sg, tmp
  INTEGER  :: iC, iS, ir, i, j

  nEv   = 1
  ev    = 0.0d0
  ev(1) = Enu
  a = Enu + U2 - U4          ! W = SQRT(P2**2+M2**2) + a
  DO iC = 1, 5
    SELECT CASE( iC )
    CASE( 1 )
      al = 0.0d0
      be = -(Mass4 + Mass3)**2
    CASE( 2 )
      al = 0.0d0
      be = -(Mass4 - Mass3)**2
    CASE( 3 )
      al = 0.0d0
      be = 0.0d0
    CASE( 4 )
      al = -2.0d0*Mass4
      be = Mass4**2 - Mass3**2
    CASE DEFAULT
      al = 2.0d0*s3*Mass3
      be = Mass3**2 - Mass4**2
    END SELECT
    ! W**2 + al*W + be = (Enu + sg*P2)**2  <=>
    ! (2a+al)*SQRT(P2**2+M2**2) = -(M2**2 + a**2 + al*a + be - Enu**2) + 2*sg*Enu*P2
    c1 = 2.0d0*a + al
    c0 = -( Mass2**2 + a**2 + al*a + be - Enu**2 )
    DO iS = 1, 2
      sg = MERGE( 1.0d0, -1.0d0, iS == 1 )
      c2 = 2.0d0*sg*Enu
      qa = c1**2 - c2**2
      qb = -2.0d0*c0*c2
      qc = c1**2*Mass2**2 - c0**2
      IF( ABS(qa) <= 1.0d-14*( c1**2 + c2**2 ) ) CYCLE
      disc = qb**2 - 4.0d0*qa*qc
      IF( disc < 0.0d0 ) CYCLE
      ! both roots, cancellation-free
      tmp  = -0.5d0*( qb + SIGN( SQRT(disc), qb ) )
      r(1) = tmp/qa
      r(2) = 0.0d0
      IF( tmp /= 0.0d0 ) r(2) = qc/tmp
      DO ir = 1, 2
        IF( r(ir) <= 0.0d0 ) CYCLE
        IF( c1*( c0 + c2*r(ir) ) < 0.0d0 ) CYCLE   ! root of the squared equation only
        nEv = nEv + 1
        ev(nEv) = r(ir)
      END DO
    END DO
  END DO

  DO i = 2, nEv   ! insertion sort
    tmp = ev(i)
    j   = i - 1
    DO WHILE( j >= 1 )
      IF( ev(j) <= tmp ) EXIT
      ev(j+1) = ev(j)
      j = j - 1
    END DO
    ev(j+1) = tmp
  END DO

END SUBROUTINE Outer_Events

!=======================================================================
!  End_Map: node u in [0,1] -> x in [0,1] with quadratic clustering at
!  the kinematic ends (x ~ u**2 there), which makes SQRT-type endpoint
!  behaviour analytic in u; plain Gauss-Legendre otherwise.
!=======================================================================
SUBROUTINE End_Map( u, mapL, mapR, xm, wm )

  REAL(DP), INTENT(IN)  :: u
  LOGICAL , INTENT(IN)  :: mapL, mapR
  REAL(DP), INTENT(OUT) :: xm, wm
#if defined(WEAKLIB_OMP_OL)
  !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
  !$ACC ROUTINE SEQ
#endif

  IF( mapL .AND. mapR ) THEN
    xm = u*u*( 3.0d0 - 2.0d0*u )
    wm = 6.0d0*u*( 1.0d0 - u )
  ELSE IF( mapL ) THEN
    xm = u*u
    wm = 2.0d0*u
  ELSE IF( mapR ) THEN
    xm = 1.0d0 - ( 1.0d0 - u )**2
    wm = 2.0d0*( 1.0d0 - u )
  ELSE
    xm = u
    wm = 1.0d0
  END IF

END SUBROUTINE End_Map

!=======================================================================
!  Breakpoints: sorted breakpoints in [lo, hi]: the ends, nK kinematic
!  points (kinks/events) and the graded points e +/- Grade*T around each
!  Fermi edge e. bk = .TRUE. marks kinematic points (lo, the kinks/events
!  and hi if hiK), which get the end clustering of End_Map.
!=======================================================================
SUBROUTINE Breakpoints( lo, hi, hiK, edge, T, nK, kink, nb, bp, bk )

  REAL(DP), INTENT(IN)  :: lo, hi, edge(nEdge), T
  LOGICAL , INTENT(IN)  :: hiK
  INTEGER , INTENT(IN)  :: nK
  REAL(DP), INTENT(IN)  :: kink(*)
  INTEGER , INTENT(OUT) :: nb
  REAL(DP), INTENT(OUT) :: bp(nBpMax)
  LOGICAL , INTENT(OUT) :: bk(nBpMax)
#if defined(WEAKLIB_OMP_OL)
  !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
  !$ACC ROUTINE SEQ
#endif

  REAL(DP) :: v, tol, tmp
  LOGICAL  :: tk
  INTEGER  :: i, j, ie, ig, is

  tol = 1.0d-9*( hi - lo ) + 1.0d-12*( ABS(lo) + ABS(hi) )
  bp = 0.0d0
  bk = .FALSE.
  nb = 1
  bp(1) = lo
  bk(1) = .TRUE.
  DO i = 1, nK
    IF( kink(i) > lo + tol .AND. kink(i) < hi - tol ) THEN
      nb = nb + 1
      bp(nb) = kink(i)
      bk(nb) = .TRUE.
    END IF
  END DO
  DO ie = 1, nEdge
    DO ig = 1, nGrade
      DO is = -1, 1, 2
        v = edge(ie) + is*Grade(ig)*T
        IF( v > lo + tol .AND. v < hi - tol ) THEN
          nb = nb + 1
          bp(nb) = v
          bk(nb) = .FALSE.
        END IF
      END DO
    END DO
  END DO
  nb = nb + 1
  bp(nb) = hi
  bk(nb) = hiK

  DO i = 3, nb - 1      ! insertion sort of the interior points
    tmp = bp(i)
    tk  = bk(i)
    j   = i - 1
    DO WHILE( j >= 2 )
      IF( bp(j) <= tmp ) EXIT
      bp(j+1) = bp(j)
      bk(j+1) = bk(j)
      j = j - 1
    END DO
    bp(j+1) = tmp
    bk(j+1) = tk
  END DO

END SUBROUTINE Breakpoints

END MODULE wlSemiLeptonicIntegration2D_Tail
