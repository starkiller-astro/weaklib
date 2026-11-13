PROGRAM wlTestInverseAntiMuonDecay

  ! Checks for inverse anti-muon decay (processes 35, 36):
  !   1. CP identity: Phout_35(T,-Mu_e,-Mu_mu) == Phout_33(T,Mu_e,Mu_mu), and
  !      36 vs 34, bit-for-bit, for CalculateAllPhout and CalculateAllPhout_3D.
  !   2. Interface Phout vs the reference CalculatePhoutPhin, processes 33..36.
  !   3. Phout * DetailedBalanceFactor vs the reference Phin, processes 33..36
  !      and one scattering process.

  USE wlKindModule, ONLY: dp
  USE wlEosConstantsModule, ONLY: ElectronMass => me, MuonMass => mmu
  USE wlGeneralLeptonScatteringModule, ONLY: &
    CalculatePhoutPhin
  USE wlGeneralLeptonScatteringModuleThornadoInterface, ONLY: &
    InitGeneralScatteringKernels, &
    FinalizeGeneralScatteringKernels, &
    CalculateAllPhout, &
    CalculateAllPhout_3D, &
    DetailedBalanceFactor, &
    ProcessRowIndex

  IMPLICIT NONE

  INTEGER , PARAMETER :: nE = 16, nTh = 24, nL = 2, nX = 4, nP = 5
  INTEGER , PARAMETER :: ProcList(nP) = (/ 33, 34, 35, 36, 1 /)
  REAL(DP), PARAMETER :: RelFloor = 1.0d-10   ! ignore entries below this * max|ref|
  REAL(DP), PARAMETER :: RelTol   = 1.0d-6

  REAL(DP) :: E(nE), T(nX), Me(nX), Mm(nX)
  REAL(DP) :: Ph_p(nP,nL), Ph_m(nP,nL)
  REAL(DP) :: Ph3_p(nE,nE,nP,nL,nX), Ph3_m(nE,nE,nP,nL,nX)
  REAL(DP) :: PhInt(nE,nE,nP,nL,nX), PhinInt(nE,nE,nP,nL,nX)
  REAL(DP) :: PhRef(nE,nE,nP,nL,nX), PhinRef(nE,nE,nP,nL,nX)
  REAL(DP) :: Phout_r(nL), Phin_r(nL), xMu2, xMu4, m2, m4, sgn
  INTEGER  :: i, j, iX, iP, iProc, r33, r34, r35, r36
  LOGICAL  :: AllOK

  DO i = 1, nE
    E(i) = 2.0d0 * 1.3d0**(i-1)
  END DO
  ! (T, Mu_e, Mu_mu) [MeV]; the last point is hot with mu_mu well below m_mu,
  ! where anti-muons are most abundant.
  T  = (/  5.0d0, 15.2d0,  30.0d0, 40.0d0 /)
  Me = (/ 20.0d0, 44.4d0, 120.0d0, 10.0d0 /)
  Mm = (/105.0d0, 37.3d0, 150.0d0, 20.0d0 /)

  CALL InitGeneralScatteringKernels( E, E, nE, nE, nTh, iProcessList_Option = ProcList )

  r33 = ProcessRowIndex(33); r34 = ProcessRowIndex(34)
  r35 = ProcessRowIndex(35); r36 = ProcessRowIndex(36)

  AllOK = .TRUE.

  ! ---------------------------------------------------------------- 1. CP
  CALL CalculateAllPhout_3D( nX, nE, nE, nP, nL, T,  Me,  Mm, Ph3_p )
  CALL CalculateAllPhout_3D( nX, nE, nE, nP, nL, T, -Me, -Mm, Ph3_m )

  DO iX = 1, nX
    DO j = 1, nE
      DO i = 1, nE
        CALL CalculateAllPhout( i, j, E(i), E(j), T(iX),  Me(iX),  Mm(iX), Ph_p, nL )
        CALL CalculateAllPhout( i, j, E(i), E(j), T(iX), -Me(iX), -Mm(iX), Ph_m, nL )
        IF ( ANY( Ph_m(r35,:) /= Ph_p(r33,:) ) .OR. &
             ANY( Ph_m(r36,:) /= Ph_p(r34,:) ) ) THEN
          AllOK = .FALSE.
          WRITE(*,*) 'CP FAIL (point) at iX,iE1,iE3 =', iX, i, j
        ENDIF
        ! The batched kernel must agree with the point kernel exactly
        IF ( ANY( Ph3_p(i,j,:,:,iX) /= Ph_p ) ) THEN
          AllOK = .FALSE.
          WRITE(*,*) 'POINT vs 3D FAIL at iX,iE1,iE3 =', iX, i, j
        ENDIF
        PhInt(i,j,:,:,iX) = Ph_p
        DO iP = 1, nP
          PhinInt(i,j,iP,:,iX) = Ph_p(iP,:) &
            * DetailedBalanceFactor( iP, E(i), E(j), T(iX), Me(iX), Mm(iX) )
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  IF ( ANY( Ph3_m(:,:,r35,:,:) /= Ph3_p(:,:,r33,:,:) ) .OR. &
       ANY( Ph3_m(:,:,r36,:,:) /= Ph3_p(:,:,r34,:,:) ) ) THEN
    AllOK = .FALSE.
    WRITE(*,*) 'CP FAIL (3D)'
  ENDIF
  WRITE(*,*) 'CP identity 35<->33, 36<->34 (point and 3D) checked'

  ! ----------------------------------------------------- 2./3. reference
  DO iP = 1, nP
    iProc = ProcList(iP)
    SELECT CASE (iProc)
    CASE (33, 34); sgn = +1.0d0; m2 = ElectronMass; m4 = MuonMass
    CASE (35, 36); sgn = -1.0d0; m2 = ElectronMass; m4 = MuonMass
    CASE (1);      sgn = +1.0d0; m2 = ElectronMass; m4 = ElectronMass
    END SELECT
    DO iX = 1, nX
      xMu2 = sgn * Me(iX)
      IF ( iProc == 1 ) THEN
        xMu4 = Me(iX)
      ELSE
        xMu4 = sgn * Mm(iX)
      ENDIF
      DO j = 1, nE
        DO i = 1, nE
          CALL CalculatePhoutPhin( E(i), E(j), T(iX), xMu2, xMu4, m2, m4, &
                                   iProc, Phout_r, Phin_r, nL, nTh )
          PhRef  (i,j,iP,:,iX) = Phout_r
          PhinRef(i,j,iP,:,iX) = Phin_r
        ENDDO
      ENDDO
    ENDDO
    CALL Report( 'Phout', iProc, PhInt  (:,:,iP,:,:), PhRef  (:,:,iP,:,:) )
    CALL Report( 'Phin ', iProc, PhinInt(:,:,iP,:,:), PhinRef(:,:,iP,:,:) )
  ENDDO

  CALL FinalizeGeneralScatteringKernels()

  IF ( AllOK ) THEN
    WRITE(*,*) 'wlTestInverseAntiMuonDecay: PASSED'
  ELSE
    WRITE(*,*) 'wlTestInverseAntiMuonDecay: FAILED'
  ENDIF

CONTAINS

  SUBROUTINE Report( Label, iProc, A, B )
    CHARACTER(*), INTENT(IN) :: Label
    INTEGER     , INTENT(IN) :: iProc
    REAL(DP)    , INTENT(IN) :: A(nE,nE,nL,nX), B(nE,nE,nL,nX)
    REAL(DP) :: ErrMax, Scale
    INTEGER  :: iL, k

    ErrMax = 0.0d0
    DO k = 1, nX
      DO iL = 1, nL
        Scale = MAXVAL( ABS( B(:,:,iL,k) ) )
        IF ( Scale <= 0.0d0 ) CYCLE
        ErrMax = MAX( ErrMax, MAXVAL( ABS( A(:,:,iL,k) - B(:,:,iL,k) ), &
                      MASK = ABS( B(:,:,iL,k) ) > RelFloor * Scale ) / Scale )
      ENDDO
    ENDDO
    WRITE(*,'(A,A,I3,A,ES10.3)') '  ', Label, iProc, '  max rel err vs reference:', ErrMax
    IF ( ErrMax > RelTol ) AllOK = .FALSE.
  END SUBROUTINE Report

END PROGRAM wlTestInverseAntiMuonDecay
