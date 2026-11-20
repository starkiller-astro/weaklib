PROGRAM wlTestKernelHalving

  ! Tests whether the general scattering kernels can be computed on the
  ! lower energy triangle only (E3 <= E1, i.e. iE3 <= iE1), for ALL processes
  ! 1..36, with the upper triangle rebuilt from a partner process:
  !
  !   NES/NMS  1..24 : Phout_p(E1,E3) = Phout_p(E3,E1) * DB_p(E3,E1)   (partner = self)
  !   LFE/LFC 25..32 : Phout_r(E1,E3) = Phout_f(E3,E1) * DB_f(E3,E1)   (f,r) = (25,26) (27,28)
  !                                                                         (29,30) (31,32)
  !   IMD     33..36 : Phout_b(E1,E3) = Phout_a(E3,E1)                 (a,b) = (33,34) (35,36)
  !
  ! with DB = DetailedBalanceFactor (Phin = Phout * DB). These follow from
  ! Phin_f(E1,E3) = Phout_r(E3,E1) for inverse pairs, and the nu_1 <-> nu_3
  ! symmetry of IMD (see wlTestInverseReactionSymmetry). For 1..24 this is
  ! what thornado's ComputeNeutrinoOpacityRates_MuonNLS already does; thornado
  ! currently computes the full nE x nE square for every process as soon as
  ! any of 25..36 is selected.
  !
  ! The test:
  !   1. computes the full square kernels with CalculateAllPhout (the scalar
  !      routine thornado calls) and times it,
  !   2. computes the lower triangle only, rebuilds the upper triangle from the
  !      partner rows, and times it,
  !   3. reports the max relative difference per process, l = 0 and l = 1.
  !
  ! The energy grid mimics thornado: nE zones, geometric zoom, nNodes Gauss
  ! nodes per zone. Thermodynamic points are (T, Mu_e, Mu_mu) in MeV.
  !
  ! Usage:  ./wlTestKernelHalving_$(MACHINE) [nTheta]      (default nTheta = 24)
  !
  ! Serial, no tables read, a few seconds to a minute with the default setup.

  USE wlKindModule, ONLY: dp
  USE wlGeneralLeptonScatteringModuleThornadoInterface, ONLY: &
    InitGeneralScatteringKernels, &
    FinalizeGeneralScatteringKernels, &
    CalculateAllPhout, &
    DetailedBalanceFactor, &
    ProcessRowIndex

  IMPLICIT NONE

  ! --- thornado-like energy grid (setup: nE=16, nNodes=2, eR=250, zoomE=1.2462)
  INTEGER , PARAMETER :: nZones = 16, nNodes = 2, nE = nZones * nNodes
  REAL(DP), PARAMETER :: eL = 0.0d0, eR = 250.0d0, ZoomE = 1.2462d0

  INTEGER , PARAMETER :: nP = 36, nL = 2, nX = 6
  REAL(DP), PARAMETER :: RelFloor = 1.0d-10  ! ignore entries below this * max|ref|
  REAL(DP), PARAMETER :: AbsFloor = 1.0d-10  ! [cm^-1] ignore entries with |Phi|*E3^3 below this
  REAL(DP), PARAMETER :: L1Reg    = 1.0d-6   ! l=1: error relative to max(|ref|, L1Reg*max|ref|)
  REAL(DP), PARAMETER :: RelTol   = 1.0d-8

  REAL(DP) :: E(nE), T(nX), Me(nX), Mm(nX)
  REAL(DP), ALLOCATABLE :: Full(:,:,:,:,:), Half(:,:,:,:,:)   ! (nE,nE,nP,nL,nX)
  REAL(DP) :: Ph(nP,nL)
  REAL(DP) :: tFull, tHalf, Err(nL), ErrMax(nL)
  INTEGER  :: nTheta, iX, iP, iL, i, j, iRow, nFail
  INTEGER  :: ProcList(nP), RowOf(nP)
  INTEGER(8) :: c0, c1, rate
  CHARACTER(32) :: arg

  ! ------------------------------------------------------------ setup
  nTheta = 24
  IF ( COMMAND_ARGUMENT_COUNT() >= 1 ) THEN
    CALL GET_COMMAND_ARGUMENT( 1, arg )
    READ( arg, * ) nTheta
  ENDIF

  CALL BuildEnergyGrid( E )

  ! (T, Mu_e, Mu_mu) [MeV]: from mildly to strongly degenerate, plus a hot
  ! point with Mu_mu well below m_mu and a cold, very degenerate one.
  T  = (/  5.0d0, 10.0d0, 15.2d0,  30.0d0, 40.0d0,   8.0d0 /)
  Me = (/ 20.0d0, 60.0d0, 44.4d0, 120.0d0, 10.0d0, 180.0d0 /)
  Mm = (/105.0d0, 95.0d0, 37.3d0, 150.0d0, 20.0d0, 175.0d0 /)

  ProcList = (/ (iP, iP = 1, nP) /)
  CALL InitGeneralScatteringKernels( E, E, nE, nE, nTheta, &
                                     iProcessList_Option = ProcList )
  DO iP = 1, nP
    RowOf(iP) = ProcessRowIndex( iP )
  ENDDO

  ALLOCATE( Full(nE,nE,nP,nL,nX), Half(nE,nE,nP,nL,nX) )
  Full = 0.0d0
  Half = 0.0d0

  WRITE(*,'(A)')        ' wlTestKernelHalving'
  WRITE(*,'(A,I4,A,I3,A,I3,A,I3)') '   nE =', nE, '   nTheta =', nTheta, &
                                   '   nProcess =', nP, '   nX =', nX
  WRITE(*,'(A,2F9.3,A)') '   E range = ', E(1), E(nE), ' MeV'

  CALL SYSTEM_CLOCK( COUNT_RATE = rate )

  ! ------------------------------------------------------------ 1. full square
  CALL SYSTEM_CLOCK( c0 )
  DO iX = 1, nX
    DO j = 1, nE          ! iE3
      DO i = 1, nE        ! iE1
        CALL CalculateAllPhout( i, j, E(i), E(j), T(iX), Me(iX), Mm(iX), Ph, nL )
        DO iP = 1, nP
          Full(i,j,iP,:,iX) = Ph(RowOf(iP),:)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  CALL SYSTEM_CLOCK( c1 )
  tFull = REAL( c1 - c0, DP ) / REAL( rate, DP )

  ! ------------------------------------------------------------ 2. lower triangle + rebuild
  CALL SYSTEM_CLOCK( c0 )
  DO iX = 1, nX
    DO j = 1, nE          ! iE3
      DO i = j, nE        ! iE1 >= iE3
        CALL CalculateAllPhout( i, j, E(i), E(j), T(iX), Me(iX), Mm(iX), Ph, nL )
        DO iP = 1, nP
          Half(i,j,iP,:,iX) = Ph(RowOf(iP),:)
        ENDDO
      ENDDO
    ENDDO
    ! upper triangle (i < j) from the partner's lower-triangle entry (j,i)
    DO iP = 1, nP
      DO j = 1, nE
        DO i = 1, j - 1
          Half(i,j,iP,:,iX) = Half(j,i,Partner(iP),:,iX) &
                              * RebuildFactor( iP, E(j), E(i), T(iX), Me(iX), Mm(iX) )
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  CALL SYSTEM_CLOCK( c1 )
  tHalf = REAL( c1 - c0, DP ) / REAL( rate, DP )

  ! ------------------------------------------------------------ 3. compare
  WRITE(*,*)
  WRITE(*,'(A)') '   process  partner    max rel err  l=0        l=1'
  nFail  = 0
  ErrMax = 0.0d0
  DO iP = 1, nP
    Err = 0.0d0
    DO iX = 1, nX
      DO iL = 1, nL
        Err(iL) = MAX( Err(iL), MaxRelErr( Half(:,:,iP,iL,iX), Full(:,:,iP,iL,iX), iL ) )
      ENDDO
    ENDDO
    ErrMax = MAX( ErrMax, Err )
    WRITE(*,'(A,I6,I9,A,2ES11.3)', ADVANCE='NO') '   ', iP, Partner(iP), '      ', Err
    IF ( ANY( Err > RelTol ) ) THEN
      nFail = nFail + 1
      WRITE(*,'(A)') '   <-- FAIL'
    ELSE
      WRITE(*,*)
    ENDIF
  ENDDO

  WRITE(*,*)
  WRITE(*,'(A,2ES11.3)') '   overall max rel err (l=0, l=1): ', ErrMax
  WRITE(*,'(A,F10.3,A)') '   full square     : ', tFull, ' s'
  WRITE(*,'(A,F10.3,A)') '   lower + rebuild : ', tHalf, ' s'
  IF ( tHalf > 0.0d0 ) &
    WRITE(*,'(A,F7.3,A,F7.3,A)') '   speedup         : ', tFull / tHalf, &
      '   (ideal 2n/(n+1) = ', 2.0d0 * nE / ( nE + 1.0d0 ), ')'

  CALL FinalizeGeneralScatteringKernels()
  DEALLOCATE( Full, Half )

  WRITE(*,*)
  IF ( nFail == 0 ) THEN
    WRITE(*,'(A)') ' wlTestKernelHalving: PASSED'
  ELSE
    WRITE(*,'(A,I3,A)') ' wlTestKernelHalving: FAILED (', nFail, ' processes)'
  ENDIF

CONTAINS

  ! Process whose (E3,E1) entry rebuilds this process's (E1,E3) entry.
  PURE INTEGER FUNCTION Partner( iProc )
    INTEGER, INTENT(IN) :: iProc
    SELECT CASE ( iProc )
    CASE ( 1:24 )
      Partner = iProc                       ! elastic in flavour: itself
    CASE ( 25, 27, 29, 31, 33, 35 )
      Partner = iProc + 1
    CASE ( 26, 28, 30, 32, 34, 36 )
      Partner = iProc - 1
    CASE DEFAULT
      Partner = -1
    END SELECT
  END FUNCTION Partner

  ! Factor multiplying Phout_partner(Ea,Eb) to give Phout_iProc(Eb,Ea).
  ! 2->2 (1..32): Phout_x(Eb,Ea) = Phin_partner(Ea,Eb) = Phout_partner(Ea,Eb) * DB_partner(Ea,Eb)
  ! IMD (33..36): pure transpose.
  REAL(DP) FUNCTION RebuildFactor( iProc, Ea, Eb, Tk, Mu_e, Mu_mu )
    INTEGER , INTENT(IN) :: iProc
    REAL(DP), INTENT(IN) :: Ea, Eb, Tk, Mu_e, Mu_mu
    IF ( iProc >= 33 ) THEN
      RebuildFactor = 1.0d0
    ELSE
      RebuildFactor = DetailedBalanceFactor( RowOf(Partner(iProc)), Ea, Eb, Tk, Mu_e, Mu_mu )
    ENDIF
  END FUNCTION RebuildFactor

  ! Max relative error of A against reference B over the whole square.
  ! l=0 (iL=1) is pointwise relative; l=1 changes sign, so its error is taken
  ! relative to max(|B|, L1Reg*max|B|). Entries negligible in absolute or
  ! relative terms are skipped (same floors as wlTestInverseReactionSymmetry).
  REAL(DP) FUNCTION MaxRelErr( A, B, iL )
    REAL(DP), INTENT(IN) :: A(nE,nE), B(nE,nE)
    INTEGER , INTENT(IN) :: iL
    REAL(DP) :: Scale, Ref, Den
    INTEGER  :: i, j

    MaxRelErr = 0.0d0
    Scale = MAXVAL( ABS( B ) )
    IF ( Scale <= 0.0d0 ) THEN
      ! reference identically zero: rebuilt kernel must be zero too
      IF ( MAXVAL( ABS( A ) ) > 0.0d0 ) MaxRelErr = HUGE( 1.0d0 )
      RETURN
    ENDIF
    DO j = 1, nE
      DO i = 1, nE
        Ref = B(i,j)
        IF ( ABS(Ref) <= RelFloor * Scale )   CYCLE
        IF ( ABS(Ref) * E(j)**3 <= AbsFloor ) CYCLE
        IF ( iL == 1 ) THEN
          Den = ABS(Ref)
        ELSE
          Den = MAX( ABS(Ref), L1Reg * Scale )
        ENDIF
        MaxRelErr = MAX( MaxRelErr, ABS( A(i,j) - Ref ) / Den )
      ENDDO
    ENDDO
  END FUNCTION MaxRelErr

  ! Node energies of a thornado-style geometric grid on [eL, eR] with
  ! nZones zones (width ratio ZoomE) and nNodes Gauss-Legendre nodes per zone.
  SUBROUTINE BuildEnergyGrid( En )
    REAL(DP), INTENT(OUT) :: En(nE)
    REAL(DP) :: w1, left, width, xi(nNodes)
    INTEGER  :: k, n

    SELECT CASE ( nNodes )
    CASE ( 1 ); xi(1) = 0.0d0
    CASE ( 2 ); xi(1) = -0.5d0 / SQRT(3.0d0)
                xi(nNodes) = 0.5d0 / SQRT(3.0d0)
    END SELECT

    w1   = ( eR - eL ) * ( ZoomE - 1.0d0 ) / ( ZoomE**nZones - 1.0d0 )
    left = eL
    DO k = 1, nZones
      width = w1 * ZoomE**(k-1)
      DO n = 1, nNodes
        En( (k-1)*nNodes + n ) = left + width * ( 0.5d0 + xi(n) )
      ENDDO
      left = left + width
    ENDDO
  END SUBROUTINE BuildEnergyGrid

END PROGRAM wlTestKernelHalving
