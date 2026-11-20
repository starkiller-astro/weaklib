PROGRAM wlTestAngularQuadrature

  ! Convergence of the angular (cos theta) quadrature used by the general
  ! scattering kernels, for ALL processes 1..36, at realistic thermodynamic
  ! conditions.
  !
  ! The kernels are computed with nThetaRef nodes (default 128) as the
  ! reference, and compared with nTheta = 12, 16 and 24 (default). A run with
  ! 2*nThetaRef checks that the reference itself is converged. The node layout
  ! (per (E1, E3, mass pair), mapped in 1 - mu, IMD range cut at its kinematic
  ! edge) is described in the header of
  ! wlGeneralLeptonScatteringModuleThornadoInterface.
  !
  ! Thermodynamic points (default): ThermoConditions.dat, rows
  !   T Rho Ye Ym Mue Mum Mun Mup ...     [MeV, g/cm^3; Mu include rest mass]
  ! keeping Rho >= RhoMin (1e10 g/cm^3). Alternatively a synthetic
  ! (T, Mu_e, Mu_mu) grid ('grid' as first argument).
  !
  ! Compared quantity: the energy-integrated "opacity" moment of process p
  !     kappa_p,l(E1) = Sum_j Phout_p,l(E1,E3_j) E3_j^2 dE3_j          [cm^-1]
  ! (no final-state blocking / incoming-occupation factor: an upper bound on
  ! what the kernel feeds into the transport). Each process is attributed to
  ! its INCOMING neutrino species s (nu_e, nu_bar_e, nu_mu, nu_bar_mu, nu_tau,
  ! nu_bar_tau, parsed from the reaction string), and
  !     kappa_tot,s(E1) = Sum_{p in s} kappa_p,0(E1)       (reference nTheta)
  ! is that species' total opacity FROM THESE KERNELS (absorption on nucleons
  ! etc. is not included, so this underestimates the true total and the
  ! measures below are conservative).
  !
  ! Three error measures, l=1 always normalized by 3x the l=0 scale (Phi_1
  ! carries Bruenn's 3/2 vs 1/2, so kappa_1 enters the transport as
  ! kappa_0 - kappa_1/3; kappa_1 changes sign, so its own value is no scale):
  !   own   : max_E1 |dkappa_p| / kappa_p,0                 (kappa_p,0   >= floor)
  !   share : max_E1 |dkappa_p| / kappa_tot,s               (kappa_tot,s >= floor)
  !           -> the error as a fraction of everything this neutrino scatters
  !              on at that energy; a 20% error in a negligible channel is small
  !   spec  : Sum_E1 w |dkappa_p| / Sum_E1 w kappa_tot,s,  w = E1^2 dE1 f_FD(E1)
  !           f_FD at the local T and the equilibrium neutrino chemical potential
  !           mu_nu_e = Mu_e - (Mu_n - Mu_p), mu_nu_mu = Mu_mu - (Mu_n - Mu_p),
  !           antineutrinos -mu, tau 0 (all 0 in 'grid' mode)
  !           -> discounts energies the local neutrino spectrum does not populate
  ! floor = KappaFloor = 1e-10 cm^-1.
  !
  ! NES (processes 1..12) is NOT tested: it can be tabulated in (E1,E3,T,eta)
  ! like weaklib's existing NES tables, so its on-the-fly quadrature accuracy
  ! does not matter. NES is still computed, because it is part of every
  ! species' total opacity in the transport and so belongs in kappa_tot,s.
  !
  ! PASS when the SHARE error of every tested process (13..36) is within
  ! OpacTol (1%) for the LARGEST tested nTheta. A warning is printed if the reference check exceeds
  ! RefTol.
  !
  ! Output:
  !   stdout : per process, share error for each nTheta (max over l), own and
  !            spec errors for the largest nTheta, and where the worst share
  !            error occurs; then a per-group summary of all three measures
  !   AngularQuadrature_ref<nRef>.dat : per thermodynamic point, for each
  !            tested nTheta, worst share then worst spec error per tested
  !            process group (NMS LFE LFC IMD), max over l; for plotting
  !
  ! Usage:  ./wlTestAngularQuadrature_$(MACHINE) [source] [nThetaRef [nTheta1 nTheta2 ...]]
  !         source = a ThermoConditions-format file (default ThermoConditions.dat)
  !                  or 'grid'
  !         default: ThermoConditions.dat 128 12 16 24
  !
  ! Thermodynamic points are OpenMP-parallel.

  USE wlKindModule, ONLY: dp
  USE wlGeneralLeptonScatteringModule, ONLY: &
    ReactionStringFromProcessIndex
  USE wlGeneralLeptonScatteringModuleThornadoInterface, ONLY: &
    InitGeneralScatteringKernels, &
    FinalizeGeneralScatteringKernels, &
    CalculateAllPhout, &
    ProcessRowIndex

  IMPLICIT NONE

  ! --- thornado-like energy grid (setup: nE=16, nNodes=2, eR=250, zoomE=1.2462)
  INTEGER , PARAMETER :: nZones = 16, nNodes = 2, nE = nZones * nNodes
  REAL(DP), PARAMETER :: eL = 0.0d0, eR = 250.0d0, ZoomE = 1.2462d0

  ! --- synthetic thermodynamic grid [MeV] ('grid' mode)
  INTEGER , PARAMETER :: nT = 8, nMe = 7, nMm = 7
  REAL(DP), PARAMETER :: T_List (nT)  = (/ 1.0d0, 2.5d0, 5.0d0, 10.0d0, 20.0d0, &
                                           30.0d0, 50.0d0, 80.0d0 /)
  REAL(DP), PARAMETER :: Me_List(nMe) = (/ 0.0d0, 5.0d0, 20.0d0, 50.0d0, 100.0d0, &
                                           200.0d0, 300.0d0 /)
  REAL(DP), PARAMETER :: Mm_List(nMm) = (/ 0.0d0, 50.0d0, 95.0d0, 105.0d0, 120.0d0, &
                                           170.0d0, 250.0d0 /)

  INTEGER , PARAMETER :: nP = 36, nL = 2, nGroup = 5, nSpec = 6, nTestMax = 8
  INTEGER , PARAMETER :: iOwn = 1, iShare = 2, iSpec = 3, nMeas = 3
  INTEGER , PARAMETER :: iGroupNES = 1          ! not tested (tabulated instead)
  REAL(DP), PARAMETER :: RhoMin     = 1.0d10   ! [g/cm^3]
  REAL(DP), PARAMETER :: KappaFloor = 1.0d-10  ! [cm^-1] ignore below this
  REAL(DP), PARAMETER :: OpacTol    = 1.0d-2
  REAL(DP), PARAMETER :: RefTol     = 1.0d-3   ! reference vs 2x reference

  CHARACTER(3), PARAMETER :: GroupName(nGroup) = (/ 'NES', 'NMS', 'LFE', 'LFC', 'IMD' /)
  CHARACTER(10), PARAMETER :: SpecName(nSpec) = &
    (/ 'nu_e      ', 'nu_bar_e  ', 'nu_mu     ', 'nu_bar_mu ', 'nu_tau    ', 'nu_bar_tau' /)
  CHARACTER(5), PARAMETER :: MeasName(nMeas) = (/ 'own  ', 'share', 'spec ' /)

  INTEGER :: nX
  REAL(DP) :: E(nE), dE(nE)
  REAL(DP), ALLOCATABLE :: T(:), Rho(:), Me(:), Mm(:), MuNu(:,:)   ! MuNu(nSpec,nX)
  REAL(DP), ALLOCATABLE :: KapRef(:,:,:,:)        ! (nE,nP,nL,nX)  reference opacities
  REAL(DP), ALLOCATABLE :: KapTot(:,:,:)          ! (nE,nSpec,nX)  reference species total, l=0
  REAL(DP), ALLOCATABLE :: Err(:,:,:,:,:)         ! (nP,nL,nX,nRun,nMeas)
  INTEGER , ALLOCATABLE :: iEErr(:,:,:,:)         ! (nP,nL,nX,nRun) E1 index of share error
  REAL(DP) :: GrpErr(nGroup,nMeas), Worst, tRun, v
  INTEGER  :: nThetaRef, nTest, nTheta(nTestMax+1), iTest, iX, iP, iL, iG, iM, iS
  INTEGER  :: nFail, iWorst, lWorst, iUnit, nArg, iArg, iArg0, iLast, ios, iTmp
  INTEGER  :: ProcList(nP), SpecOf(nP)
  LOGICAL  :: GridMode, Tested(nP)
  CHARACTER(256) :: arg, Source
  CHARACTER(64)  :: FileName

  ! ------------------------------------------------------------ arguments
  Source    = 'ThermoConditions.dat'
  nThetaRef = 128
  nTest     = 3
  nTheta(1:3) = (/ 12, 16, 24 /)
  nArg  = COMMAND_ARGUMENT_COUNT()
  iArg0 = 1
  IF ( nArg >= 1 ) THEN
    CALL GET_COMMAND_ARGUMENT( 1, arg )
    READ( arg, *, IOSTAT = ios ) iTmp
    IF ( ios /= 0 ) THEN          ! not a number: thermodynamic source
      Source = arg
      iArg0  = 2
    ENDIF
  ENDIF
  IF ( nArg >= iArg0 ) THEN
    CALL GET_COMMAND_ARGUMENT( iArg0, arg ); READ( arg, * ) nThetaRef
  ENDIF
  IF ( nArg >= iArg0 + 1 ) THEN
    nTest = MIN( nArg - iArg0, nTestMax )
    DO iArg = 1, nTest
      CALL GET_COMMAND_ARGUMENT( iArg0 + iArg, arg ); READ( arg, * ) nTheta(iArg)
    ENDDO
  ENDIF
  nTheta(nTest+1) = 2 * nThetaRef   ! reference check
  iLast = nTest                     ! the run PASS/FAIL is based on
  GridMode = ( TRIM(Source) == 'grid' )

  ! ------------------------------------------------------------ setup
  CALL BuildEnergyGrid( E, dE )
  IF ( GridMode ) THEN
    CALL SetupGrid()
  ELSE
    CALL ReadThermoConditions( TRIM(Source) )
  ENDIF

  ProcList = (/ (iP, iP = 1, nP) /)
  DO iP = 1, nP
    SpecOf(iP) = SpeciesOfProcess( iP )
    Tested(iP) = ( GroupOf(iP) /= iGroupNES )
  ENDDO

  WRITE(*,'(A)') ' wlTestAngularQuadrature'
  WRITE(*,'(A,A)') '   thermodynamic points: ', TRIM(Source)
  IF ( .NOT. GridMode ) WRITE(*,'(A,ES8.1,A)') '   keeping Rho >= ', RhoMin, ' g/cm^3'
  WRITE(*,'(A,I4,A,I4,A,*(I4))') '   nE =', nE, '   reference nTheta =', nThetaRef, &
    '   tested nTheta =', nTheta(1:nTest)
  WRITE(*,'(A,I4,A,I4,A,I5)') '   nProcess =', nP, '   tested =', COUNT(Tested), &
    ' (NES only in species totals)   nX =', nX
  WRITE(*,'(A,2F9.3,A)')   '   E range  = ', E(1), E(nE), ' MeV'
  WRITE(*,'(A,2F9.2,A)')   '   T range  = ', MINVAL(T),  MAXVAL(T),  ' MeV'
  WRITE(*,'(A,2F9.2,A)')   '   Mu_e     = ', MINVAL(Me), MAXVAL(Me), ' MeV'
  WRITE(*,'(A,2F9.2,A)')   '   Mu_mu    = ', MINVAL(Mm), MAXVAL(Mm), ' MeV'
  IF ( .NOT. GridMode ) &
    WRITE(*,'(A,2ES10.2,A)') '   Rho      = ', MINVAL(Rho), MAXVAL(Rho), ' g/cm^3'

  ALLOCATE( KapRef(nE,nP,nL,nX), KapTot(nE,nSpec,nX) )
  ALLOCATE( Err(nP,nL,nX,nTest+1,nMeas), iEErr(nP,nL,nX,nTest+1) )

  ! ------------------------------------------------------------ reference
  CALL ComputeOpacities( nThetaRef, KapRef, tRun )
  WRITE(*,'(A,I4,A,F8.2,A)') '   nTheta =', nThetaRef, ' (reference): ', tRun, ' s'
  KapTot = 0.0d0
  DO iX = 1, nX
    DO iP = 1, nP
      KapTot(:,SpecOf(iP),iX) = KapTot(:,SpecOf(iP),iX) + KapRef(:,iP,1,iX)
    ENDDO
  ENDDO

  ! ------------------------------------------------------------ tested + check
  DO iTest = 1, nTest + 1
    CALL CompareToReference( nTheta(iTest), Err(:,:,:,iTest,:), iEErr(:,:,:,iTest), tRun )
    WRITE(*,'(A,I4,A,F8.2,A)') '   nTheta =', nTheta(iTest), '            : ', tRun, ' s'
  ENDDO

  ! ------------------------------------------------------------ report per process
  WRITE(*,*)
  WRITE(*,'(A,I0,A)') '   opacity error vs nTheta = ', nThetaRef, ', max over l and points'
  WRITE(*,'(A)',ADVANCE='NO') '   proc grp species      share:'
  DO iTest = 1, nTest
    WRITE(*,'(A,I3)',ADVANCE='NO') '  nTh', nTheta(iTest)
  ENDDO
  WRITE(*,'(A,I3,A,I3,A)') '  chk', nTheta(nTest+1), '   own', nTheta(iLast), &
    '  spec   worst share at: T     Rho       Mu_e  Mu_mu   E1  l'

  nFail = 0
  DO iP = 1, nP
    IF ( .NOT. Tested(iP) ) CYCLE
    WRITE(*,'(A,I5,A,A3,A,A10,A)',ADVANCE='NO') '  ', iP, '  ', GroupName(GroupOf(iP)), &
      ' ', SpecName(SpecOf(iP)), '       '
    DO iTest = 1, nTest + 1
      WRITE(*,'(ES9.1)',ADVANCE='NO') MAXVAL( Err(iP,:,:,iTest,iShare) )
    ENDDO
    WRITE(*,'(2ES9.1)',ADVANCE='NO') MAXVAL( Err(iP,:,:,iLast,iOwn) ), &
                                     MAXVAL( Err(iP,:,:,iLast,iSpec) )
    iWorst = 1; lWorst = 1; Worst = -1.0d0
    DO iX = 1, nX
      DO iL = 1, nL
        IF ( Err(iP,iL,iX,iLast,iShare) > Worst ) THEN
          Worst = Err(iP,iL,iX,iLast,iShare); iWorst = iX; lWorst = iL
        ENDIF
      ENDDO
    ENDDO
    WRITE(*,'(A,F7.1,ES10.2,2F7.1,F7.1,I3)',ADVANCE='NO') '    ', T(iWorst), Rho(iWorst), &
      Me(iWorst), Mm(iWorst), E(iEErr(iP,lWorst,iWorst,iLast)), lWorst-1
    IF ( Worst > OpacTol ) THEN
      nFail = nFail + 1
      WRITE(*,'(A)') '  <-- FAIL'
    ELSE
      WRITE(*,*)
    ENDIF
  ENDDO

  ! ------------------------------------------------------------ per-group summary
  WRITE(*,*)
  WRITE(*,'(A)') '   per group, max over processes, l and points:   own / share / spec'
  WRITE(*,'(A)',ADVANCE='NO') '   grp'
  DO iTest = 1, nTest + 1
    WRITE(*,'(A,I4,A)',ADVANCE='NO') '            nTh', nTheta(iTest), '        '
  ENDDO
  WRITE(*,*)
  DO iG = 1, nGroup
    IF ( iG == iGroupNES ) CYCLE
    WRITE(*,'(A,A3)',ADVANCE='NO') '   ', GroupName(iG)
    DO iTest = 1, nTest + 1
      DO iM = 1, nMeas
        v = 0.0d0
        DO iP = 1, nP
          IF ( GroupOf(iP) == iG ) v = MAX( v, MAXVAL( Err(iP,:,:,iTest,iM) ) )
        ENDDO
        WRITE(*,'(ES9.1)',ADVANCE='NO') v
      ENDDO
      WRITE(*,'(A)',ADVANCE='NO') '   '
    ENDDO
    WRITE(*,*)
  ENDDO

  ! ------------------------------------------------------------ per-point file
  WRITE( FileName, '(A,I0,A)' ) 'AngularQuadrature_ref', nThetaRef, '.dat'
  OPEN( NEWUNIT = iUnit, FILE = TRIM(FileName), STATUS = 'REPLACE', ACTION = 'WRITE' )
  WRITE(iUnit,'(A,I0,A,A)') '# opacity error vs nTheta = ', nThetaRef, ', points from ', TRIM(Source)
  WRITE(iUnit,'(A,*(I5))') '# blocks of 10 columns for nTheta =', nTheta(1:nTest+1)
  WRITE(iUnit,'(A)') '# T[MeV] Rho[g/cm3] Mu_e[MeV] Mu_mu[MeV]  then per block (max over l):' // &
    ' share_NMS LFE LFC IMD  spec_NMS LFE LFC IMD'
  DO iX = 1, nX
    WRITE(iUnit,'(F9.3,ES11.3,2F9.3)',ADVANCE='NO') T(iX), Rho(iX), Me(iX), Mm(iX)
    DO iTest = 1, nTest + 1
      GrpErr = 0.0d0
      DO iP = 1, nP
        IF ( .NOT. Tested(iP) ) CYCLE
        iG = GroupOf(iP)
        DO iM = 1, nMeas
          GrpErr(iG,iM) = MAX( GrpErr(iG,iM), MAXVAL( Err(iP,:,iX,iTest,iM) ) )
        ENDDO
      ENDDO
      WRITE(iUnit,'(8ES11.3)',ADVANCE='NO') GrpErr(2:,iShare), GrpErr(2:,iSpec)
    ENDDO
    WRITE(iUnit,*)
  ENDDO
  CLOSE( iUnit )

  WRITE(*,*)
  DO iTest = 1, nTest + 1
    WRITE(*,'(A,I4,A,3ES10.2)') '   overall max error (tested), nTheta =', nTheta(iTest), &
      ' (own, share, spec): ', ( MaxTested( iTest, iM ), iM = 1, nMeas )
  ENDDO
  WRITE(*,'(A,A)') '   per-point errors written to ', TRIM(FileName)
  IF ( MaxTested( nTest+1, iOwn ) > RefTol ) &
    WRITE(*,'(A,I0,A,I0,A,ES8.1)') '   WARNING: reference not converged: nTheta = ', &
      nThetaRef, ' vs ', nTheta(nTest+1), ' differ by more than ', RefTol

  WRITE(*,*)
  IF ( nFail == 0 ) THEN
    WRITE(*,'(A,I0,A,ES8.1,A)') ' wlTestAngularQuadrature: PASSED (nTheta = ', &
      nTheta(iLast), ', share error tol ', OpacTol, ')'
  ELSE
    WRITE(*,'(A,I3,A,I0,A,ES8.1,A)') ' wlTestAngularQuadrature: FAILED (', nFail, &
      ' processes above share error tol at nTheta = ', nTheta(iLast), ', tol ', OpacTol, ')'
  ENDIF

CONTAINS

  ! max of Err(:,:,:,iTest,iM) over the tested processes
  REAL(DP) FUNCTION MaxTested( iTest, iM )
    INTEGER, INTENT(IN) :: iTest, iM
    INTEGER :: iP
    MaxTested = 0.0d0
    DO iP = 1, nP
      IF ( Tested(iP) ) MaxTested = MAX( MaxTested, MAXVAL( Err(iP,:,:,iTest,iM) ) )
    ENDDO
  END FUNCTION MaxTested

  PURE INTEGER FUNCTION GroupOf( iProc )
    INTEGER, INTENT(IN) :: iProc
    SELECT CASE ( iProc )
    CASE ( 1:12  ); GroupOf = 1
    CASE ( 13:24 ); GroupOf = 2
    CASE ( 25:28 ); GroupOf = 3
    CASE ( 29:32 ); GroupOf = 4
    CASE DEFAULT  ; GroupOf = 5
    END SELECT
  END FUNCTION GroupOf

  ! Incoming neutrino species = first token of the reaction string.
  INTEGER FUNCTION SpeciesOfProcess( iProc )
    INTEGER, INTENT(IN) :: iProc
    CHARACTER(128) :: r
    INTEGER :: k
    CALL ReactionStringFromProcessIndex( iProc, r )
    r = ADJUSTL( r )
    r = r( 1 : INDEX( r, ' ' ) - 1 )
    SpeciesOfProcess = 0
    DO k = 1, nSpec
      IF ( TRIM(r) == TRIM(SpecName(k)) ) SpeciesOfProcess = k
    ENDDO
    IF ( SpeciesOfProcess == 0 ) THEN
      WRITE(*,*) 'SpeciesOfProcess: cannot parse species of process ', iProc, ': ', TRIM(r)
      STOP
    ENDIF
  END FUNCTION SpeciesOfProcess

  ! Points from a ThermoConditions-format file, Rho >= RhoMin.
  SUBROUTINE ReadThermoConditions( File )
    CHARACTER(*), INTENT(IN) :: File
    REAL(DP), ALLOCATABLE :: Row(:,:)
    INTEGER :: nRow, i, j, u
    OPEN( NEWUNIT = u, FILE = File, STATUS = 'OLD', ACTION = 'READ' )
    READ(u,*) nRow
    READ(u,*)
    ALLOCATE( Row(8,nRow) )
    DO i = 1, nRow
      READ(u,*) Row(:,i)          ! T Rho Ye Ym Mue Mum Mun Mup (rest ignored)
    ENDDO
    CLOSE( u )
    nX = COUNT( Row(2,:) >= RhoMin )
    ALLOCATE( T(nX), Rho(nX), Me(nX), Mm(nX), MuNu(nSpec,nX) )
    j = 0
    DO i = 1, nRow
      IF ( Row(2,i) < RhoMin ) CYCLE
      j = j + 1
      T(j) = Row(1,i); Rho(j) = Row(2,i); Me(j) = Row(5,i); Mm(j) = Row(6,i)
      ! chemical equilibrium: mu_nu_l = mu_l - (mu_n - mu_p)
      MuNu(1,j) =   Me(j) - ( Row(7,i) - Row(8,i) )
      MuNu(2,j) = - MuNu(1,j)
      MuNu(3,j) =   Mm(j) - ( Row(7,i) - Row(8,i) )
      MuNu(4,j) = - MuNu(3,j)
      MuNu(5,j) = 0.0d0
      MuNu(6,j) = 0.0d0
    ENDDO
    DEALLOCATE( Row )
  END SUBROUTINE ReadThermoConditions

  SUBROUTINE SetupGrid()
    INTEGER :: iT, iMe, iMm, j
    nX = nT * nMe * nMm
    ALLOCATE( T(nX), Rho(nX), Me(nX), Mm(nX), MuNu(nSpec,nX) )
    j = 0
    DO iMm = 1, nMm
      DO iMe = 1, nMe
        DO iT = 1, nT
          j = j + 1
          T(j) = T_List(iT); Me(j) = Me_List(iMe); Mm(j) = Mm_List(iMm)
        ENDDO
      ENDDO
    ENDDO
    Rho  = 0.0d0      ! not defined in grid mode
    MuNu = 0.0d0
  END SUBROUTINE SetupGrid

  ! Kernels with nTh angular nodes -> opacities Kap(E1,iP,iL,iX).
  SUBROUTINE ComputeOpacities( nTh, Kap, tSec )
    INTEGER , INTENT(IN)  :: nTh
    REAL(DP), INTENT(OUT) :: Kap(nE,nP,nL,nX)
    REAL(DP), INTENT(OUT) :: tSec
    REAL(DP) :: Ph(nP,nL), Phi(nE,nE,nP,nL)
    INTEGER  :: RowOf(nP), i, j, iP, iL, iX
    INTEGER(8) :: c0, c1, rate

    CALL InitGeneralScatteringKernels( E, E, nE, nE, nTh, &
                                       iProcessList_Option = ProcList )
    DO iP = 1, nP
      RowOf(iP) = ProcessRowIndex( iP )
    ENDDO

    CALL SYSTEM_CLOCK( c0, rate )
    !$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE( i, j, iP, iL, Ph, Phi )
    DO iX = 1, nX
      DO j = 1, nE
        DO i = 1, nE
          CALL CalculateAllPhout( i, j, E(i), E(j), T(iX), Me(iX), Mm(iX), Ph, nL )
          DO iP = 1, nP
            Phi(i,j,iP,:) = Ph(RowOf(iP),:)
          ENDDO
        ENDDO
      ENDDO
      DO iL = 1, nL
        DO iP = 1, nP
          Kap(:,iP,iL,iX) = Opacity( Phi(:,:,iP,iL) )
        ENDDO
      ENDDO
    ENDDO
    !$OMP END PARALLEL DO
    CALL SYSTEM_CLOCK( c1 )
    tSec = REAL( c1 - c0, DP ) / REAL( rate, DP )

    CALL FinalizeGeneralScatteringKernels()
  END SUBROUTINE ComputeOpacities

  SUBROUTINE CompareToReference( nTh, ErrOut, iEOut, tSec )
    INTEGER , INTENT(IN)  :: nTh
    REAL(DP), INTENT(OUT) :: ErrOut(nP,nL,nX,nMeas)
    INTEGER , INTENT(OUT) :: iEOut(nP,nL,nX)
    REAL(DP), INTENT(OUT) :: tSec
    REAL(DP), ALLOCATABLE :: Kap(:,:,:,:)
    REAL(DP) :: dK(nE), w(nE), Norm, e_i
    INTEGER  :: iP, iL, iX, iS, i

    ALLOCATE( Kap(nE,nP,nL,nX) )
    CALL ComputeOpacities( nTh, Kap, tSec )
    ErrOut = 0.0d0
    iEOut  = 1
    DO iX = 1, nX
      DO iP = 1, nP
        iS = SpecOf(iP)
        ! thermal spectrum of the incoming species: E^2 dE f_FD
        DO i = 1, nE
          w(i) = E(i)**2 * dE(i) / ( EXP( MIN( ( E(i) - MuNu(iS,iX) ) / T(iX), 700.0d0 ) ) + 1.0d0 )
        ENDDO
        DO iL = 1, nL
          Norm = REAL( 2*iL-1, DP )
          dK   = ABS( Kap(:,iP,iL,iX) - KapRef(:,iP,iL,iX) )
          DO i = 1, nE
            IF ( ABS( KapRef(i,iP,1,iX) ) >= KappaFloor ) &
              ErrOut(iP,iL,iX,iOwn) = MAX( ErrOut(iP,iL,iX,iOwn), &
                                           dK(i) / ( Norm * ABS( KapRef(i,iP,1,iX) ) ) )
            IF ( KapTot(i,iS,iX) >= KappaFloor ) THEN
              e_i = dK(i) / ( Norm * KapTot(i,iS,iX) )
              IF ( e_i > ErrOut(iP,iL,iX,iShare) ) THEN
                ErrOut(iP,iL,iX,iShare) = e_i
                iEOut (iP,iL,iX)        = i
              ENDIF
            ENDIF
          ENDDO
          IF ( SUM( w * KapTot(:,iS,iX) ) > 0.0d0 ) &
            ErrOut(iP,iL,iX,iSpec) = SUM( w * dK ) / ( Norm * SUM( w * KapTot(:,iS,iX) ) )
        ENDDO
      ENDDO
    ENDDO
    DEALLOCATE( Kap )
  END SUBROUTINE CompareToReference

  ! kappa(E1) = Sum_j Phi(E1,E3_j) E3_j^2 dE3_j   [cm^-1]
  PURE FUNCTION Opacity( Phi ) RESULT( Kap )
    REAL(DP), INTENT(IN) :: Phi(nE,nE)
    REAL(DP) :: Kap(nE)
    INTEGER  :: i
    DO i = 1, nE
      Kap(i) = SUM( Phi(i,:) * E(:)**2 * dE(:) )
    ENDDO
  END FUNCTION Opacity

  ! Node energies and widths (Gauss weight * zone width) of a thornado-style
  ! geometric grid on [eL, eR]: nZones zones (ratio ZoomE), nNodes nodes per zone.
  SUBROUTINE BuildEnergyGrid( En, dEn )
    REAL(DP), INTENT(OUT) :: En(nE), dEn(nE)
    REAL(DP) :: w1, left, width, xi(nNodes), wq(nNodes)
    INTEGER  :: k, n

    SELECT CASE ( nNodes )
    CASE ( 1 ); xi(1) = 0.0d0
                wq(1) = 1.0d0
    CASE ( 2 ); xi(1) = -0.5d0 / SQRT(3.0d0)
                xi(nNodes) = 0.5d0 / SQRT(3.0d0)
                wq = 0.5d0
    END SELECT

    w1   = ( eR - eL ) * ( ZoomE - 1.0d0 ) / ( ZoomE**nZones - 1.0d0 )
    left = eL
    DO k = 1, nZones
      width = w1 * ZoomE**(k-1)
      DO n = 1, nNodes
        En ( (k-1)*nNodes + n ) = left + width * ( 0.5d0 + xi(n) )
        dEn( (k-1)*nNodes + n ) = width * wq(n)
      ENDDO
      left = left + width
    ENDDO
  END SUBROUTINE BuildEnergyGrid

END PROGRAM wlTestAngularQuadrature
