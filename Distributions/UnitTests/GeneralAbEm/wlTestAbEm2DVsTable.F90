PROGRAM wlTestAbEm2DVsTable

  !===========================================================================
  !  Compare the 2D semi-leptonic absorption opacity
  !  (wlSemiLeptonicOpacityModule2D, Guo/Lohs formalism) with the
  !  tabulated EmAb opacity on free nucleons (Reddy et al. 1998 x
  !  Horowitz weak-magnetism factor, see TableCreator/OpacitySource/
  !  FullWeakPhysics/wlCreateOpacityTable.F90 and Chimera/
  !  abem_nucleons_Reddy_et_al_1998_weaklib.f90).
  !
  !  What is stored in the table (per species, units 1/cm):
  !    Value = ab + em = kappa_a * ( 1 + exp(eta) )
  !  i.e. the absorption opacity corrected for stimulated absorption, with
  !    nu_e    : eta = ( Mue + Mup - Mun - E ) / T
  !    nu_e_bar: eta = ( Mun - Mup - Mue - E ) / T
  !  The 2D module returns kappa_a only, so we multiply it by (1+exp(eta)).
  !
  !  Chemical potential convention (identical to the table creator):
  !    The weaklib EOS stores cmpn, cmpp WITHOUT rest mass and shifted by
  !    -dmnp (TableCreator/EOSSource/TableOperations/TableFinishing/
  !    wlTableFinisher.f90:224-225), cmpe INCLUDING the electron rest mass.
  !    The Reddy routine rebuilds the relativistic chemical potentials as
  !      Mun = cmpn + dmnp + mn ,  Mup = cmpp + dmnp + mp ,  Mue = cmpe
  !    (abem_nucleons_Reddy_et_al_1998_weaklib.f90, nu_type branches) with
  !    the Chimera constants of physcnst_module (mn, mp, dmnp). We use the
  !    same constants here. Bare masses, no mean-field potentials (U = 0).
  !
  !  Note on couplings: the Reddy routine uses ga = 1.21 and a hard-coded
  !  G2 = 6.79d-10 MeV^-5 m^-1 (abem_module_weaklib.f90, cc_difcs), while
  !  the 2D module uses ga = 1.2723 and (Gw_MeV*Vud)^2. The leading-order
  !  ratio of the two normalizations is printed and used to give a
  !  "coupling-matched" error estimate as well.
  !
  !  Two modes:
  !    UseExactTablePoints = .TRUE.  : evaluate on table grid points
  !                                    (no interpolation error).
  !    UseExactTablePoints = .FALSE. : evaluate at the (rho, T, Ye) of
  !                                    ThermoConditions.dat, interpolating
  !                                    the opacity and EOS tables.
  !===========================================================================

  USE wlKindModule,         ONLY: dp
  USE wlEosConstantsModule, ONLY: &
    me, kmev, kmev_inv, Gw_MeV, Vud, ga, gv, hbarc
  USE wlInterpolationModule, ONLY: &
    LogInterpolateSingleVariable_1D3D_Custom_Point, &
    LogInterpolateSingleVariable_3D_Custom_Point
  USE wlOpacityFieldsModule, ONLY: &
    iNu_e, iNu_e_bar
  USE wlOpacityTableModule, ONLY: &
    OpacityTableType, &
    DeAllocateOpacityTable
  USE wlIOModuleHDF, ONLY: &
    InitializeHDF, &
    FinalizeHDF
  USE wlOpacityTableIOModuleHDF, ONLY: &
    ReadOpacityTableHDF
  USE wlSemiLeptonicOpacityModule2D, ONLY: &
    Opacity_CC_2D_PrecomputedNodes, &
    gauleg
  USE wlSemiLeptonicOpacityIntegrationModule2D, ONLY: &
    NuAbsorptionOnNeutrons, &
    NuBarAbsorptionOnProtons
#if defined(WEAKLIB_OMP)
  USE omp_lib, ONLY: &
    omp_get_max_threads
#endif

  IMPLICIT NONE

  !--------- test parameters ------------------------------------------------
  INTEGER,  PARAMETER :: nE_quad             = 16      ! GL nodes per dimension
  INTEGER,  PARAMETER :: WhichCorrection     = 1       ! 0 LO, 1 +WM (as table), 3 full
  LOGICAL,  PARAMETER :: UseExactTablePoints = .FALSE.
  LOGICAL,  PARAMETER :: DoIntegrationModule = .TRUE.  ! also time/check the 2-phase wrappers
  INTEGER,  PARAMETER :: nRepeat_Table       = 100     ! table timing is averaged over repeats

  ! Exact-table-point selection
  INTEGER,  PARAMETER :: nYeSel              = 5
  REAL(dp), PARAMETER :: Ye_Targets(nYeSel)  = [ 0.1d0, 0.2d0, 0.3d0, 0.4d0, 0.5d0 ]
  INTEGER,  PARAMETER :: iRho_Stride         = 1
  REAL(dp), PARAMETER :: Rho_Min             = 1.0d10  ! g/cm^3
  INTEGER,  PARAMETER :: iT_Stride           = 1
  REAL(dp), PARAMETER :: T_Min_MeV           = 1.0d0

  ! Statistics filters
  REAL(dp), PARAMETER :: Kappa_Floor         = 1.0d-10 ! 1/cm
  REAL(dp), PARAMETER :: E_Stat_Min          = 2.0d0   ! MeV
  REAL(dp), PARAMETER :: E_Stat_Max          = 300.0d0 ! MeV
  ! "Clean" statistics exclude regimes where the Reddy routine used for the
  ! table is known to break down (see cc_difcs in abem_module_weaklib.f90):
  !  - eta > EtaCut: final lepton strongly blocked; cc_difcs sets the
  !    integrand to zero for (mue - ep)/T > 20 (efac ~ eta), and the
  !    truncated kappa_a is then amplified by (1+exp(eta)).
  !  - (Mu_N - m_N)/T > EtaNucCut: strongly degenerate nucleons; the
  !    response uses differences of approximate polylog fits (POLYLOG) and
  !    loses all accuracy.
  REAL(dp), PARAMETER :: EtaCut              = 10.0d0
  REAL(dp), PARAMETER :: EtaNucCut           = 20.0d0

  ! Chimera constants used by the table creator (physcnst_module.f90)
  REAL(dp), PARAMETER :: mn_Chimera   = 939.565d0
  REAL(dp), PARAMETER :: mp_Chimera   = 938.272d0
  REAL(dp), PARAMETER :: dmnp_Chimera = 1.29333d0
  REAL(dp), PARAMETER :: ga_Reddy     = 1.21d0
  REAL(dp), PARAMETER :: G2_Reddy     = 6.79d-10       ! MeV^-5 m^-1
  REAL(dp), PARAMETER :: EtaMax       = 300.0d0        ! same cap as fexp

  CHARACTER(LEN=256)  :: OpTableName          = 'wl-Op-SFHo-15-25-50-E40-EmAb.h5'
  CHARACTER(LEN=256)  :: EosTableName         = 'wl-EOS-SFHo-15-25-50.h5'
  CHARACTER(LEN=256)  :: ThermoConditionsName = 'ThermoConditions.dat'

  !--------- tables ---------------------------------------------------------
  TYPE(OpacityTableType) :: OpacityTable
  INTEGER  :: icmpe, icmpn, icmpp
  INTEGER  :: iRho, iT, iYe, iEOS_Rho, iEOS_T, iEOS_Ye
  REAL(dp) :: Offset_cmpe, Offset_cmpn, Offset_cmpp, Offset_EmAb(2)

  !--------- thermodynamic points ------------------------------------------
  INTEGER  :: nPoints, nE
  REAL(dp), ALLOCATABLE :: Energy(:), LogEs(:), LogDs(:), LogTs(:)
  REAL(dp), ALLOCATABLE :: P_Rho(:), P_T(:), P_Ye(:)             ! T in K
  REAL(dp), ALLOCATABLE :: P_Mue(:), P_Mun(:), P_Mup(:)          ! MeV, incl. rest mass
  REAL(dp), ALLOCATABLE :: Kappa_Table(:,:,:)                    ! (nE, nPoints, 2)
  REAL(dp), ALLOCATABLE :: Kappa_2D(:,:,:), KappaA_2D(:,:,:), KappaA_Int(:,:,:)
  REAL(dp), ALLOCATABLE :: EtaNu(:,:,:)                            ! degeneracy of the eq. neutrino

  !--------- misc -----------------------------------------------------------
  REAL(dp) :: xa(nE_quad), wxa(nE_quad)
  INTEGER  :: IncludeCorrections(3)
  INTEGER  :: i, iE, iS, j, k, l, ierr, nRead, iYeSel(nYeSel)
  REAL(dp) :: TMeV, eta, buffer, rho_in, T_in, Ye_in, cmp(3)
  REAL(dp) :: CouplingRatio, MaxDiffInt, GridErr
  INTEGER  :: count_start, count_end, count_rate
  REAL(dp) :: t_2D, t_Int, t_Table, MaxDiffInterp, nEval
  INTEGER  :: iRep, nThreads
  REAL(dp), ALLOCATABLE :: Kappa_Interp(:,:,:), Cmp_Interp(:,:)
  INTEGER  :: nNaN(2)
  REAL(dp), ALLOCATABLE :: Rho_tmp(:), T_tmp(:), Ye_tmp(:)

  CALL MPI_INIT( ierr )

  !------------------------------------------------------
  !    read in the opacity table (+ its EOS table)
  !------------------------------------------------------
  CALL InitializeHDF( )
  CALL ReadOpacityTableHDF( OpacityTable,         &
         FileName_EmAb_Option = TRIM(OpTableName), &
         EquationOfStateTableName_Option = TRIM(EosTableName), &
         Verbose_Option = .TRUE. )
  CALL FinalizeHDF( )

  iRho     = OpacityTable % TS % Indices % iRho
  iT       = OpacityTable % TS % Indices % iT
  iYe      = OpacityTable % TS % Indices % iYe
  iEOS_Rho = OpacityTable % EOSTable % TS % Indices % iRho
  iEOS_T   = OpacityTable % EOSTable % TS % Indices % iT
  iEOS_Ye  = OpacityTable % EOSTable % TS % Indices % iYe

  icmpe = OpacityTable % EOSTable % DV % Indices % iElectronChemicalPotential
  icmpn = OpacityTable % EOSTable % DV % Indices % iNeutronChemicalPotential
  icmpp = OpacityTable % EOSTable % DV % Indices % iProtonChemicalPotential
  Offset_cmpe = OpacityTable % EOSTable % DV % Offsets(icmpe)
  Offset_cmpn = OpacityTable % EOSTable % DV % Offsets(icmpn)
  Offset_cmpp = OpacityTable % EOSTable % DV % Offsets(icmpp)
  Offset_EmAb = OpacityTable % EmAb % Offsets

  WRITE(*,'(A,2ES12.4)') ' EmAb offsets (nue, nuebar): ', Offset_EmAb
  WRITE(*,'(A,3ES12.4)') ' EOS offsets (cmpe,cmpn,cmpp): ', &
    Offset_cmpe, Offset_cmpn, Offset_cmpp

  nE = OpacityTable % EnergyGrid % nPoints
  ALLOCATE( Energy(nE), LogEs(nE) )
  Energy = OpacityTable % EnergyGrid % Values
  LogEs  = LOG10( Energy )

  ! The table creator reads the EOS DVs with the opacity grid indices,
  ! so the two thermo-state grids must coincide.
  GridErr = 0.0d0
  GridErr = MAX( GridErr, MaxRelDiff( OpacityTable % TS % States(iRho) % Values, &
                                      OpacityTable % EOSTable % TS % States(iEOS_Rho) % Values ) )
  GridErr = MAX( GridErr, MaxRelDiff( OpacityTable % TS % States(iT) % Values, &
                                      OpacityTable % EOSTable % TS % States(iEOS_T) % Values ) )
  GridErr = MAX( GridErr, MaxRelDiff( OpacityTable % TS % States(iYe) % Values, &
                                      OpacityTable % EOSTable % TS % States(iEOS_Ye) % Values ) )
  WRITE(*,'(A,ES12.4)') ' Max rel. difference opacity vs EOS TS grids: ', GridErr
  IF( GridErr > 1.0d-6 )THEN
    WRITE(*,*) 'ERROR: opacity and EOS thermo-state grids do not match.'
    STOP
  END IF

  !------------------------------------------------------
  !    build the list of thermodynamic points
  !------------------------------------------------------
  ASSOCIATE &
  ( Rhos => OpacityTable % TS % States(iRho) % Values, &
    Ts   => OpacityTable % TS % States(iT)   % Values, &
    Yes  => OpacityTable % TS % States(iYe)  % Values, &
    EOS_cmpe => OpacityTable % EOSTable % DV % Variables(icmpe) % Values, &
    EOS_cmpn => OpacityTable % EOSTable % DV % Variables(icmpn) % Values, &
    EOS_cmpp => OpacityTable % EOSTable % DV % Variables(icmpp) % Values, &
    TableNue    => OpacityTable % EmAb % Opacity(iNu_e)     % Values, &
    TableNuebar => OpacityTable % EmAb % Opacity(iNu_e_bar) % Values )

  ALLOCATE( LogDs(SIZE(Rhos)), LogTs(SIZE(Ts)) )
  LogDs = LOG10( Rhos )
  LogTs = LOG10( Ts )

  IF( UseExactTablePoints )THEN

    DO l = 1, nYeSel
      iYeSel(l) = MINLOC( ABS( Yes - Ye_Targets(l) ), DIM = 1 )
    END DO

    nPoints = 0
    DO l = 1, nYeSel
      DO k = 1, SIZE(Ts), iT_Stride
        IF( Ts(k) * kmev < T_Min_MeV ) CYCLE
        DO j = 1, SIZE(Rhos), iRho_Stride
          IF( Rhos(j) < Rho_Min ) CYCLE
          nPoints = nPoints + 1
        END DO
      END DO
    END DO

    CALL AllocatePointArrays( nPoints, nE )

    i = 0
    DO l = 1, nYeSel
      DO k = 1, SIZE(Ts), iT_Stride
        IF( Ts(k) * kmev < T_Min_MeV ) CYCLE
        DO j = 1, SIZE(Rhos), iRho_Stride
          IF( Rhos(j) < Rho_Min ) CYCLE
          i = i + 1
          P_Rho(i) = Rhos(j)
          P_T  (i) = Ts(k)
          P_Ye (i) = Yes(iYeSel(l))
          ! Undo the log10(value + offset) storage exactly (no interpolation)
          cmp(1) = 10.0d0**EOS_cmpe(j,k,iYeSel(l)) - Offset_cmpe
          cmp(2) = 10.0d0**EOS_cmpn(j,k,iYeSel(l)) - Offset_cmpn
          cmp(3) = 10.0d0**EOS_cmpp(j,k,iYeSel(l)) - Offset_cmpp
          CALL SetChemicalPotentials( i, cmp )
          Kappa_Table(:,i,iNu_e)     = 10.0d0**TableNue   (:,j,k,iYeSel(l)) - Offset_EmAb(iNu_e)
          Kappa_Table(:,i,iNu_e_bar) = 10.0d0**TableNuebar(:,j,k,iYeSel(l)) - Offset_EmAb(iNu_e_bar)
        END DO
      END DO
    END DO

  ELSE

    ! ThermoConditions.dat: T[MeV] Rho Ye Ym Mue Mum Mun Mup Un Up M*n M*p
    ! Only (T, Rho, Ye) are used; chemical potentials are interpolated from
    ! the table's own EOS to stay consistent with the tabulated opacity.
    OPEN( UNIT = 123, FILE = TRIM(ThermoConditionsName), STATUS = 'OLD', ACTION = 'READ' )
    READ(123,*) nRead
    READ(123,*)
    ALLOCATE( Rho_tmp(nRead), T_tmp(nRead), Ye_tmp(nRead) )
    nPoints = 0
    DO i = 1, nRead
      READ(123,*) T_in, rho_in, Ye_in, buffer
      IF( rho_in < MAX( Rho_Min, Rhos(1) ) .OR. rho_in > Rhos(SIZE(Rhos)) ) CYCLE
      IF( T_in < MAX( T_Min_MeV, Ts(1)*kmev ) .OR. T_in > Ts(SIZE(Ts))*kmev ) CYCLE
      IF( Ye_in < Yes(1) .OR. Ye_in > Yes(SIZE(Yes)) ) CYCLE
      nPoints = nPoints + 1
      Rho_tmp(nPoints) = rho_in
      T_tmp  (nPoints) = T_in * kmev_inv
      Ye_tmp (nPoints) = Ye_in
    END DO
    CLOSE(123)
    WRITE(*,'(A,I6,A,I6)') ' ThermoConditions points inside table range: ', nPoints, ' / ', nRead

    CALL AllocatePointArrays( nPoints, nE )
    P_Rho = Rho_tmp(1:nPoints)
    P_T   = T_tmp  (1:nPoints)
    P_Ye  = Ye_tmp (1:nPoints)
    DEALLOCATE( Rho_tmp, T_tmp, Ye_tmp )

  END IF

  !------------------------------------------------------
  !    table approach (timed): interpolate EOS chemical
  !    potentials + EmAb opacities at every point
  !------------------------------------------------------
  ! This is what a transport code does when using the table, so it is
  ! timed for both point selections. At exact table points it must
  ! reproduce the stored values; otherwise it provides Kappa_Table.
  ALLOCATE( Kappa_Interp(nE,nPoints,2), Cmp_Interp(3,nPoints) )

  WRITE(*,*) 'Interpolating EOS + EmAb tables ...'
  CALL SYSTEM_CLOCK( count_start, count_rate = count_rate )
  DO iRep = 1, nRepeat_Table
#if defined(WEAKLIB_OMP)
    !$OMP PARALLEL DO SCHEDULE(static)
#endif
    DO i = 1, nPoints
      CALL LogInterpolateSingleVariable_3D_Custom_Point &
             ( P_Rho(i), P_T(i), P_Ye(i), Rhos, Ts, Yes, Offset_cmpe, EOS_cmpe, Cmp_Interp(1,i) )
      CALL LogInterpolateSingleVariable_3D_Custom_Point &
             ( P_Rho(i), P_T(i), P_Ye(i), Rhos, Ts, Yes, Offset_cmpn, EOS_cmpn, Cmp_Interp(2,i) )
      CALL LogInterpolateSingleVariable_3D_Custom_Point &
             ( P_Rho(i), P_T(i), P_Ye(i), Rhos, Ts, Yes, Offset_cmpp, EOS_cmpp, Cmp_Interp(3,i) )
      CALL LogInterpolateSingleVariable_1D3D_Custom_Point &
             ( LogEs, LOG10(P_Rho(i)), LOG10(P_T(i)), P_Ye(i), LogEs, LogDs, LogTs, Yes, &
               Offset_EmAb(iNu_e), TableNue, Kappa_Interp(:,i,iNu_e) )
      CALL LogInterpolateSingleVariable_1D3D_Custom_Point &
             ( LogEs, LOG10(P_Rho(i)), LOG10(P_T(i)), P_Ye(i), LogEs, LogDs, LogTs, Yes, &
               Offset_EmAb(iNu_e_bar), TableNuebar, Kappa_Interp(:,i,iNu_e_bar) )
    END DO
#if defined(WEAKLIB_OMP)
    !$OMP END PARALLEL DO
#endif
  END DO
  CALL SYSTEM_CLOCK( count_end )
  t_Table = REAL( count_end - count_start, dp ) / REAL( count_rate, dp ) &
          / REAL( nRepeat_Table, dp )

  MaxDiffInterp = 0.0d0
  IF( UseExactTablePoints )THEN
    DO iS = 1, 2
      DO i = 1, nPoints
        DO iE = 1, nE
          IF( Kappa_Table(iE,i,iS) <= Kappa_Floor ) CYCLE
          MaxDiffInterp = MAX( MaxDiffInterp, &
            ABS( Kappa_Interp(iE,i,iS) - Kappa_Table(iE,i,iS) ) / Kappa_Table(iE,i,iS) )
        END DO
      END DO
    END DO
  ELSE
    Kappa_Table = Kappa_Interp
    DO i = 1, nPoints
      CALL SetChemicalPotentials( i, Cmp_Interp(:,i) )
    END DO
  END IF

  END ASSOCIATE

  CALL DeAllocateOpacityTable( OpacityTable )

  WRITE(*,'(A,I6,A,I3,A)') ' Number of thermodynamic points: ', nPoints, &
                           ' x ', nE, ' energies'

  !------------------------------------------------------
  !    2D semi-leptonic opacity
  !------------------------------------------------------
  CALL gauleg( 0.0d0, 1.0d0, xa, wxa, nE_quad )

  SELECT CASE( WhichCorrection )
  CASE (0) ; IncludeCorrections = [ 0, 0, 0 ]
  CASE (1) ; IncludeCorrections = [ 1, 0, 0 ]
  CASE (2) ; IncludeCorrections = [ 1, 1, 0 ]
  CASE DEFAULT ; IncludeCorrections = [ 1, 1, 1 ]
  END SELECT

  WRITE(*,*) 'Computing 2D opacities (Opacity_CC_2D_PrecomputedNodes) ...'
  CALL SYSTEM_CLOCK( count_start, count_rate = count_rate )
#if defined(WEAKLIB_OMP)
  !$OMP PARALLEL DO SCHEDULE(dynamic) COLLAPSE(2) PRIVATE( TMeV )
#endif
  DO i = 1, nPoints
    DO iE = 1, nE
      TMeV = P_T(i) * kmev
      CALL Opacity_CC_2D_PrecomputedNodes &
             ( WhichCorrection, 1, Energy(iE), KappaA_2D(iE,i,iNu_e), TMeV, &
               P_Mue(i), P_Mun(i), P_Mup(i), me, mn_Chimera, mp_Chimera, &
               0.0d0, 0.0d0, nE_quad, xa, wxa )
      CALL Opacity_CC_2D_PrecomputedNodes &
             ( WhichCorrection, 2, Energy(iE), KappaA_2D(iE,i,iNu_e_bar), TMeV, &
               P_Mue(i), P_Mun(i), P_Mup(i), me, mn_Chimera, mp_Chimera, &
               0.0d0, 0.0d0, nE_quad, xa, wxa )
    END DO
  END DO
#if defined(WEAKLIB_OMP)
  !$OMP END PARALLEL DO
#endif
  CALL SYSTEM_CLOCK( count_end )
  t_2D = REAL( count_end - count_start, dp ) / REAL( count_rate, dp )

  t_Int      = 0.0d0
  MaxDiffInt = 0.0d0
  IF( DoIntegrationModule )THEN
    WRITE(*,*) 'Computing 2D opacities (integration module wrappers) ...'
    CALL SYSTEM_CLOCK( count_start, count_rate = count_rate )
#if defined(WEAKLIB_OMP)
    !$OMP PARALLEL DO SCHEDULE(dynamic) COLLAPSE(2) PRIVATE( TMeV )
#endif
    DO i = 1, nPoints
      DO iE = 1, nE
        TMeV = P_T(i) * kmev
        CALL NuAbsorptionOnNeutrons &
               ( Energy(iE), TMeV, P_Mue(i), P_Mun(i), P_Mup(i), me, &
                 mn_Chimera, mp_Chimera, 0.0d0, 0.0d0, IncludeCorrections, &
                 nE_quad, xa, wxa, xa, wxa, KappaA_Int(iE,i,iNu_e) )
        CALL NuBarAbsorptionOnProtons &
               ( Energy(iE), TMeV, P_Mue(i), P_Mun(i), P_Mup(i), me, &
                 mn_Chimera, mp_Chimera, 0.0d0, 0.0d0, IncludeCorrections, &
                 nE_quad, xa, wxa, xa, wxa, KappaA_Int(iE,i,iNu_e_bar) )
      END DO
    END DO
#if defined(WEAKLIB_OMP)
    !$OMP END PARALLEL DO
#endif
    CALL SYSTEM_CLOCK( count_end )
    t_Int = REAL( count_end - count_start, dp ) / REAL( count_rate, dp )

    DO iS = 1, 2
      DO i = 1, nPoints
        DO iE = 1, nE
          IF( IsBad( KappaA_2D(iE,i,iS) ) .OR. IsBad( KappaA_Int(iE,i,iS) ) ) CYCLE
          IF( KappaA_2D(iE,i,iS) <= 0.0d0 ) CYCLE
          MaxDiffInt = MAX( MaxDiffInt, &
            ABS( KappaA_Int(iE,i,iS) - KappaA_2D(iE,i,iS) ) / KappaA_2D(iE,i,iS) )
        END DO
      END DO
    END DO
  END IF

  !------------------------------------------------------
  !    stimulated absorption: kappa* = kappa_a * (1 + exp(eta))
  !------------------------------------------------------
  nNaN = 0
  DO i = 1, nPoints
    TMeV = P_T(i) * kmev
    DO iE = 1, nE
      DO iS = 1, 2
        IF( iS == iNu_e )THEN
          eta = ( P_Mue(i) + P_Mup(i) - P_Mun(i) - Energy(iE) ) / TMeV
        ELSE
          eta = ( P_Mun(i) - P_Mup(i) - P_Mue(i) - Energy(iE) ) / TMeV
        END IF
        EtaNu(iE,i,iS) = eta
        eta = MIN( EtaMax, MAX( -EtaMax, eta ) )
        IF( IsBad( KappaA_2D(iE,i,iS) ) )THEN
          nNaN(iS) = nNaN(iS) + 1
          Kappa_2D(iE,i,iS) = -1.0d0          ! flag invalid sample
        ELSE
          Kappa_2D(iE,i,iS) = KappaA_2D(iE,i,iS) * ( 1.0d0 + EXP( eta ) )
        END IF
      END DO
    END DO
  END DO

  !------------------------------------------------------
  !    output
  !------------------------------------------------------
  OPEN( UNIT = 66, FILE = 'AbEm2D_Table.dat',   STATUS = 'REPLACE', ACTION = 'WRITE' )
  OPEN( UNIT = 77, FILE = 'AbEm2D_General.dat', STATUS = 'REPLACE', ACTION = 'WRITE' )
  WRITE(66,'(A)') '# iE, E[MeV], rho[g/cc], T[MeV], Ye, kappa*_nue[1/cm], kappa*_nuebar[1/cm]  (table)'
  WRITE(77,'(A)') '# iE, E[MeV], rho[g/cc], T[MeV], Ye, kappa*_nue[1/cm], kappa*_nuebar[1/cm]  (2D, -1 = NaN)'
  DO i = 1, nPoints
    DO iE = 1, nE
      WRITE(66,'(I4,1X,6(ES17.8E3,1X))') iE, Energy(iE), P_Rho(i), P_T(i)*kmev, P_Ye(i), &
        Kappa_Table(iE,i,iNu_e), Kappa_Table(iE,i,iNu_e_bar)
      WRITE(77,'(I4,1X,6(ES17.8E3,1X))') iE, Energy(iE), P_Rho(i), P_T(i)*kmev, P_Ye(i), &
        Kappa_2D(iE,i,iNu_e), Kappa_2D(iE,i,iNu_e_bar)
    END DO
  END DO
  CLOSE(66)
  CLOSE(77)

  !------------------------------------------------------
  !    summary
  !------------------------------------------------------
  ! Leading-order (non-relativistic) normalization ratio Reddy-table / 2D:
  !   G2 (gv^2 + 3 ga_R^2)  vs  (Gw Vud)^2/(hbar c) (gv^2 + 3 ga^2)
  CouplingRatio = G2_Reddy * ( 1.0d0 + 3.0d0 * ga_Reddy**2 ) &
                / ( ( Gw_MeV * Vud )**2 / ( hbarc * 1.0d-2 ) * ( gv**2 + 3.0d0 * ga**2 ) )

  WRITE(*,*)
  WRITE(*,'(A)') ' ================= wlTestAbEm2DVsTable summary ================='
  WRITE(*,'(A,L2,A,I2,A,I4)') ' UseExactTablePoints =', UseExactTablePoints, &
    '   WhichCorrection =', WhichCorrection, '   nE_quad =', nE_quad
  nThreads = 1
#if defined(WEAKLIB_OMP)
  nThreads = omp_get_max_threads()
#endif
  ! one evaluation = one (point, energy, species) opacity
  nEval = REAL( nPoints, dp ) * REAL( nE, dp ) * 2.0d0
  WRITE(*,'(A,I4,A,ES10.3)') ' Timing (OpenMP threads =', nThreads, &
    '), evaluations (points x energies x species) = ', nEval
  WRITE(*,'(A)') '                                              wall [s]   per eval [s]   / table'
  WRITE(*,'(A,ES12.4,ES14.4,ES12.4)') ' Table interpolation (EOS mu + EmAb)    : ', &
    t_Table, t_Table / nEval, 1.0d0
  WRITE(*,'(A,ES12.4,ES14.4,ES12.4)') ' 2D (Opacity_CC_2D_PrecomputedNodes)    : ', &
    t_2D, t_2D / nEval, t_2D / MAX( t_Table, TINY(1.0d0) )
  IF( DoIntegrationModule )THEN
    WRITE(*,'(A,ES12.4,ES14.4,ES12.4)') ' 2D (integration module wrappers)       : ', &
      t_Int, t_Int / nEval, t_Int / MAX( t_Table, TINY(1.0d0) )
    WRITE(*,'(A,ES12.4)')   ' Max rel. diff integration module vs 2D  : ', MaxDiffInt
  END IF
  WRITE(*,'(A,I4,A)') ' (table time averaged over ', nRepeat_Table, ' repeats)'
  IF( UseExactTablePoints ) &
    WRITE(*,'(A,ES12.4)') ' Max rel. diff table interpolation vs stored values at nodes: ', &
      MaxDiffInterp
  WRITE(*,'(A,F8.4)') ' LO coupling ratio Reddy-table / 2D (ga=1.21,G2 vs ga=1.2723,GwVud): ', &
    CouplingRatio
  WRITE(*,'(A,ES9.2,A,F6.1,A,F6.1,A)') ' Stats for kappa_table > ', Kappa_Floor, &
    ' /cm and ', E_Stat_Min, ' <= E <= ', E_Stat_Max, ' MeV'
  WRITE(*,'(A,F5.1,A,F5.1,A)') ' clean: drop samples with eta > ', EtaCut, &
    ' or (Mu_N-m_N)/T > ', EtaNucCut, ' (Reddy routine unreliable, see header)'

  DO iS = 1, 2
    IF( iS == iNu_e )THEN
      WRITE(*,'(/,A,I8)') ' --- nu_e  (nu_e + n -> e- + p) ---   NaN samples: ', nNaN(iS)
    ELSE
      WRITE(*,'(/,A,I8)') ' --- nu_e_bar (nu_e_bar + p -> e+ + n) ---   NaN samples: ', nNaN(iS)
    END IF
    CALL PrintStatistics( iS, 1.0d0,         .FALSE., 'all            ' )
    CALL PrintStatistics( iS, 1.0d0,         .TRUE.,  'clean (eta, eta_N cuts)' )
    CALL PrintStatistics( iS, CouplingRatio, .TRUE.,  'clean, 2D x LO coupling ratio' )
    CALL PrintDecadeTable( iS )
  END DO

  CALL PrintSamples( 1.0d13, 10.0d0 )
  CALL PrintSamples( 1.0d11,  3.0d0 )

  DEALLOCATE( Energy, LogEs, LogDs, LogTs )
  DEALLOCATE( P_Rho, P_T, P_Ye, P_Mue, P_Mun, P_Mup )
  DEALLOCATE( Kappa_Table, Kappa_2D, KappaA_2D, KappaA_Int, EtaNu )
  DEALLOCATE( Kappa_Interp, Cmp_Interp )

  CALL MPI_FINALIZE( ierr )

CONTAINS

  SUBROUTINE AllocatePointArrays( nP, nEn )
    INTEGER, INTENT(IN) :: nP, nEn
    ALLOCATE( P_Rho(nP), P_T(nP), P_Ye(nP), P_Mue(nP), P_Mun(nP), P_Mup(nP) )
    ALLOCATE( Kappa_Table(nEn,nP,2), Kappa_2D(nEn,nP,2) )
    ALLOCATE( KappaA_2D(nEn,nP,2), KappaA_Int(nEn,nP,2), EtaNu(nEn,nP,2) )
    KappaA_Int = 0.0d0
  END SUBROUTINE AllocatePointArrays


  ! cmp = (cmpe, cmpn, cmpp) in the weaklib/Chimera EOS convention
  SUBROUTINE SetChemicalPotentials( ip, cmpIn )
    INTEGER,  INTENT(IN) :: ip
    REAL(dp), INTENT(IN) :: cmpIn(3)
    P_Mue(ip) = cmpIn(1)
    P_Mun(ip) = cmpIn(2) + dmnp_Chimera + mn_Chimera
    P_Mup(ip) = cmpIn(3) + dmnp_Chimera + mp_Chimera
  END SUBROUTINE SetChemicalPotentials


  LOGICAL FUNCTION IsBad( x )
    REAL(dp), INTENT(IN) :: x
    IsBad = ( x /= x ) .OR. ( ABS(x) > HUGE(x) )
  END FUNCTION IsBad


  REAL(dp) FUNCTION MaxRelDiff( a, b )
    REAL(dp), INTENT(IN) :: a(:), b(:)
    IF( SIZE(a) /= SIZE(b) )THEN
      MaxRelDiff = HUGE(1.0d0)
    ELSE
      MaxRelDiff = MAXVAL( ABS( a - b ) / MAX( ABS(b), TINY(1.0d0) ) )
    END IF
  END FUNCTION MaxRelDiff


  LOGICAL FUNCTION UseSample( iEn, ip, iSp, CutEta )
    INTEGER, INTENT(IN) :: iEn, ip, iSp
    LOGICAL, INTENT(IN) :: CutEta
    UseSample = Energy(iEn) >= E_Stat_Min .AND. Energy(iEn) <= E_Stat_Max &
          .AND. Kappa_Table(iEn,ip,iSp) > Kappa_Floor &
          .AND. Kappa_2D(iEn,ip,iSp) >= 0.0d0
    IF( CutEta ) UseSample = UseSample .AND. EtaNu(iEn,ip,iSp) <= EtaCut &
      .AND. MAX( P_Mun(ip) - mn_Chimera, P_Mup(ip) - mp_Chimera ) / ( P_T(ip) * kmev ) <= EtaNucCut
  END FUNCTION UseSample


  ! Stats of |kappa_2D * Scale / kappa_table - 1| over the filtered samples
  SUBROUTINE PrintStatistics( iSp, Scale, CutEta, Label )
    INTEGER,          INTENT(IN) :: iSp
    REAL(dp),         INTENT(IN) :: Scale
    LOGICAL,          INTENT(IN) :: CutEta
    CHARACTER(LEN=*), INTENT(IN) :: Label
    REAL(dp), ALLOCATABLE :: AbsErr(:), Ratio(:)
    INTEGER :: n, ip, iEn
    n = 0
    ALLOCATE( AbsErr(nE*nPoints), Ratio(nE*nPoints) )
    DO ip = 1, nPoints
      DO iEn = 1, nE
        IF( .NOT. UseSample( iEn, ip, iSp, CutEta ) ) CYCLE
        n = n + 1
        Ratio (n) = Kappa_2D(iEn,ip,iSp) * Scale / Kappa_Table(iEn,ip,iSp)
        AbsErr(n) = ABS( Ratio(n) - 1.0d0 )
      END DO
    END DO
    IF( n == 0 )THEN
      WRITE(*,'(2X,A,A)') Label, ': no valid samples'
    ELSE
      WRITE(*,'(2X,A)') Label
      WRITE(*,'(4X,A,I7,A,ES10.3,A,ES10.3,A,ES10.3,A,F8.4)') &
        'N =', n, '  max|rel| =', MAXVAL(AbsErr(1:n)), &
        '  median|rel| =', Median(AbsErr(1:n)), '  mean|rel| =', SUM(AbsErr(1:n))/n, &
        '  median(2D/tab) =', Median(Ratio(1:n))
    END IF
    DEALLOCATE( AbsErr, Ratio )
  END SUBROUTINE PrintStatistics


  ! Relative error binned by density decade and by energy band
  ! ("clean" samples; nCut = samples removed by the eta / eta_N cuts)
  SUBROUTINE PrintDecadeTable( iSp )
    INTEGER, INTENT(IN) :: iSp
    INTEGER, PARAMETER  :: nBand = 3
    REAL(dp), PARAMETER :: EBand(0:nBand) = [ 2.0d0, 10.0d0, 40.0d0, 300.0d0 ]
    REAL(dp), ALLOCATABLE :: AbsErr(:), Ratio(:)
    INTEGER :: iDec, iB, n, nCut, ip, iEn, iDecMin, iDecMax
    iDecMin = FLOOR( LOG10( MINVAL(P_Rho) ) )
    iDecMax = FLOOR( LOG10( MAXVAL(P_Rho) ) )
    ALLOCATE( AbsErr(nE*nPoints), Ratio(nE*nPoints) )
    WRITE(*,'(2X,A)') 'log10(rho)  E-band[MeV]      N   nCut  median(2D/tab)  median|rel|   max|rel|'
    DO iDec = iDecMin, iDecMax
      DO iB = 1, nBand
        n = 0
        nCut = 0
        DO ip = 1, nPoints
          IF( FLOOR( LOG10( P_Rho(ip) ) ) /= iDec ) CYCLE
          DO iEn = 1, nE
            IF( Energy(iEn) < EBand(iB-1) .OR. Energy(iEn) >= EBand(iB) ) CYCLE
            IF( .NOT. UseSample( iEn, ip, iSp, .FALSE. ) ) CYCLE
            IF( .NOT. UseSample( iEn, ip, iSp, .TRUE. ) )THEN
              nCut = nCut + 1
              CYCLE
            END IF
            n = n + 1
            Ratio (n) = Kappa_2D(iEn,ip,iSp) / Kappa_Table(iEn,ip,iSp)
            AbsErr(n) = ABS( Ratio(n) - 1.0d0 )
          END DO
        END DO
        IF( n == 0 ) CYCLE
        WRITE(*,'(4X,I2,A,I3,A,F5.0,A,F5.0,2I7,F14.4,2ES13.3)') iDec, '-', iDec+1, &
          '   ', EBand(iB-1), '-', EBand(iB), n, nCut, Median(Ratio(1:n)), &
          Median(AbsErr(1:n)), MAXVAL(AbsErr(1:n))
      END DO
    END DO
    DEALLOCATE( AbsErr, Ratio )
  END SUBROUTINE PrintDecadeTable


  ! Print a few energies at the point closest to (RhoTarget, TTarget[MeV])
  SUBROUTINE PrintSamples( RhoTarget, TTarget )
    REAL(dp), INTENT(IN) :: RhoTarget, TTarget
    INTEGER  :: ip, iEn, iSp, k2
    INTEGER, PARAMETER :: nSample = 3
    REAL(dp), PARAMETER :: ETargets(nSample) = [ 5.0d0, 15.0d0, 40.0d0 ]
    REAL(dp) :: TM
    ip = MINLOC( ABS( LOG10( P_Rho / RhoTarget ) ) + ABS( LOG10( P_T * kmev / TTarget ) ), DIM = 1 )
    TM = P_T(ip) * kmev
    WRITE(*,'(/,A,ES11.4,A,F8.3,A,F6.3)') ' Sample point: rho =', P_Rho(ip), &
      '  T[MeV] =', TM, '  Ye =', P_Ye(ip)
    WRITE(*,'(A,3F12.4)') '   Mue, Mun-mn, Mup-mp [MeV] = ', &
      P_Mue(ip), P_Mun(ip) - mn_Chimera, P_Mup(ip) - mp_Chimera
    WRITE(*,'(A)') '   species  E[MeV]      eta     kappa_a(2D)  1+exp(eta)   kappa*(2D)   kappa*(table)  2D/table'
    DO iSp = 1, 2
      DO k2 = 1, nSample
        iEn = MINLOC( ABS( Energy - ETargets(k2) ), DIM = 1 )
        IF( Kappa_2D(iEn,ip,iSp) < 0.0d0 .OR. KappaA_2D(iEn,ip,iSp) <= 0.0d0 )THEN
          WRITE(*,'(4X,I2,F11.3,A)') iSp, Energy(iEn), '   invalid (NaN or zero)'
          CYCLE
        END IF
        WRITE(*,'(4X,I2,F11.3,F10.3,4ES13.4,F10.4)') iSp, Energy(iEn), EtaNu(iEn,ip,iSp), KappaA_2D(iEn,ip,iSp), &
          Kappa_2D(iEn,ip,iSp) / KappaA_2D(iEn,ip,iSp), Kappa_2D(iEn,ip,iSp), &
          Kappa_Table(iEn,ip,iSp), Kappa_2D(iEn,ip,iSp) / Kappa_Table(iEn,ip,iSp)
      END DO
    END DO
  END SUBROUTINE PrintSamples


  REAL(dp) FUNCTION Median( x )
    REAL(dp), INTENT(IN) :: x(:)
    REAL(dp) :: y(SIZE(x)), tmp
    INTEGER  :: n, gap, ii, jj
    n = SIZE(x)
    y = x
    ! Shell sort
    gap = n / 2
    DO WHILE( gap > 0 )
      DO ii = gap + 1, n
        tmp = y(ii)
        jj  = ii
        DO WHILE( jj > gap )
          IF( y(jj-gap) <= tmp ) EXIT
          y(jj) = y(jj-gap)
          jj = jj - gap
        END DO
        y(jj) = tmp
      END DO
      gap = gap / 2
    END DO
    IF( MOD(n,2) == 1 )THEN
      Median = y( (n+1)/2 )
    ELSE
      Median = 0.5d0 * ( y(n/2) + y(n/2+1) )
    END IF
  END FUNCTION Median

END PROGRAM wlTestAbEm2DVsTable
