PROGRAM wlTestEmAb2DParallel

  ! Timing and consistency of the MOMENTUM integration method
  ! (OpacitySource/GeneralEmAb) in serial and with OpenMP: Opacity_CC_2D
  ! after InitOpacity_CC_2D(nE_2D, iEmAb2D_Momentum), the per-reaction
  ! wrappers and the collapsed two-phase loop of
  ! wlSemiLeptonicOpacityIntegrationModule2D (the GPU driver of that
  ! method). All of them read the quadrature rules computed once by Init.
  ! For electrons and muons it checks that the OpenMP, wrapper and
  ! collapsed results agree with the serial Opacity_CC_2D (relative
  ! difference |a-b|/max(|a|,|b|) above 1d-8 is reported; the serial,
  ! OpenMP and wrapper results are expected to be bit-identical, the
  ! collapsed sum only to round-off since it sums in another order).

  USE omp_lib
  USE wlKindModule, ONLY: dp
  USE wlSemiLeptonicOpacityModule2D, ONLY: &
    Opacity_CC_2D, InitOpacity_CC_2D, FinalizeOpacity_CC_2D, iEmAb2D_Momentum
  USE wlSemiLeptonicOpacityIntegrationModule2D, ONLY: &
    NuAbsorptionOnNeutrons, NuBarAbsorptionOnProtons, &
    InverseNeutronDecay, NeutronEmissivity, &
    CC2D_CallData, CC2D_RowData, &
    CC2D_PrepareCall, CC2D_RowSetup, CC2D_NodeKernel, &
    CC2D_ApplyPrefactor, nE3Max
  USE wlEosConstantsModule, ONLY: &
   pi, Gw_MeV, ga, gv, mn, mp, me, mmu, mpi, &
   Vud, massA, massV, gamma_p, gamma_n

  implicit none

  INTEGER, PARAMETER :: NP = 60, nOp = 4, nApprox = 4, nE_2D = 40
  INTEGER, PARAMETER :: nPointsMax = 200 ! up to 1000 in ThermoConditions.dat (~5 min on 4 threads)
  REAL(DP), PARAMETER :: masse = me , massm = mmu

  REAL(DP) :: EnuA(NP)
  REAL(DP) :: xTem, cheml, chemn, chemp, xUn, xUp, massl, xMassn, xMassp
  INTEGER :: i, j, k, l, nThermoPoints, iE2, iE3, iLep
  INTEGER :: nBadPar, nBadWrap, nBadColl, nExactPar, nExactWrap
  REAL(DP), PARAMETER :: RelTol = 1.0d-8
  CHARACTER(LEN=8) :: LepName

  REAL(DP), allocatable :: OpaA_2D(:,:,:,:), OpaA_2D_parallel(:,:,:,:)
  REAL(DP), allocatable :: OpaA_2D_wrapped(:,:,:,:), OpaA_2D_collapsed(:,:,:,:)
  REAL(DP), ALLOCATABLE :: T(:), Rho(:), Ye(:), Ym(:), Mue(:), Mum(:), Mul(:), &
      Mun(:), Mup(:), Un(:), Up(:), EffMassn(:), EffMassp(:)
  TYPE(CC2D_CallData) :: CD
  TYPE(CC2D_RowData)  :: RD(nE_2D)

  REAL(DP), allocatable :: Store2DIntegral(:,:,:)
  REAL(DP) :: OpaA_dummy
  REAL(DP) :: t_tot, t_start, t_end, t_serial
  REAL(DP) :: t_2D, t_2D_serial, t_2D_wrapped, t_2D_collapsed
  INTEGER :: IncludeCorrections(3)
  INTEGER :: count_start, count_interm, count_end, count_rate

  ! Initialize neutrino energies
  DO l = 1, NP
     EnuA(l) = l * 5.d0
  END DO

  ! Read thermodynamic data
  OPEN(unit=123, file=trim(adjustl('ThermoConditions.dat')), status='old', action='read')
  read(123,*) nThermoPoints
  read(123,*)

  ! You can also set nThermoPoints to a smaller value for quick checks
  nThermoPoints = MIN(nThermoPoints, nPointsMax)
  ALLOCATE(OpaA_2D(NP, nThermoPoints, nApprox, nOp))
  ALLOCATE(OpaA_2D_parallel(NP, nThermoPoints, nApprox, nOp))
  ALLOCATE(OpaA_2D_wrapped(NP, nThermoPoints, nApprox, nOp))
  ALLOCATE(OpaA_2D_collapsed(NP, nThermoPoints, nApprox, nOp))
  ALLOCATE(T(nThermoPoints), Rho(nThermoPoints), Ye(nThermoPoints), Ym(nThermoPoints))
  ALLOCATE(Mue(nThermoPoints), Mum(nThermoPoints), Mul(nThermoPoints), Mun(nThermoPoints), Mup(nThermoPoints))
  ALLOCATE(Un(nThermoPoints), Up(nThermoPoints), EffMassn(nThermoPoints), EffMassp(nThermoPoints))
  ALLOCATE(Store2DIntegral(nE_2D, nE3Max, nOp))

  ! Columns: T Rho Ye Ym Mue Mum Mun Mup Un Up M*n M*p
  DO i = 1, nThermoPoints
     read(123,*) T(i), Rho(i), Ye(i), Ym(i), Mue(i), Mum(i), Mun(i), Mup(i), Un(i), Up(i), &
        EffMassn(i), EffMassp(i)
  END DO
  CLOSE(123)

  ! Quadrature rules of the momentum method, computed once
  CALL InitOpacity_CC_2D( nE_2D, iEmAb2D_Momentum )

  DO iLep = 1, 2

  IF (iLep == 2) THEN
    Mul(:) = Mum(:)
    massl = massm
    LepName = 'muons'
  ELSE
    Mul(:) = Mue(:)
    massl = masse
    LepName = 'electrons'
  END IF

  OpaA_2D = 0.d0

  CALL SYSTEM_CLOCK(count_start, count_rate=count_rate)
  DO j = 1, nApprox
  DO k = 1, nOp
  DO i = 1, nThermoPoints
  DO l=1, NP
    CALL Opacity_CC_2D(j-1, k, EnuA(l), OpaA_2D(l, i, j, k), &
          T(i), Mul(i), Mun(i), Mup(i), massl, EffMassn(i), EffMassp(i), Un(i), Up(i))
  END DO
  END DO
  END DO
  END DO
  CALL SYSTEM_CLOCK(count_end)
  t_2D_serial = REAL(count_end - count_start) / REAL(count_rate)

  CALL SYSTEM_CLOCK(count_interm)
  !$OMP PARALLEL DO COLLAPSE(4) SCHEDULE(DYNAMIC) &
  !$OMP PRIVATE(i, j, k, l)
  DO j = 1, nApprox
  DO k = 1, nOp
  DO i = 1, nThermoPoints
  DO l=1, NP
    CALL Opacity_CC_2D(j-1, k, EnuA(l), OpaA_2D_parallel(l, i, j, k), &
          T(i), Mul(i), Mun(i), Mup(i), massl, EffMassn(i), EffMassp(i), Un(i), Up(i))
  END DO
  END DO
  END DO
  END DO
  !$OMP END PARALLEL DO
  CALL SYSTEM_CLOCK(count_end)
  t_2D = REAL(count_end - count_interm) / REAL(count_rate)

  CALL SYSTEM_CLOCK(count_interm)
  !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(DYNAMIC) &
  !$OMP PRIVATE(i, j, k, l, IncludeCorrections)
  DO j = 1, nApprox
  DO i = 1, nThermoPoints
  DO l = 1, NP
    SELECT CASE(j-1)
    CASE (0)
      IncludeCorrections(1) = 0
      IncludeCorrections(2) = 0
      IncludeCorrections(3) = 0
    CASE (1)
      IncludeCorrections(1) = 1
      IncludeCorrections(2) = 0
      IncludeCorrections(3) = 0
    CASE (2)
      IncludeCorrections(1) = 1
      IncludeCorrections(2) = 1
      IncludeCorrections(3) = 0
    CASE (3)
      IncludeCorrections(1) = 1
      IncludeCorrections(2) = 1
      IncludeCorrections(3) = 1
    END SELECT

    ! NOTE: wrappers now take physical arguments (as Opacity_CC_2D); the
    ! per-reaction particle mapping happens inside the integration module
    CALL NuAbsorptionOnNeutrons(EnuA(l), T(i), Mul(i), Mun(i), Mup(i), massl, &
        EffMassn(i), EffMassp(i), Un(i), Up(i), &
        IncludeCorrections, OpaA_2D_wrapped(l, i, j, 1))
    CALL NuBarAbsorptionOnProtons(EnuA(l), T(i), Mul(i), Mun(i), Mup(i), massl, &
        EffMassn(i), EffMassp(i), Un(i), Up(i), &
        IncludeCorrections, OpaA_2D_wrapped(l, i, j, 2))
    CALL InverseNeutronDecay(EnuA(l), T(i), Mul(i), Mun(i), Mup(i), massl, &
        EffMassn(i), EffMassp(i), Un(i), Up(i), &
        IncludeCorrections, OpaA_2D_wrapped(l, i, j, 3))
    CALL NeutronEmissivity(EnuA(l), T(i), Mul(i), Mun(i), Mup(i), massl, &
        EffMassn(i), EffMassp(i), Un(i), Up(i), &
        IncludeCorrections, OpaA_2D_wrapped(l, i, j, 4))
  END DO
  END DO
  END DO
  !$OMP END PARALLEL DO
  CALL SYSTEM_CLOCK(count_end)
  t_2D_wrapped = REAL(count_end - count_interm) / REAL(count_rate)

  CALL SYSTEM_CLOCK(count_interm)
  OpaA_2D_collapsed = 0.0d0
  !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(DYNAMIC) &
  !$OMP PRIVATE(i, j, k, l, iE2, iE3, Store2DIntegral, IncludeCorrections, CD, RD)
  DO i = 1, nThermoPoints
  DO j = 1, nApprox
  DO l = 1, NP

    SELECT CASE(j-1)
    CASE (0)
      IncludeCorrections(1) = 0
      IncludeCorrections(2) = 0
      IncludeCorrections(3) = 0
    CASE (1)
      IncludeCorrections(1) = 1
      IncludeCorrections(2) = 0
      IncludeCorrections(3) = 0
    CASE (2)
      IncludeCorrections(1) = 1
      IncludeCorrections(2) = 1
      IncludeCorrections(3) = 0
    CASE (3)
      IncludeCorrections(1) = 1
      IncludeCorrections(2) = 1
      IncludeCorrections(3) = 1
    END SELECT

    ! Two-phase structure: per-call setup, per-row setup, then the
    ! (iE2,iE3) node kernel exposed as a collapsed loop (GPU pattern)
    DO k = 1, nOp
      CALL CC2D_PrepareCall(k, EnuA(l), T(i), Mul(i), Mun(i), Mup(i), massl, &
          EffMassn(i), EffMassp(i), Un(i), Up(i), IncludeCorrections, CD)

      DO iE2 = 1, nE_2D
        CALL CC2D_RowSetup(CD, iE2, RD(iE2))
      END DO

      Store2DIntegral(:,:,k) = 0.0d0
      !$OMP PARALLEL DO COLLAPSE(2) &
      !$OMP PRIVATE( iE2, iE3 )
      DO iE2=1,nE_2D
      DO iE3=1,nE3Max
        CALL CC2D_NodeKernel(CD, RD(iE2), iE3, Store2DIntegral(iE2,iE3,k))
      END DO
      END DO
      !$OMP END PARALLEL DO

      OpaA_2D_collapsed(l, i, j, k) = SUM(Store2DIntegral(:,:,k))
      CALL CC2D_ApplyPrefactor(CD, OpaA_2D_collapsed(l, i, j, k))
    END DO

  END DO
  END DO
  END DO
  !$OMP END PARALLEL DO
  CALL SYSTEM_CLOCK(count_end)
  t_2D_collapsed = REAL(count_end - count_interm) / REAL(count_rate)

  t_tot = REAL(count_end - count_start) / REAL(count_rate)

   WRITE(*,'(/,2A)')      ' Leptons           : ', TRIM(LepName)
   WRITE(*,'(A,f10.3)')   'Total wall‑clock time :', t_tot
   WRITE(*,'(A,f10.3)')   ' Serial 2D         :', t_2D_serial
   WRITE(*,'(A,f10.3)')   ' PARALLEL 2D       :', t_2D
   WRITE(*,'(A,f10.3)')   ' PARALLEL 2D wrap  :', t_2D_wrapped
   WRITE(*,'(A,f10.3)')   ' PARALLEL 2D coll  :', t_2D_collapsed
   WRITE(*,'(A,i10)')     ' Threads           :', omp_get_max_threads()
   WRITE(*,'(A,f10.3)')   ' Ratio 2D          :', t_2D_serial/t_2D

  nBadPar = 0; nBadWrap = 0; nBadColl = 0; nExactPar = 0; nExactWrap = 0
  DO j = 1, nApprox
  DO i = 1, nThermoPoints
  DO l = 1, NP
  DO k = 1, nOp

    IF (OpaA_2D(l, i, j, k) == OpaA_2D_parallel(l, i, j, k)) nExactPar = nExactPar + 1
    IF (OpaA_2D(l, i, j, k) == OpaA_2D_wrapped(l, i, j, k))  nExactWrap = nExactWrap + 1
    IF (RelDiff(OpaA_2D(l, i, j, k), OpaA_2D_parallel(l, i, j, k)) > RelTol) THEN
      nBadPar = nBadPar + 1
      WRITE(*,*) l,i,j,k,'Parallel wrong', OpaA_2D(l, i, j, k), OpaA_2D_parallel(l, i, j, k)
    ENDIF
    IF (RelDiff(OpaA_2D(l, i, j, k), OpaA_2D_wrapped(l, i, j, k)) > RelTol) THEN
      nBadWrap = nBadWrap + 1
      WRITE(*,*) l,i,j,k,'Wrapped wrong', OpaA_2D(l, i, j, k), OpaA_2D_wrapped(l, i, j, k)
    ENDIF
    IF (RelDiff(OpaA_2D(l, i, j, k), OpaA_2D_collapsed(l, i, j, k)) > RelTol) THEN
      nBadColl = nBadColl + 1
      WRITE(*,*) l,i,j,k,'Collapsed wrong', OpaA_2D(l, i, j, k), OpaA_2D_collapsed(l, i, j, k)
    ENDIF
  ENDDO
  ENDDO
  ENDDO
  ENDDO
  WRITE(*,'(A,I9)')      ' Samples           :', NP*nThermoPoints*nApprox*nOp
  WRITE(*,'(A,3I9)')     ' Mismatches > tol  (parallel, wrapped, collapsed):', nBadPar, nBadWrap, nBadColl
  WRITE(*,'(A,2I9)')     ' Bit-identical     (parallel, wrapped)           :', nExactPar, nExactWrap

  END DO   ! iLep

  !-------------------------------------------------------------
  ! Persist results
  !-------------------------------------------------------------
  ! OPEN(200, file='OpaA_2D.bin', status='replace', form='unformatted', access='stream')
  ! WRITE(200) OpaA_2D
  ! CLOSE(200)

  DEALLOCATE(OpaA_2D, OpaA_2D_parallel, OpaA_2D_wrapped, OpaA_2D_collapsed)
  DEALLOCATE(Store2DIntegral)
  CALL FinalizeOpacity_CC_2D()
  DEALLOCATE(T, Rho, Ye, Ym, Mue, Mum, Mul, Mun, Mup, Un, Up, EffMassn, EffMassp)

CONTAINS

  ! |a-b|/max(|a|,|b|), 0 if both are 0
  REAL(DP) FUNCTION RelDiff( a, b )
    REAL(DP), INTENT(IN) :: a, b
    IF( MAX( ABS(a), ABS(b) ) > 0.0d0 ) THEN
      RelDiff = ABS( a - b )/MAX( ABS(a), ABS(b) )
    ELSE
      RelDiff = 0.0d0
    END IF
  END FUNCTION RelDiff

END PROGRAM wlTestEmAb2DParallel

