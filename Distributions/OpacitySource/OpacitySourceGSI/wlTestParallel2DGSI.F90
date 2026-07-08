PROGRAM wlTestParallel2DGSI

  USE omp_lib
  USE wlKindModule, ONLY: dp
  USE wlSemiLeptonicOpacityModule2D, ONLY: &
    Opacity_CC_2D
  USE wlSemiLeptonicOpacityIntegrationModule2D, ONLY: &
    NuAbsorptionOnNeutrons, NuBarAbsorptionOnProtons, &
    InverseNeutronDecay, NeutronEmissivity, &
    CC2D_CallData, CC2D_RowData, &
    CC2D_PrepareCall, CC2D_RowSetup, CC2D_NodeKernel, CC2D_ApplyPrefactor, &
    gauleg
  USE wlEosConstantsModule, ONLY: &
   pi, Gw_MeV, ga, gv, mn, mp, me, mmu, mpi, &
   Vud, massA, massV, gamma_p, gamma_n

  implicit none

  INTEGER, PARAMETER :: NP = 60, nOp = 4, nApprox = 4, nE_2D = 40
  REAL(DP), PARAMETER :: masse = me , massm = mmu

  REAL(DP) :: EnuA(NP)
  REAL(DP) :: xTem, cheml, chemn, chemp, xUn, xUp, massl, xMassn, xMassp
  INTEGER :: i, j, k, l, nThermoPoints, iE2, iE3
  LOGICAL, PARAMETER :: DoMuons = .true.

  REAL(DP), allocatable :: OpaA_2D(:,:,:,:), OpaA_2D_parallel(:,:,:,:)
  REAL(DP), allocatable :: OpaA_2D_wrapped(:,:,:,:), OpaA_2D_collapsed(:,:,:,:)
  REAL(DP), ALLOCATABLE :: T(:), Rho(:), Ye(:), Ym(:), Mue(:), Mum(:), Mul(:), &
      Mun(:), Mup(:), Un(:), Up(:), EffMassn(:), EffMassp(:)
  TYPE(CC2D_CallData) :: CD
  TYPE(CC2D_RowData)  :: RD(nE_2D)

  REAL(DP), allocatable :: xa(:), wxa(:)
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
  nThermoPoints = 1000
  ALLOCATE(OpaA_2D(NP, nThermoPoints, nApprox, nOp))
  ALLOCATE(OpaA_2D_parallel(NP, nThermoPoints, nApprox, nOp))
  ALLOCATE(OpaA_2D_wrapped(NP, nThermoPoints, nApprox, nOp))
  ALLOCATE(OpaA_2D_collapsed(NP, nThermoPoints, nApprox, nOp))
  ALLOCATE(T(nThermoPoints), Rho(nThermoPoints), Ye(nThermoPoints), Ym(nThermoPoints))
  ALLOCATE(Mue(nThermoPoints), Mum(nThermoPoints), Mul(nThermoPoints), Mun(nThermoPoints), Mup(nThermoPoints))
  ALLOCATE(Un(nThermoPoints), Up(nThermoPoints), EffMassn(nThermoPoints), EffMassp(nThermoPoints))
  ALLOCATE(Store2DIntegral(nE_2D, nE_2D, nOp))
  ALLOCATE(xa(nE_2D))
  ALLOCATE(wxa(nE_2D))

  DO i = 1, nThermoPoints
     read(123,*) T(i), Rho(i), Ye(i), Ym(i), Mue(i), Mum(i), Mun(i), Mup(i), Un(i), Up(i)
  END DO
  CLOSE(123)

  IF (DoMuons) THEN
    Mul(:) = Mum(:)
    massl = massm
  ELSE
    Mul(:) = Mue(:)
    massl = masse
  END IF

  OpaA_2D = 0.d0

  CALL SYSTEM_CLOCK(count_start, count_rate=count_rate)
  DO j = 1, nApprox
  DO k = 1, nOp
  DO i = 1, nThermoPoints
  DO l=1, NP
    CALL Opacity_CC_2D(j-1, k, EnuA(l), OpaA_2D(l, i, j, k), &
          T(i), Mul(i), Mun(i), Mup(i), massl, EffMassn(i), EffMassp(i), Un(i), Up(i), nE_2D)
  END DO
  END DO
  END DO
  END DO
  CALL SYSTEM_CLOCK(count_end)
  t_2D_serial = REAL(count_end - count_start) / REAL(count_rate)

  CALL SYSTEM_CLOCK(count_interm)
  !$OMP PARALLEL DO COLLAPSE(4) &
  !$OMP PRIVATE(i, j, k, l)
  DO j = 1, nApprox
  DO k = 1, nOp
  DO i = 1, nThermoPoints
  DO l=1, NP
    CALL Opacity_CC_2D(j-1, k, EnuA(l), OpaA_2D_parallel(l, i, j, k), &
          T(i), Mul(i), Mun(i), Mup(i), massl, EffMassn(i), EffMassp(i), Un(i), Up(i), nE_2D)
  END DO
  END DO
  END DO
  END DO
  !$OMP END PARALLEL DO
  CALL SYSTEM_CLOCK(count_end)
  t_2D = REAL(count_end - count_interm) / REAL(count_rate)

  CALL gauleg(0.d0, 1.d0, xa, wxa, nE_2D)

  CALL SYSTEM_CLOCK(count_interm)
  !$OMP PARALLEL DO COLLAPSE(3) &
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
        IncludeCorrections, nE_2D, xa, wxa, xa, wxa, OpaA_2D_wrapped(l, i, j, 1))
    CALL NuBarAbsorptionOnProtons(EnuA(l), T(i), Mul(i), Mun(i), Mup(i), massl, &
        EffMassn(i), EffMassp(i), Un(i), Up(i), &
        IncludeCorrections, nE_2D, xa, wxa, xa, wxa, OpaA_2D_wrapped(l, i, j, 2))
    CALL InverseNeutronDecay(EnuA(l), T(i), Mul(i), Mun(i), Mup(i), massl, &
        EffMassn(i), EffMassp(i), Un(i), Up(i), &
        IncludeCorrections, nE_2D, xa, wxa, xa, wxa, OpaA_2D_wrapped(l, i, j, 3))
    CALL NeutronEmissivity(EnuA(l), T(i), Mul(i), Mun(i), Mup(i), massl, &
        EffMassn(i), EffMassp(i), Un(i), Up(i), &
        IncludeCorrections, nE_2D, xa, wxa, xa, wxa, OpaA_2D_wrapped(l, i, j, 4))
  END DO
  END DO
  END DO
  !$OMP END PARALLEL DO
  CALL SYSTEM_CLOCK(count_end)
  t_2D_wrapped = REAL(count_end - count_interm) / REAL(count_rate)

  CALL SYSTEM_CLOCK(count_interm)
  OpaA_2D_collapsed = 0.0d0
  !$OMP PARALLEL DO COLLAPSE(3) &
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
        CALL CC2D_RowSetup(CD, xa(iE2), wxa(iE2), RD(iE2))
      END DO

      Store2DIntegral(:,:,k) = 0.0d0
      !$OMP PARALLEL DO COLLAPSE(2) &
      !$OMP PRIVATE( iE2, iE3 )
      DO iE2=1,nE_2D
      DO iE3=1,nE_2D
        CALL CC2D_NodeKernel(CD, RD(iE2), xa(iE3), wxa(iE3), Store2DIntegral(iE2,iE3,k))
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

   WRITE(*,'(/,A,f10.3)') 'Total wall‑clock time :', t_tot
   WRITE(*,'(A,f10.3)')   ' Serial 2D         :', t_2D_serial
   WRITE(*,'(A,f10.3)')   ' PARALLEL 2D       :', t_2D
   WRITE(*,'(A,f10.3)')   ' PARALLEL 2D wrap  :', t_2D_wrapped
   WRITE(*,'(A,f10.3)')   ' PARALLEL 2D coll  :', t_2D_collapsed
   WRITE(*,'(A,i10)')     ' Threads           :', omp_get_max_threads()
   WRITE(*,'(A,f10.3)')   ' Ratio 2D          :', t_2D_serial/t_2D

  DO j = 1, nApprox
  DO i = 1, nThermoPoints
  DO l = 1, NP
  DO k = 1, nOp

    IF (OpaA_2D(l, i, j, k) .ne. OpaA_2D_wrapped(l, i, j, k)) THEN
      IF ( ABS(OpaA_2D(l, i, j, k) - OpaA_2D_wrapped(l, i, j, k)) &
          / OpaA_2D(l, i, j, k) > 1.0d-8 ) THEN
        WRITE(*,*) l,i,j,'Wrapped wrong'
        WRITE(*,*) OpaA_2D(l, i, j, k), OpaA_2D_wrapped(l, i, j, k), &
                  OpaA_2D_wrapped(l, i, j, k) / OpaA_2D(l, i, j, k)
      ENDIF
    ENDIF
    IF (OpaA_2D(l, i, j, k) .ne. OpaA_2D_collapsed(l, i, j, k)) THEN
      IF ( ABS(OpaA_2D(l, i, j, k) - OpaA_2D_collapsed(l, i, j, k)) &
          / OpaA_2D(l, i, j, k) > 1.0d-8 ) THEN
        WRITE(*,*) l,i,j,'Collapsed wrong'
        WRITE(*,*) OpaA_2D(l, i, j, k), OpaA_2D_collapsed(l, i, j, k), &
                  OpaA_2D_collapsed(l, i, j, k) / OpaA_2D(l, i, j, k)
      ENDIF
    ENDIF
  ENDDO
  ENDDO
  ENDDO
  ENDDO

  !-------------------------------------------------------------
  ! Persist results
  !-------------------------------------------------------------
  ! OPEN(200, file='OpaA_2D.bin', status='replace', form='unformatted', access='stream')
  ! WRITE(200) OpaA_2D
  ! CLOSE(200)

  DEALLOCATE(OpaA_2D, OpaA_2D_collapsed)
  DEALLOCATE(T, Rho, Ye, Ym, Mue, Mum, Mul, Mun, Mup, Un, Up, EffMassn, EffMassp)

END PROGRAM wlTestParallel2DGSI

