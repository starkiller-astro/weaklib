PROGRAM wlTestParallel2DGSI

  USE omp_lib
  USE wlKindModule, ONLY: dp
  USE wlSemiLeptonicOpacityModule2D, ONLY: &
    Opacity_CC_2D
  USE wlSemiLeptonicOpacityIntegrationModule2D, ONLY: &
    NuAbsorptionOnNeutrons, NuBarAbsorptionOnProtons, &
    InverseNeutronDecay, NeutronEmissivity, &
    NuCapture2DIntegration, NuBarCapture2DIntegration, &
    InverseNeutronDecay2DIntegration, NeutronEmissivity2DIntegration, &
    gauleg
  USE wlEosConstantsModule, ONLY: &
   pi, Gw_MeV, ga, gv, mn, mp, me, mmu, mpi, &
   Vud, massA, massV, gamma_p, gamma_n

  implicit none

  INTEGER, PARAMETER :: NP = 60, nOp = 4, nApprox = 4, nE_2D = 40
  REAL(DP), PARAMETER :: masse = me , massm = mmu
  REAL(DP), PARAMETER :: massn = mn, massp = mp

  REAL(DP) :: EnuA(NP)
  REAL(DP) :: xTem, cheml, chemn, chemp, xUn, xUp, massl
  INTEGER :: i, j, k, l, nThermoPoints, iE2, iE3
  LOGICAL, PARAMETER :: DoMuons = .true.

  REAL(DP), allocatable :: OpaA_2D(:,:,:,:), OpaA_2D_parallel(:,:,:,:)
  REAL(DP), allocatable :: OpaA_2D_wrapped(:,:,:,:), OpaA_2D_collapsed(:,:,:,:)
  REAL(DP), allocatable :: T(:), Rho(:), Ye(:), Ym(:), &
                         Mue(:), Mum(:), Mul(:), Mun(:), Mup(:), Un(:), Up(:)
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
  ALLOCATE(Un(nThermoPoints), Up(nThermoPoints))
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
          T(i), Mul(i), Mun(i), Mup(i), massl, massn, massp, Un(i), Up(i), nE_2D)
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
          T(i), Mul(i), Mun(i), Mup(i), massl, massn, massp, Un(i), Up(i), nE_2D)
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

    CALL NuAbsorptionOnNeutrons(EnuA(l), T(i), massn, massl, massp, &
        Un(i), Up(i), Mun(i), Mul(i), Mup(i), &
        IncludeCorrections, nE_2D, xa, wxa, xa, wxa, OpaA_2D_wrapped(l, i, j, 1))
    CALL NuBarAbsorptionOnProtons(EnuA(l), T(i), massp, massl, massn, &
        Up(i), Un(i), Mup(i), Mul(i), Mun(i), &
        IncludeCorrections, nE_2D, xa, wxa, xa, wxa, OpaA_2D_wrapped(l, i, j, 2))
    CALL InverseNeutronDecay(EnuA(l), T(i), massp, massl, massn, &
        Up(i), Un(i), Mup(i), Mul(i), Mun(i), &
        IncludeCorrections, nE_2D, xa, wxa, xa, wxa, OpaA_2D_wrapped(l, i, j, 3))
    CALL NeutronEmissivity(EnuA(l), T(i), massp, massl, massn, &
        Up(i), Un(i), Mup(i), Mul(i), Mun(i), &
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
  !$OMP PRIVATE(i, j, k, l, iE2, iE3, Store2DIntegral, IncludeCorrections)
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

    Store2DIntegral(:,:,:) = 0.0d0
    !$OMP PARALLEL DO COLLAPSE(2) &
    !$OMP PRIVATE( iE2, iE3 )
    DO iE2=1,nE_2D
    DO iE3=1,nE_2D
      
      CALL NuCapture2DIntegration(EnuA(l), T(i), massn, massl, massp, &
          Un(i), Up(i), Mun(i), Mul(i), Mup(i), &
          IncludeCorrections, xa(iE2), wxa(iE2), xa(iE3), wxa(iE3), Store2DIntegral(iE2,iE3,1))
      CALL NuBarCapture2DIntegration(EnuA(l), T(i), massp, massl, massn, &
          Up(i), Un(i), Mup(i), Mul(i), Mun(i), &
          IncludeCorrections, xa(iE2), wxa(iE2), xa(iE3), wxa(iE3), Store2DIntegral(iE2,iE3,2))
      CALL InverseNeutronDecay2DIntegration(EnuA(l), T(i), massp, massl, massn, &
          Up(i), Un(i), Mup(i), Mul(i), Mun(i), &
          IncludeCorrections, xa(iE2), wxa(iE2), xa(iE3), wxa(iE3), Store2DIntegral(iE2,iE3,3))
      CALL NeutronEmissivity2DIntegration(EnuA(l), T(i), massp, massl, massn, &
          Up(i), Un(i), Mup(i), Mul(i), Mun(i), &
          IncludeCorrections, xa(iE2), wxa(iE2), xa(iE3), wxa(iE3), Store2DIntegral(iE2,iE3,4))
    END DO
    END DO
    !$OMP END PARALLEL DO
    OpaA_2D_collapsed(l, i, j, 1) = SUM(Store2DIntegral(:,:,1))
    OpaA_2D_collapsed(l, i, j, 2) = SUM(Store2DIntegral(:,:,2))
    OpaA_2D_collapsed(l, i, j, 3) = SUM(Store2DIntegral(:,:,3))
    OpaA_2D_collapsed(l, i, j, 4) = SUM(Store2DIntegral(:,:,4))
  
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
    k = 1

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

  !-------------------------------------------------------------
  ! Persist results
  !-------------------------------------------------------------
  ! OPEN(200, file='OpaA_2D.bin', status='replace', form='unformatted', access='stream')
  ! WRITE(200) OpaA_2D
  ! CLOSE(200)

  DEALLOCATE(OpaA_2D, OpaA_2D_collapsed)
  DEALLOCATE(T, Rho, Ye, Ym, Mue, Mum, Mul, Mun, Mup, Un, Up)

END PROGRAM wlTestParallel2DGSI

