PROGRAM wlTestAbEm4DVsGSI

  ! Compares the weaklib 2D (Opacity_CC_2D) and 4D CUBA (Opacity_CC_4D)
  ! kernels with the original GSI 2D/4D codes. Requires USE_CUBA=TRUE.
  !
  ! Units: weaklib kernels return cm^-1 (reactions 1-3) and
  ! (s cm^3 MeV)^-1 (reaction 4). The GSI codes return km^-1 for
  ! reactions 1-3 (converted here) and (s cm^3 MeV)^-1 for reaction 4
  ! (no conversion).
  !
  ! NOTE: Integration_D in wlSemiLeptonicOpacityModule4D currently applies
  ! an extra 1/1.0d5 to the reaction-4 emissivity, so 4D-vs-GSI4D for
  ! reaction 4 is expected to be off by a factor 1e-5 until that is fixed.
  !
  ! The GSI codes return G <= 0 for reactions 3,4 (|M|^2 clipped before the
  ! crossing sign). Sign mismatches are counted and reported; relative
  ! errors are computed on absolute values where both are nonzero.

  USE wlKindModule, ONLY: dp
  USE wlSemiLeptonicOpacityModule2D, ONLY: &
    Opacity_CC_2D, InitOpacity_CC_2D
  USE wlSemiLeptonicOpacityModule4D, ONLY: &
    Opacity_CC_4D
  USE wlEosConstantsModule, ONLY: &
   pi, Gw_MeV, ga, gv, mn, mp, me, mmu, mpi, &
   Vud, massA, massV, gamma_p, gamma_n

  IMPLICIT NONE

  INTEGER, PARAMETER :: NP = 60, nOp = 4, nApprox = 4, nE_2D = 50
  REAL(DP), PARAMETER :: masse = me , massm = mmu

  REAL(DP) :: EnuA(NP), Error(NP)
  REAL(DP) :: xTem, cheml, chemn, chemp, xUn, xUp, massl, xMassn, xMassp
  INTEGER :: i, j, k, l, nThermoPoints
  INTEGER :: nDone, total, lastPrint
  LOGICAL, PARAMETER :: DoMuons = .true.

  REAL(DP), ALLOCATABLE :: OpaA_2D(:,:,:,:), OpaA_4D(:,:,:,:)
  REAL(DP), ALLOCATABLE :: OpaA_2D_OLD(:,:,:,:), OpaA_4D_OLD(:,:,:,:)
  REAL(DP), ALLOCATABLE :: T(:), Rho(:), Ye(:), Ym(:), Mue(:), Mum(:), &
      Mun(:), Mup(:), Un(:), Up(:), EffMassn(:), EffMassp(:)

  REAL(DP) :: t1, t2, t_start, t_end, t_2D, t_4D
  REAL(DP) :: t_2D_OLD, t_4D_OLD, err
  REAL(DP) :: errMax2D(nApprox, nOp), errMax4D(nApprox, nOp)
  INTEGER  :: nSign2D(nApprox, nOp), nSign4D(nApprox, nOp)
  INTEGER  :: nNeg2D(nApprox, nOp)
  REAL(DP), PARAMETER :: km_to_cm_inv = 1.0d-5

  ! GSI Constants, they will change the results!
  REAL(DP) , PARAMETER :: pi_GSI = 3.1415927d0, Gw_MeV_GSI = 1.166d-11, &
      gA_GSI = 1.2723d0, gV_GSI=1.d0, Mnp_GSI = 938.919d0, &
      Dnp_GSI = 1.293d0, massA_GSI = 1.0d3, massV_GSI = 840.0d0, &
      F2wm0_GSI = 3.706d0

  WRITE(*,*) 'CONSTANTS REL ERR'
  WRITE(*,*) 'pi', ABS(pi - pi_GSI)/pi_GSI
  WRITE(*,*) 'Gw_MeV', ABS(Gw_MeV - Gw_MeV_GSI)/Gw_MeV_GSI
  WRITE(*,*) 'ga', ABS(ga - gA_GSI)/gA_GSI
  WRITE(*,*) 'Mnp', ABS(Mnp_GSI - (mn + mp)/2.0d0)/Mnp_GSI
  WRITE(*,*) 'Dnp', ABS(Dnp_GSI - (mn - mp))/Dnp_GSI
  WRITE(*,*) 'massA', ABS(massA_GSI - massA)/massA_GSI
  WRITE(*,*) 'massV', ABS(massV_GSI - massV)/massV_GSI
  WRITE(*,*) 'F2wm0', ABS(F2wm0_GSI - (gamma_p - gamma_n - 1.0d0))/F2wm0_GSI

  ! Initialize neutrino energies
  DO l = 1, NP
    EnuA(l) = l * 5.d0
  END DO

  ! READ thermodynamic data
  OPEN(unit=123, file=trim(adjustl('ThermoConditions.dat')), status='old', action='READ')
  READ(123,*) nThermoPoints
  READ(123,*)

  ! You can also set nThermoPoints to a smaller value for quick checks
  nThermoPoints = 1
  ALLOCATE(OpaA_2D(NP, nThermoPoints, nApprox, nOp))
  ALLOCATE(OpaA_4D(NP, nThermoPoints, nApprox, nOp))
  ALLOCATE(OpaA_2D_OLD(NP, nThermoPoints, nApprox, nOp))
  ALLOCATE(OpaA_4D_OLD(NP, nThermoPoints, nApprox, nOp))
  ALLOCATE(T(nThermoPoints), Rho(nThermoPoints), Ye(nThermoPoints), Ym(nThermoPoints))
  ALLOCATE(Mue(nThermoPoints), Mum(nThermoPoints), Mun(nThermoPoints), Mup(nThermoPoints))
  ALLOCATE(Un(nThermoPoints), Up(nThermoPoints), EffMassn(nThermoPoints), EffMassp(nThermoPoints))

  DO i = 1, nThermoPoints
    READ(123,*) T(i), Rho(i), Ye(i), Ym(i), Mue(i), Mum(i), &
    Mun(i), Mup(i), Un(i), Up(i), EffMassn(i), EffMassp(i)
  END DO
  CLOSE(123)

  CALL CPU_TIME(t_start)

  OpaA_2D = 0.d0
  OpaA_2D_OLD = 0.d0
  OpaA_4D = 0.d0
  OpaA_4D_OLD = 0.d0

  errMax2D = 0.0d0
  errMax4D = 0.0d0
  nSign2D  = 0
  nSign4D  = 0
  nNeg2D   = 0

  nDone = 0
  lastPrint = -1
  total = 4 * 4 * nThermoPoints

  t_2D = 0.0d0
  t_4D = 0.0d0
  t_2D_OLD = 0.0d0
  t_4D_OLD = 0.0d0
  CALL InitOpacity_CC_2D( nE_2D )   ! 2D quadrature nodes, computed once
  DO j = 1, nApprox
    DO k = 1, nOp
      DO i = 1, nThermoPoints
        xTem = T(i)
        if (DoMuons) then
          cheml = Mum(i)
          massl = massm
        else
          cheml = Mue(i)
          massl = masse
        endif
        chemn = Mun(i)
        chemp = Mup(i)
        xUn = Un(i)
        xUp = Up(i)
        xMassn = EffMassn(i)
        xMassp = EffMassp(i)

        CALL CPU_TIME(t1)
        DO l=1, NP
          call Opacity_CC_2D(j-1, k, EnuA(l), OpaA_2D(l, i, j, k), &
                xTem, cheml, chemn, chemp, massl, xMassn, xMassp, xUn, xUp)
        END DO
        CALL CPU_TIME(t2)
        t_2D = t_2D + t2 - t1

        CALL CPU_TIME(t1)
        call Opacity_CC_2D_GSI(j-1, k, NP, EnuA, OpaA_2D_OLD(:, i, j, k), &
                        xTem, cheml, chemn, chemp, massl, xMassn, xMassp, xUn, xUp)
        CALL CPU_TIME(t2)
        ! 1/km to 1/cm for opacities; emissivity (k=4) is already in (s cm^3 MeV)^-1
        IF (k <= 3) OpaA_2D_OLD(:, i, j, k) = OpaA_2D_OLD(:, i, j, k) * km_to_cm_inv
        t_2D_OLD = t_2D_OLD + t2 - t1

        CALL CPU_TIME(t1)
        DO l=1, NP
          CALL Opacity_CC_4D(j-1, k, EnuA(l), OpaA_4D(l, i, j, k), &
                  xTem, cheml, chemn, chemp, massl, xMassn, xMassp, xUn, xUp)
        END DO
        CALL CPU_TIME(t2)
        t_4D = t_4D + t2 - t1

        CALL CPU_TIME(t1)
        call Opacity_CC_4D_GSI(j-1, k, NP, EnuA, OpaA_4D_OLD(:, i, j, k), &
                          xTem, cheml, chemn, chemp, massl, xMassn, xMassp, xUn, xUp)
        CALL CPU_TIME(t2)
        ! 1/km to 1/cm for opacities; emissivity (k=4) is already in (s cm^3 MeV)^-1
        IF (k <= 3) OpaA_4D_OLD(:, i, j, k) = OpaA_4D_OLD(:, i, j, k) * km_to_cm_inv
        t_4D_OLD = t_4D_OLD + t2 - t1

        ! Check 2D
        Error = ABS(OpaA_2D_OLD(:, i, j, k)) - ABS(OpaA_2D(:, i, j, k))
        err = 0.0d0
        DO l=1,NP
          IF (OpaA_2D_OLD(l, i, j, k) /= 0.0d0 .and. OpaA_2D(l, i, j, k) /= 0.0d0) THEN
            err = MAX( err, ABS(Error(l))/ABS(OpaA_2D_OLD(l, i, j, k)) )
          ENDIF
          IF ( SignMismatch(OpaA_2D(l, i, j, k), OpaA_2D_OLD(l, i, j, k)) ) &
            nSign2D(j, k) = nSign2D(j, k) + 1
          IF ( OpaA_2D(l, i, j, k) < 0.0d0 ) nNeg2D(j, k) = nNeg2D(j, k) + 1
        END DO
        errMax2D(j, k) = MAX( errMax2D(j, k), err )

        IF (err > 1.0d-5) THEN
          WRITE(*,*)
          WRITE(*,*) '2D Error large', err
          WRITE(*,*) '2D OLD', OpaA_2D_OLD(:, i, j, k)
          WRITE(*,*) '2D NEW', OpaA_2D(:, i, j, k)
          WRITE(*,*) i, j, k
        ENDIF

        IF ( ANY(OpaA_2D(:, i, j, k) < 0.0d0) ) THEN
          WRITE(*,*)
          WRITE(*,*) '2D NEW negative'
          WRITE(*,*) '2D NEW', OpaA_2D(:, i, j, k)
          WRITE(*,*) i, j, k
        ENDIF

        ! Check 4D
        Error = ABS(OpaA_4D_OLD(:, i, j, k)) - ABS(OpaA_4D(:, i, j, k))
        err = 0.0d0
        DO l=1,NP
          IF (OpaA_4D_OLD(l, i, j, k) /= 0.0d0 .and. OpaA_4D(l, i, j, k) /= 0.0d0) THEN
            err = MAX( err, ABS(Error(l))/ABS(OpaA_4D_OLD(l, i, j, k)) )
          ENDIF
          IF ( SignMismatch(OpaA_4D(l, i, j, k), OpaA_4D_OLD(l, i, j, k)) ) &
            nSign4D(j, k) = nSign4D(j, k) + 1
        END DO
        errMax4D(j, k) = MAX( errMax4D(j, k), err )

      END DO
    END DO
  END DO

  CALL CPU_TIME(t_end)

  WRITE(*,'(/,A,f10.3)') 'Total wall‑clock time :', t_end - t_start
  WRITE(*,'(A,f10.3)')   ' OLD 2D kernels       :', t_2D_old
  WRITE(*,'(A,f10.3)')   ' NEW 2D kernels       :', t_2D
  WRITE(*,'(A,f10.3)')   ' OLD 4D kernels       :', t_4D_old
  WRITE(*,'(A,f10.3)')   ' NEW 4D kernels       :', t_4D

  WRITE(*,'(/,A)') ' Max rel. error on |values| (both nonzero) and sign mismatches'
  WRITE(*,'(A,I8)') ' (sign mismatch = W>0 & G<=0, W<=0 & G>0, or opposite signs); samples per row:', &
    NP * nThermoPoints
  WRITE(*,'(A)') ' Reac Corr |  2D-vs-GSI2D  SignMis  W2D<0 |  4D-vs-GSI4D  SignMis'
  DO k = 1, nOp
    DO j = 1, nApprox
      WRITE(*,'(2I5,A,ES13.4,2I8,A,ES13.4,I8)') k, j-1, ' |', errMax2D(j, k), &
        nSign2D(j, k), nNeg2D(j, k), ' |', errMax4D(j, k), nSign4D(j, k)
    END DO
  END DO
  
  OPEN(200, file='OpaA_2D.bin', status='replace', form='unformatted', access='stream')
  WRITE(200) OpaA_2D
  CLOSE(200)

  OPEN(200, file='OpaA_2D_OLD.bin', status='replace', form='unformatted', access='stream')
  WRITE(200) OpaA_2D_OLD
  CLOSE(200)

  OPEN(201, file='OpaA_4D.bin', status='replace', form='unformatted', access='stream')
  WRITE(201) OpaA_4D
  CLOSE(201)

  OPEN(201, file='OpaA_4D_OLD.bin', status='replace', form='unformatted', access='stream')
  WRITE(201) OpaA_4D_OLD
  CLOSE(201)


  DEALLOCATE(OpaA_2D, OpaA_2D_OLD, OpaA_4D, OpaA_4D_OLD)
  DEALLOCATE(T, Rho, Ye, Ym, Mue, Mum, Mun, Mup, Un, Up, EffMassn, EffMassp)

CONTAINS

  LOGICAL FUNCTION SignMismatch(a, b)
    REAL(DP), INTENT(IN) :: a, b
    SignMismatch = ( a * b < 0.0d0 ) .OR. &
                   ( a > 0.0d0 .AND. b <= 0.0d0 ) .OR. &
                   ( a <= 0.0d0 .AND. b > 0.0d0 )
  END FUNCTION SignMismatch

END PROGRAM wlTestAbEm4DVsGSI


