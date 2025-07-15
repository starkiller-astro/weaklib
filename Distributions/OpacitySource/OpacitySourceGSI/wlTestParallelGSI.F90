PROGRAM wlTestParallelGSI

  USE omp_lib
  USE wlKindModule, ONLY: dp
  USE wlSemiLeptonicOpacityModule2D, ONLY: &
    Opacity_CC_2D
  USE wlSemiLeptonicOpacityModule4D, ONLY: &
    Opacity_CC_4D
  USE wlEosConstantsModule, ONLY: &
   pi, Gw_MeV, ga, gv, mn, mp, me, mmu, mpi, &
   Vud, massA, massV, gamma_p, gamma_n

  implicit none

  INTEGER, PARAMETER :: NP = 60, nOp = 4, nApprox = 4
  REAL(DP), PARAMETER :: masse = me , massm = mmu
  REAL(DP), PARAMETER :: massn = mn, massp = mp

  REAL(DP) :: EnuA(NP)
  REAL(DP) :: xTem, cheml, chemn, chemp, xUn, xUp, massl
  INTEGER :: i, j, k, l, nThermoPoints
  LOGICAL, PARAMETER :: DoMuons = .true.

  REAL(DP), allocatable :: OpaA_2D(:,:,:,:), OpaA_4D(:,:,:,:)
  REAL(DP), allocatable :: OpaA_2D_OLD(:,:,:,:), OpaA_4D_OLD(:,:,:,:)
  REAL(DP), allocatable :: T(:), Rho(:), Ye(:), Ym(:), &
                         Mue(:), Mum(:), Mun(:), Mup(:), Un(:), Up(:)

  REAL(DP) :: t_tot, t_start, t_end, t_serial
  REAL(DP) :: t1, t2
  REAL(DP) :: t_2D, t_4D
  REAL(DP) :: t_2D_serial, t_4D_serial

  ! Initialize neutrino energies
  DO l = 1, NP
     EnuA(l) = l * 5.d0
  END DO

  ! Read thermodynamic data
  OPEN(unit=123, file=trim(adjustl('ThermoConditions.dat')), status='old', action='read')
  read(123,*) nThermoPoints
  read(123,*)

  ! You can also set nThermoPoints to a smaller value for quick checks
  ! nThermoPoints = 100
  ALLOCATE(OpaA_2D(NP, nThermoPoints, nApprox, nOp))
  ALLOCATE(OpaA_4D(NP, nThermoPoints, nApprox, nOp))
  ALLOCATE(T(nThermoPoints), Rho(nThermoPoints), Ye(nThermoPoints), Ym(nThermoPoints))
  ALLOCATE(Mue(nThermoPoints), Mum(nThermoPoints), Mun(nThermoPoints), Mup(nThermoPoints))
  ALLOCATE(Un(nThermoPoints), Up(nThermoPoints))

  DO i = 1, nThermoPoints
     read(123,*) T(i), Rho(i), Ye(i), Ym(i), Mue(i), Mum(i), Mun(i), Mup(i), Un(i), Up(i)
  END DO
  CLOSE(123)

  CALL CPU_TIME(t_start)

  OpaA_2D = 0.d0
  OpaA_4D = 0.d0

  t_2D_serial = 0.0d0
  t_4D_serial = 0.0d0
  DO j = 1, nApprox
    DO k = 1, nOp
      DO i = 1, nThermoPoints
        DO l=1, NP

          ! Thread‑local thermodynamic state copy
          xTem  = T(i)
          if (DoMuons) then
              cheml = Mum(i)
              massl = massm
          else
              cheml = Mue(i)
              massl = masse
          END if
          chemn = Mun(i)
          chemp = Mup(i)
          xUn   = Un(i)
          xUp   = Up(i)

          CALL cpu_time(t1)
          CALL Opacity_CC_2D(j-1, k, EnuA(l), OpaA_2D(l, i, j, k), &
                xTem, cheml, chemn, chemp, massl, massn, massp, xUn, xUp)
          CALL cpu_time(t2)
          t_2D_serial = t_2D_serial + (t2 - t1)

          !---------------- 4‑D integral ----------------
          CALL cpu_time(t1)
          CALL Opacity_CC_4D(j-1, k, EnuA(l), OpaA_4D(l, i, j, k), &
                xTem, cheml, chemn, chemp, massl, massn, massp, xUn, xUp)
          CALL cpu_time(t2)
          t_4D_serial = t_4D_serial + (t2 - t1)

        END DO
      END DO
    END DO
  END DO

  CALL CPU_TIME(t_serial)

  t_2D = 0.0d0
  t_4D = 0.0d0
  !$omp parallel private(i, j, k, l, xTem, cheml, massl, chemn, chemp, xUn, xUp, t1, t2) &
  !$omp reduction(+:t_2D, t_4D)
  !$omp DO schedule(dynamic,1) collapse(4)
  DO j = 1, nApprox
    DO k = 1, nOp
      DO i = 1, nThermoPoints
        DO l=1, NP

          ! Thread‑local thermodynamic state copy
          xTem  = T(i)
          if (DoMuons) then
              cheml = Mum(i)
              massl = massm
          else
              cheml = Mue(i)
              massl = masse
          END if
          chemn = Mun(i)
          chemp = Mup(i)
          xUn   = Un(i)
          xUp   = Up(i)

          t1 = omp_get_wtime()
          CALL Opacity_CC_2D(j-1, k, EnuA(l), OpaA_2D(l, i, j, k), &
                xTem, cheml, chemn, chemp, massl, massn, massp, xUn, xUp)
          t2 = omp_get_wtime()
          t_2D = t_2D + (t2 - t1)

          !---------------- 4‑D integral ----------------
          t1 = omp_get_wtime()
          CALL Opacity_CC_4D(j-1, k, EnuA(l), OpaA_4D(l, i, j, k), &
                xTem, cheml, chemn, chemp, massl, massn, massp, xUn, xUp)
          t2 = omp_get_wtime()
          t_4D = t_4D + (t2 - t1)

        END DO
      END DO
    END DO
  END DO
  !$omp END DO
  !$omp END parallel

  CALL CPU_TIME(t_end)
  t_tot = t_end - t_start

   WRITE(*,'(/,A,f10.3)') 'Total wall‑clock time :', t_tot
   WRITE(*,'(A,f10.3)')   ' Serial 2D    :', t_2D_serial
   WRITE(*,'(A,f10.3)')   ' Parallel 2D  :', t_2D
   WRITE(*,'(A,f10.3)')   ' Serial 4D    :', t_4D_serial
   WRITE(*,'(A,f10.3)')   ' Parallel 4D  :', t_4D
   WRITE(*,'(A,i10)')     ' Threads      :', omp_get_max_threads()
   WRITE(*,'(A,f10.3)')   ' Ratio 2D     :', t_2D/t_2D_serial
   WRITE(*,'(A,f10.3)')   ' Ratio 4D     :', t_4D/t_4D_serial
   WRITE(*,'(A,f10.3)')   ' Ratio Serial/Parallel :', &
      (t_serial - t_start)/(t_end - t_serial)

  !-------------------------------------------------------------
  ! Persist results
  !-------------------------------------------------------------
  OPEN(200, file='OpaA_2D.bin', status='replace', form='unformatted', access='stream')
  WRITE(200) OpaA_2D
  CLOSE(200)

  OPEN(201, file='OpaA_4D.bin', status='replace', form='unformatted', access='stream')
  WRITE(201) OpaA_4D
  CLOSE(201)

  DEALLOCATE(OpaA_2D, OpaA_4D)
  DEALLOCATE(T, Rho, Ye, Ym, Mue, Mum, Mun, Mup, Un, Up)

END PROGRAM wlTestParallelGSI

