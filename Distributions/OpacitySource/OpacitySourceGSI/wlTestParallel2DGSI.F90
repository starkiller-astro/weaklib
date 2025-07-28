PROGRAM wlTestParallel2DGSI

  USE omp_lib
  USE wlKindModule, ONLY: dp
  USE wlSemiLeptonicOpacityModule2D, ONLY: &
    Opacity_CC_2D
  USE wlEosConstantsModule, ONLY: &
   pi, Gw_MeV, ga, gv, mn, mp, me, mmu, mpi, &
   Vud, massA, massV, gamma_p, gamma_n

  implicit none

  INTEGER, PARAMETER :: NP = 60, nOp = 4, nApprox = 4, nE_2D = 40
  REAL(DP), PARAMETER :: masse = me , massm = mmu
  REAL(DP), PARAMETER :: massn = mn, massp = mp

  REAL(DP) :: EnuA(NP)
  REAL(DP) :: xTem, cheml, chemn, chemp, xUn, xUp, massl
  INTEGER :: i, j, k, l, nThermoPoints
  LOGICAL, PARAMETER :: DoMuons = .true.

  REAL(DP), allocatable :: OpaA_2D(:,:,:,:)
  REAL(DP), allocatable :: T(:), Rho(:), Ye(:), Ym(:), &
                         Mue(:), Mum(:), Mun(:), Mup(:), Un(:), Up(:)

  REAL(DP) :: OpaA_dummy
  REAL(DP) :: t_tot, t_start, t_end, t_serial
  REAL(DP) :: t_2D, t_2D_serial
  INTEGER :: count_start, count_end1, count_end2, count_rate

  ! Initialize neutrino energies
  DO l = 1, NP
     EnuA(l) = l * 5.d0
  END DO

  ! Read thermodynamic data
  OPEN(unit=123, file=trim(adjustl('ThermoConditions.dat')), status='old', action='read')
  read(123,*) nThermoPoints
  read(123,*)

  ! You can also set nThermoPoints to a smaller value for quick checks
  nThermoPoints = 100
  ALLOCATE(OpaA_2D(NP, nThermoPoints, nApprox, nOp))
  ALLOCATE(T(nThermoPoints), Rho(nThermoPoints), Ye(nThermoPoints), Ym(nThermoPoints))
  ALLOCATE(Mue(nThermoPoints), Mum(nThermoPoints), Mun(nThermoPoints), Mup(nThermoPoints))
  ALLOCATE(Un(nThermoPoints), Up(nThermoPoints))

  DO i = 1, nThermoPoints
     read(123,*) T(i), Rho(i), Ye(i), Ym(i), Mue(i), Mum(i), Mun(i), Mup(i), Un(i), Up(i)
  END DO
  CLOSE(123)

  OpaA_2D = 0.d0

  CALL SYSTEM_CLOCK(count_start, count_rate=count_rate)
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

          CALL Opacity_CC_2D(j-1, k, EnuA(l), OpaA_2D(l, i, j, k), &
                xTem, cheml, chemn, chemp, massl, massn, massp, xUn, xUp, nE_2D)

        END DO
      END DO
    END DO
  END DO
  CALL SYSTEM_CLOCK(count_end1)
  t_2D_serial = REAL(count_end1 - count_start) / REAL(count_rate)

 !$omp parallel do collapse(4) default(shared) &
 !$OMP private(i, j, k, l, xTem, cheml, massl, chemn, chemp, xUn, xUp, OpaA_dummy)
  DO i = 1, nThermoPoints
    DO j = 1, nApprox
      DO k = 1, nOp
        DO l = 1, NP

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

          CALL Opacity_CC_2D(j-1, k, EnuA(l), OpaA_2D(l, i, j, k), &
                xTem, cheml, chemn, chemp, massl, massn, massp, xUn, xUp, nE_2D)

        END DO
      END DO
    END DO
  END DO
 !$omp END PARALLEL DO
  CALL SYSTEM_CLOCK(count_end2)
  t_2D = REAL(count_end2 - count_end1) / REAL(count_rate)
  t_tot = REAL(count_end2 - count_start) / REAL(count_rate)

   WRITE(*,'(/,A,f10.3)') 'Total wall‑clock time :', t_tot
   WRITE(*,'(A,f10.3)')   ' Serial 2D    :', t_2D_serial
   WRITE(*,'(A,f10.3)')   ' Parallel 2D  :', t_2D
   WRITE(*,'(A,i10)')     ' Threads      :', omp_get_max_threads()
   WRITE(*,'(A,f10.3)')   ' Ratio 2D     :', t_2D_serial/t_2D

  !-------------------------------------------------------------
  ! Persist results
  !-------------------------------------------------------------
  OPEN(200, file='OpaA_2D.bin', status='replace', form='unformatted', access='stream')
  WRITE(200) OpaA_2D
  CLOSE(200)

  DEALLOCATE(OpaA_2D)
  DEALLOCATE(T, Rho, Ye, Ym, Mue, Mum, Mun, Mup, Un, Up)

END PROGRAM wlTestParallel2DGSI

