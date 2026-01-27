PROGRAM wlTestGSI

  USE wlKindModule, ONLY: dp
  USE wlSemiLeptonicOpacityModule2D, ONLY: &
    Opacity_CC_2D
  USE wlEosConstantsModule, ONLY: &
   pi, Gw_MeV, ga, gv, mn, mp, me, mmu, mpi, &
   Vud, massA, massV, gamma_p, gamma_n, kmev

  IMPLICIT NONE

  INTEGER, PARAMETER :: NP = 60, nOp = 4, nApprox = 4, nE_2D = 50
  REAL(DP), PARAMETER :: masse = me , massm = mmu

  REAL(DP) :: EnuA(NP), Error(NP)
  REAL(DP) :: xTem, cheml, chemn, chemp, xUn, xUp, massl, xMassn, xMassp
  INTEGER :: i, j, k, l, nThermoPoints
  INTEGER :: nDone, total, lastPrint
  LOGICAL, PARAMETER :: DoMuons = .true.

  REAL(DP), ALLOCATABLE :: OpaA_2D(:,:,:,:)
  REAL(DP), ALLOCATABLE :: T(:), Rho(:), Ye(:), Ym(:), Mue(:), Mum(:), &
      Mun(:), Mup(:), Un(:), Up(:), EffMassn(:), EffMassp(:)

  REAL(DP) :: t1, t2, t_start, t_end, t_2D

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
  ALLOCATE(T(nThermoPoints), Rho(nThermoPoints), Ye(nThermoPoints), Ym(nThermoPoints))
  ALLOCATE(Mue(nThermoPoints), Mum(nThermoPoints), Mun(nThermoPoints), Mup(nThermoPoints))
  ALLOCATE(Un(nThermoPoints), Up(nThermoPoints), EffMassn(nThermoPoints), EffMassp(nThermoPoints))

  DO i = 1, nThermoPoints
    READ(123,*) T(i), Rho(i), Ye(i), Ym(i), &
        Mue(i), Mum(i), Mun(i), Mup(i), Un(i), Up(i), EffMassn(i), EffMassp(i)
  END DO
  CLOSE(123)

  ! ! Force single point to compare to Table II of Fischer 2020
  ! Rho(1) = 1.0d13
  ! T(1)  = 10.0d0
  ! Ye(1) = 0.2d0 !Does not matter...
  ! Ym(1) = 1.0d-4
  ! Mue(1) = 59.8987d0
  ! Mum(1) = 35.48134000d0
  ! Mun(1) = -15.50788137687d0 + mn
  ! Mup(1) = -38.7855487990624d0 + mp
  ! Un(1)  = 5.6396315248d0
  ! Up(1)  = 3.58702109103d0
  ! EffMassn(1)  = 930.881301628124999999999999999999980d0
  ! EffMassp(1)  = 929.588001628125000000000000000000005d0

  CALL CPU_TIME(t_start)

  OpaA_2D = 0.d0

  nDone = 0
  lastPrint = -1
  total = 4 * 4 * nThermoPoints

  t_2D = 0.0d0
  DO j = 4, nApprox
    DO k = 1, 2
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
        xmassn = EffMassn(i)
        xmassp = EffMassp(i)

        CALL CPU_TIME(t1)
        DO l=1, NP
          call Opacity_CC_2D(j-1, k, EnuA(l), OpaA_2D(l, i, j, k), &
                xTem, cheml, chemn, chemp, massl, xmassn, xmassp, xUn, xUp, nE_2D)
        END DO
        
        CALL CPU_TIME(t2)
        t_2D = t_2D + t2 - t1

        nDone = nDone + 1
        if ((nDone * 100 / total) > lastPrint) then
          lastPrint = nDone * 100 / total
          print*, 'Progress: ', lastPrint, '% completed'
        endif

      END DO
    END DO
  END DO

  CALL CPU_TIME(t_end)

  WRITE(*,'(/,A,f10.3)') 'Total wall‑clock time :', t_end - t_start
  WRITE(*,'(A,f10.3)')   ' NEW 2D kernels       :', t_2D

  OPEN(200, file='OpaA_2D.bin', status='replace', form='unformatted', access='stream')
  WRITE(200) OpaA_2D
  CLOSE(200)

  DEALLOCATE(OpaA_2D)
  DEALLOCATE(T, Rho, Ye, Ym, Mue, Mum, Mun, Mup, Un, Up, EffMassn, EffMassp)

END PROGRAM wlTestGSI