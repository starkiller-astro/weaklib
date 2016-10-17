MODULE GreyVariables
!-----------------------------------------------------------------------
!
!    File:         GreyVariables
!    Module:       GreyVariables
!    Type:         Module
!    Author:       R. Chu, Dept of Physics, UTK
!
!    Date:         8/21/16
!
!    Purpose:
!      Subroutines needed for GreyMoment_Number/Energy
!                        and  GreyOpacity_Number/Energy
!    Mudules used: 
!      wlKindModule
!      gaquad
!-----------------------------------------------------------------------

  USE wlKindModule, ONLY: dp
  USE B85

  PUBLIC :: &
    GreyMomentWithGaussianQuadrature,&
    GreyOpacityWithGaussianQuadrature

CONTAINS

  SUBROUTINE GreyMomentWithGaussianQuadrature&
                  ( nquad, bb, outcome, func, debug )

  INTEGER, INTENT(in)        :: nquad
  REAL(dp), INTENT(in)       :: bb
  CHARACTER(18), INTENT(in)  :: func
  LOGICAL, INTENT(in)        :: debug

  REAL(dp), INTENT(out)      :: outcome

  INTEGER, PARAMETER         :: npiece = 47
  REAL(dp), DIMENSION(nquad) :: roots, weights, FD
  REAL(dp), DIMENSION(npiece):: lim
  INTEGER                    :: ii, jj
  REAL(dp)                   :: llim, ulim

  IF (debug) THEN
    PRINT*,"Calculating ", func
    PRINT*,"with bb = ",bb
  END IF

  outcome = 0.0_dp
  lim(1) = 0.0_dp
 
  DO jj = 1,10
     lim(jj+1) = lim(jj) + 0.1_dp
  END DO
  DO jj = 12, 20
     lim(jj) = lim(jj-1) + 1.0_dp
  END DO
  DO jj = 21, 29
     lim(jj) = lim(jj-1) + 10.0_dp
  END DO
  DO jj = 30, 38
     lim(jj) = lim(jj-1) + 100.0_dp
  END DO
  DO jj = 39, 47
     lim(jj) = lim(jj-1) + 1000.0_dp
  END DO


  DO jj = 1,(npiece-1)

    llim = lim(jj)
    ulim = lim(jj+1)
    
    CALL gaquad( nquad, roots, weights, llim, ulim )
    
    DO ii = 1, nquad
      
      IF (func == "GreyMoment_Number  ") THEN
      
        FD(ii) = Number_FD( roots(ii), bb)
        outcome = outcome + weights(ii) * &
                  Number_FD( roots(ii), bb)

      ELSE IF (func == "GreyMoment_Energy  ") THEN

        FD(ii) = Energy_FD( roots(ii), bb)
        outcome = outcome + weights(ii) * &
                  Energy_FD( roots(ii), bb)

      END IF

    END DO ! ii (nquad)

  END DO ! jj 

  IF (debug) THEN
    WRITE(*,'(A24,F7.1,A2,F7.1,A1)'),&
        " Integration domain is [",lim(1),", ",lim(47),"]"
    print*,"divided into", lim
    PRINT*, "with nquad=", nquad
    PRINT*, "and the outcome is", outcome
  END IF

  END SUBROUTINE GreyMomentWithGaussianQuadrature

  SUBROUTINE GreyOpacityWithGaussianQuadrature&
                  ( nquad, bb, &
                    rho, T, Z, A, chem_e, chem_n,&
                    chem_p, xheavy, xn, xp,&
                    outcome, func, debug )
  INTEGER, INTENT(in)        :: nquad
  REAL(dp), INTENT(in)       :: bb, rho, T, Z, A, &
                                chem_e, chem_n, &
                                chem_p, xheavy, xn,&
                                xp
  CHARACTER(18), INTENT(in)  :: func
  LOGICAL, INTENT(in)        :: debug

  REAL(dp), INTENT(out)      :: outcome

  INTEGER, PARAMETER         :: npiece = 47
  REAL(dp), DIMENSION(nquad) :: roots, weights, opacity
  REAL(dp), DIMENSION(npiece):: lim
  INTEGER                    :: ii, jj
  REAL(dp)                   :: llim, ulim

  IF (debug) THEN
    PRINT*,"Calculating ", func
    PRINT*,"with bb =", bb
  END IF

  outcome = 0.0_dp
  lim(1) = 0.0_dp

  DO jj = 1,10
     lim(jj+1) = lim(jj) + 0.1_dp
  END DO
  DO jj = 12, 20
     lim(jj) = lim(jj-1) + 1.0_dp
  END DO
  DO jj = 21, 29
     lim(jj) = lim(jj-1) + 10.0_dp
  END DO
  DO jj = 30, 38
     lim(jj) = lim(jj-1) + 100.0_dp
  END DO
  DO jj = 39, 47
     lim(jj) = lim(jj-1) + 1000.0_dp
  END DO

  DO jj = 1,(npiece-1)

    llim = lim(jj)
    ulim = lim(jj+1)

    CALL gaquad( nquad, roots, weights, llim, ulim )

    DO ii = 1, nquad

      opacity(ii) = totalECapEm(roots(ii), rho, T, Z, A,&
                      chem_e, chem_n, chem_p, &
                      xheavy, xn, xp )

      IF (func == "GreyOpacity_Number ") THEN

        outcome = outcome + weights(ii) * &
                  Number_FD( roots(ii), bb) * opacity(ii)

      ELSE IF (func == "GreyOpacity_Energy ") THEN

        outcome = outcome + weights(ii) * &
                  Energy_FD( roots(ii), bb) * opacity(ii)
                                                           
      END IF
    END DO ! ii (nquad)

  END DO ! jj 

  IF (debug) THEN
    WRITE(*,'(A24,F7.1,A2,F7.1,A1)'),&
        " Integration domain is [",lim(1),", ",lim(47),"]"
    print*,"divided into", lim
    PRINT*, "with nquad=", nquad
   PRINT*, "and the outcome is", outcome
  END IF

  END SUBROUTINE GreyOpacityWithGaussianQuadrature


  SUBROUTINE GreyOpacityWithGaussianQuadrature_scattIso&
                  ( nquad, bb, &
                    rho, T, xh, A, Z, xn, xp, l,&
                    outcome, func, debug )
  INTEGER, INTENT(in)        :: nquad, l
  REAL(dp), INTENT(in)       :: bb, rho, T, xh, A, Z, xn, xp
  CHARACTER(18), INTENT(in)  :: func
  LOGICAL, INTENT(in)        :: debug

  REAL(dp), INTENT(out)      :: outcome

  INTEGER, PARAMETER         :: npiece = 47
  REAL(dp), DIMENSION(nquad) :: roots, weights, opacity
  REAL(dp), DIMENSION(npiece):: lim
  INTEGER                    :: ii, jj
  REAL(dp)                   :: llim, ulim

  IF (debug) THEN
    PRINT*,"Calculating ", func
    PRINT*,"with bb =", bb
  END IF

  outcome = 0.0_dp
  lim(1) = 0.0_dp

  DO jj = 1,10
     lim(jj+1) = lim(jj) + 0.1_dp
  END DO
  DO jj = 12, 20
     lim(jj) = lim(jj-1) + 1.0_dp
  END DO
  DO jj = 21, 29
     lim(jj) = lim(jj-1) + 10.0_dp
  END DO
  DO jj = 30, 38
     lim(jj) = lim(jj-1) + 100.0_dp
  END DO
  DO jj = 39, 47
     lim(jj) = lim(jj-1) + 1000.0_dp
  END DO

  DO jj = 1,(npiece-1)

    llim = lim(jj)
    ulim = lim(jj+1)

    CALL gaquad( nquad, roots, weights, llim, ulim )

    DO ii = 1, nquad

      opacity(ii) = totalElasticScatteringKernel&
     ( roots(ii), rho, T, xh, A, Z, xn, xp, l )

      IF (func == "GreyOpacity_Number ") THEN

        outcome = outcome + weights(ii) * &
                  Number_FD( roots(ii), bb) * opacity(ii)

      ELSE IF (func == "GreyOpacity_Energy ") THEN

        outcome = outcome + weights(ii) * &
                  Energy_FD( roots(ii), bb) * opacity(ii)
                                                           
      END IF
    END DO ! ii (nquad)

  END DO ! jj 

  IF (debug) THEN
    WRITE(*,'(A24,F7.1,A2,F7.1,A1)'),&
        " Integration domain is [",lim(1),", ",lim(47),"]"
    print*,"divided into", lim
    PRINT*, "with nquad=", nquad
   PRINT*, "and the outcome is", outcome
  END IF

  END SUBROUTINE GreyOpacityWithGaussianQuadrature_scattIso

!-------------------------------------------------------
!    Declear GrayFunctions
!-------------------------------------------------------

  FUNCTION Number_FD( x, b )

  USE wlKindModule, ONLY: dp
   REAL(dp), INTENT(in)  :: x,b    ! function argument
   REAL(dp)              :: fexp, Number_FD

   Number_FD =  ( x*x )/( fexp(x-b) + 1.0 )

  END FUNCTION Number_FD

  FUNCTION Energy_FD( x, b )

  USE wlKindModule, ONLY: dp

   REAL(dp), INTENT(in)  :: x,b    ! function argument
   REAL(dp)              :: fexp, Energy_FD

   Energy_FD =  ( x*x*x )/( fexp(x-b) + 1.0 )
                                                     
  END FUNCTION Energy_FD

END MODULE GreyVariables
