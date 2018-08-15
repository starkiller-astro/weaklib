MODULE B85_scattTP

  USE wlKindModule, ONLY: dp
  USE wlExtPhysicalConstantsModule, ONLY: &
    h, kMeV, therm1, therm2, dmnp, me, mbG, mp, mn, cvel_inv, cvel, ergmev,&
      cv_p, cv_n, ca_p, ca_n, gf, hbarc, cv, ca
  USE wlExtNumericalModule, ONLY: pi, half, twpi, zero, one

  IMPLICIT NONE

  PUBLIC TPKernel

CONTAINS

  SUBROUTINE TPKernel( energygrid, TMeV, chem_e, nquad, l, TPK, species )
!------------------------------------------------------------------------------
! Purpose:
!   To compute the zero and first legendre coefs for the thermal production
!   annihilation kernel for electron-type neutrino and antineutrino 
!------------------------------------------------------------------------------
  IMPLICIT NONE

  REAL(dp), DIMENSION(:), INTENT(in)    :: energygrid
  REAL(dp), INTENT(in)                  :: TMeV, chem_e
  INTEGER , INTENT(in)                  :: nquad, l, species
  REAL(dp), DIMENSION(:,:), INTENT(out) :: TPK

  REAL(dp), DIMENSION(:), ALLOCATABLE :: TPK_Ee
  INTEGER                    :: ii, jj, kk, zz, nPointsE
  INTEGER, PARAMETER         :: npiece = 47
  REAL(dp), DIMENSION(nquad) :: roots, weights, midFD
  REAL(dp), DIMENSION(npiece):: lim
  REAL(dp)                   :: UnitConvertConstant, outcome, llim, ulim

  UnitConvertConstant = cvel_inv**4.0 / h**3.0 !! double check 
  nPointsE = SIZE(energygrid)

!  ALLOCATE( TPK_Ee( npiece ) )
!
!  lim(1) = zero
!
!  DO jj = 1,10
!     lim(jj+1) = lim(jj) + 0.1_dp
!  END DO
!  DO jj = 12, 20
!     lim(jj) = lim(jj-1) + 1.0_dp
!  END DO
!  DO jj = 21, 29
!     lim(jj) = lim(jj-1) + 10.0_dp
!  END DO
!  DO jj = 30, 38
!     lim(jj) = lim(jj-1) + 100.0_dp
!  END DO
!  DO jj = 39, 47
!     lim(jj) = lim(jj-1) + 1000.0_dp !!!!!!!!!!!!!!!!!!!!!!!
!  END DO
!
  DO jj = 1, nPointsE  ! ep
    DO kk = 1, nPointsE ! e

       outcome = zero

!       DO ii = 1, (npiece-1)
!
!         llim = lim(ii)
!         ulim = lim(ii+1)
!
!         CALL gaquad( nquad, roots, weights, llim, ulim )
!         CALL TPAbsorptionKernelWithElectronEnergy&
!              ( energygrid(kk), energygrid(jj), roots, TMeV, chem_e, l, TPK_Ee )
!         CALL TPAK_midFD &
!              ( energygrid(kk), energygrid(jj), roots, TMeV, chem_e, midFD )
!  
!         DO zz = 1, nquad
!
!           outcome = outcome + weights(zz) * TPK_Ee(zz) * midFD(zz)
!
!         END DO
!
!       END DO

      TPK( jj, kk ) = UnitConvertConstant !* outcome 

    END DO
  END DO

  END SUBROUTINE TPKernel


  SUBROUTINE TPAbsorptionKernelWithElectronEnergy &
             ( e, ep, Ee, TMeV, chem_e, l, TPKarr )

  IMPLICIT NONE

    REAL(dp), DIMENSION(:), INTENT(in) :: Ee
    REAL(dp), INTENT(in)               :: e, ep, TMeV, chem_e
    INTEGER,  INTENT(in)               :: l
    REAL(dp), DIMENSION(:), INTENT(out) :: TPKarr

    REAL(dp)                   :: Coeff1, Coeff2
    INTEGER                    :: nEe
    REAL(dp), DIMENSION(:), ALLOCATABLE :: JI, JII

    nEe    = SIZE(Ee)
    Coeff1 = ( cv + ca ) * ( cv + ca )
    Coeff2 = ( cv - ca ) * ( cv - ca )

    ALLOCATE( JI(nEe), JII(nEe) )

!    CALL ComputeJfuncForTPK( e, ep, Ee, JI )
!    CALL ComputeJfuncForTPK( ep, e, Ee, JII )

    TPKarr = Coeff1 * JI + Coeff2 * JII

    DEALLOCATE( JI, JII)

  END SUBROUTINE TPAbsorptionKernelWithElectronEnergy


  SUBROUTINE ComputeJfuncForTPK( ep, e, Ee, TMeV, chem_e, J0I, J0II, J1I, J1II )

  IMPLICIT NONE
 
    REAL(dp), INTENT(in)  :: ep, e, TMeV, chem_e
    REAL(dp), INTENT(in)  :: Ee
    !REAL(dp), DIMENSION(:), INTENT(in) :: Ee
    REAL(dp), INTENT(out) :: J0I, J0II, J1I, J1II
   ! REAL(dp), DIMENSION(:), INTENT(out) :: J0I, J0II, J1I, J1II

    REAL(dp)  :: w, w2, w3, w4, w5, wp, wp2, wp3, wp4, wp5, wwp, &
                 w_p_wp, w_p_wp2, w_p_wp3, beta, eta
    REAL(dp)  :: coef
    
      coef     = 1.0_dp !!!!
      w        = e
      w2       = w * w
      w3       = w2 * w
      w4       = w3 * w
      w5       = w4 * w
      wp       = ep
      wp2      = wp * wp
      wp3      = wp2 * wp
      wp4      = wp3 * wp
      wp5      = wp4 * wp
      wwp      = w * wp
      w_p_wp   = w + wp
      w_p_wp2  = w_p_wp * w_p_wp
      w_p_wp3  = w_p_wp2 * w_p_wp
      beta     = one/TMeV
      eta      = beta * chem_e

      J0I      = zero
      J0II     = zero
      J1I      = zero
      J1II     = zero
     
!-----------------------------------------------------------------------
!  w > wp
!-----------------------------------------------------------------------

      IF ( w >= wp ) THEN
      
        CALL pr_w_gt_wp_e_lt_wp( zero, wp, eta, beta, j0i, j0ii, j1i, j1ii )
        CALL pr_w_gt_wp_lt_e_lt_w( wp, w, eta, beta, j0i, j0ii, j1i, j1ii )
        CALL pr_w_gt_wp_w_lt_e( w, w_p_wp, eta, beta, j0i, j0ii, j1i, j1ii )
        
      ELSE
      
      !-----------------------------------------------------------------------
      !  w < wp
      !-----------------------------------------------------------------------
      
        CALL pr_w_lt_wp_e_lt_w( zero, w, eta, beta, j0i, j0ii, j1i, j1ii )
        CALL pr_w_lt_wp_e_lt_wp( w, wp, eta, beta, j0i, j0ii, j1i, j1ii )
        CALL pr_w_lt_wp_wp_lt_e( wp, w_p_wp, eta, beta, j0i, j0ii, j1i, j1ii )
        
      END IF
      
      !-----------------------------------------------------------------------
      !  .Final touches
      !-----------------------------------------------------------------------
      
      j0i                = j0i /( e * e * ep * ep )
      j0ii               = j0ii/( e * e * ep * ep )
      j1i                = j1i /( e * e * ep * ep )
      j1ii               = j1ii/( e * e * ep * ep )
      
      j0i                = coef * j0i
      j0ii               = coef * j0ii
      j1i                = coef * j1i
      j1ii               = coef * j1ii
      
      RETURN

  END SUBROUTINE ComputeJfuncForTPK


  SUBROUTINE TPAK_midFD( e1, e2, ee, TMeV, chem_e, midFD ) 

  IMPLICIT NONE

  REAL(dp), INTENT(in)                :: e1, e2, TMeV, chem_e
  REAL(dp), DIMENSION(:), INTENT(in)  :: ee
  REAL(dp), DIMENSION(:), INTENT(out) :: midFD

  REAL(dp) :: fexp

  midFD = FEXP( (e1+e2)/TMeV ) &
         / ( (FEXP( (ee-chem_e) / TMeV ) + 1.0_dp) &
           * (FEXP( (e1+e2-ee+chem_e)/TMeV ) + 1.0_dp) )

  END SUBROUTINE TPAK_midFD  


  SUBROUTINE pr_w_gt_wp_e_lt_wp( xl, xu, eta, beta, j0i, j0ii, j1i, j1ii )

  IMPLICIT NONE

    REAL(dp), INTENT(in)  :: xl, xu, eta, beta
    REAL(dp), INTENT(out) :: j0i, j0ii, j1i, j1ii

    INTEGER  :: i
    REAL(dp) :: h0i, h0ii, su0i, su0ii, su1i, su1ii, e_mid, e_del,&
                e_var, arg1, arg2, exp1, exp2, eblock, FEXP

    REAL(dp) :: r4_15, r4_3, r8_3

    r4_15 = 4.0_dp / 15.0_dp
    r4_3  = 4.0_dp / 3.0_dp
    r8_3  = 8.0_dp / 3.0_dp

      IF ( xl == xu ) RETURN
      
      !--------------------------------------------------------------
      !  Initialize for integration
      !--------------------------------------------------------------
      
      su0i               = zero
      su0ii              = zero
      su1i               = zero
      su1ii              = zero
      
      e_mid              = half * ( xu + xl )
      e_del              = half * ( xu - xl )
      
!      !-----------------------------------------------------------------------
!      !  Integrate
!      !-----------------------------------------------------------------------
!      
!      DO i = 1,nleg
!      
!        e_var            = xe(i) * e_del
!        e                = e_mid + e_var
!      
!        arg1             = beta * e - eta
!        exp1             = fexp(arg1)
!        arg2             = beta * ( w_p_wp - e ) + eta
!        exp2             = fexp(arg2)
!        eblock           = ( exp1/( exp1 + one ) ) * ( exp2/( exp2 + one ) )
!      
!        h0i              = r4_15 * e**5 - r4_3 * e**4 * wp + r8_3 * e**3 * wp2
!        h0ii             = r4_15 * e**5 - r4_3 * e**4 * w + r8_3 * e**3 * w2
!      
!        su0i             = su0i  + eblock * h0i  * wte(i)
!        su0ii            = su0ii + eblock * h0ii * wte(i)
!      
!      END DO
!      
      su0i               = e_del * su0i
      su0ii              = e_del * su0ii
      su1i               = e_del * su1i
      su1ii              = e_del * su1ii
      
      
      j0i                = j0i  + su0i
      j0ii               = j0ii + su0ii
      j1i                = j1i  + su1i
      j1ii               = j1ii + su1ii
      
      RETURN
      
  END SUBROUTINE pr_w_gt_wp_e_lt_wp


  SUBROUTINE pr_w_gt_wp_lt_e_lt_w( xl, xu, eta, beta, j0i, j0ii, j1i, j1ii )

  IMPLICIT NONE

    REAL(dp), INTENT(in)  :: xl, xu, eta, beta
    REAL(dp), INTENT(out) :: j0i, j0ii, j1i, j1ii

  END SUBROUTINE pr_w_gt_wp_lt_e_lt_w


  SUBROUTINE pr_w_gt_wp_w_lt_e( xl, xu, eta, beta, j0i, j0ii, j1i, j1ii )

  IMPLICIT NONE

    REAL(dp), INTENT(in)  :: xl, xu, eta, beta
    REAL(dp), INTENT(out) :: j0i, j0ii, j1i, j1ii

  END SUBROUTINE pr_w_gt_wp_w_lt_e


  SUBROUTINE pr_w_lt_wp_e_lt_w( xl, xu, eta, beta, j0i, j0ii, j1i, j1ii )

  IMPLICIT NONE

    REAL(dp), INTENT(in)  :: xl, xu, eta, beta
    REAL(dp), INTENT(out) :: j0i, j0ii, j1i, j1ii

  END SUBROUTINE pr_w_lt_wp_e_lt_w


  SUBROUTINE pr_w_lt_wp_e_lt_wp( xl, xu, eta, beta, j0i, j0ii, j1i, j1ii )

  IMPLICIT NONE

    REAL(dp), INTENT(in)  :: xl, xu, eta, beta
    REAL(dp), INTENT(out) :: j0i, j0ii, j1i, j1ii

  END SUBROUTINE pr_w_lt_wp_e_lt_wp


  SUBROUTINE pr_w_lt_wp_wp_lt_e( xl, xu, eta, beta, j0i, j0ii, j1i, j1ii )

  IMPLICIT NONE

    REAL(dp), INTENT(in)  :: xl, xu, eta, beta
    REAL(dp), INTENT(out) :: j0i, j0ii, j1i, j1ii

  END SUBROUTINE pr_w_lt_wp_wp_lt_e
  
END MODULE B85_scattTP
