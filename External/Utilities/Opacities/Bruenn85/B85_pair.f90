MODULE B85_pair
!--------------------------------------------------------------------
!
!    File:         B85_pair.f90
!    Module:       B85_pair
!    Type:         Module w/ Routines
!    Author:       R. Chu, Dept. Phys. & Astronomy
!                  U. Tennesee, Knoxville
!
!    Created:      8/16/18
!
!    WeakLib ver:  weaklib/External/Utilities/Opacities/Bruenn85
!
!    Purpose:      The zero and first legendre coefs for the n-type 
!                  neutrino pair annihilation kernel are computed 
!                  considering the physics in Bruenn 85.
!
!    CONTAINS:     Routine  TPKernel:
!                         gives thermal production kernel
!
!    Input arguments:
!                  energygrid : neutrino energy array [MeV]
!                  TMeV       : matter temperature [MeV]
!                  chem_e     : electron chemical potential 
!                                      (including rest mass) [MeV]
!                  nquad      : number of quadarture       
!                  l          : order of the legendre coefs/moment
!                  species    : neutrino flavor index    
!
!    Output arguments:
!                  TPK    
!
!    Modules used:
!                  wlKindModule
!                  wlExtPhysicalConstantsModule
!                  wlExtNumericalModule
!                  ( function fexp is called )
!
!--------------------------------------------------------------------
  USE wlKindModule, ONLY: dp
  USE wlExtPhysicalConstantsModule, ONLY: &
        h, kMeV, therm1, therm2, therm3, dmnp, me, mbG, mp, mn, cvel_inv, &
        cvel, ergmev, cv_p, cv_n, ca_p, ca_n, gf, hbarc, cv, ca, &
        Gw, hbar
  USE wlExtNumericalModule, ONLY: &
        pi, half, twpi, zero, one

  IMPLICIT NONE

  PUBLIC TPKernel

CONTAINS

  SUBROUTINE TPKernel &
               ( energygrid, TMeV, chem_e, nquad, l, TPK, species )
!---------------------------------------------------------------------
! Purpose:
!     To compute the zero and first legendre coefs for the thermal 
!     production annihilation kernel for electron-type neutrino and
!     antineutrino 
!
! Output:
!     TPK(ep,e) for a given TMeV, chem_e and moment order l
!---------------------------------------------------------------------
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
  REAL(dp)                   :: outcome, ep, e

  nPointsE = SIZE(energygrid)

  TPK = zero

  DO jj = 1, nPointsE  ! ep
    
    ep = energygrid(jj)

    DO kk = 1, nPointsE ! e
    
      e = energygrid(kk)

      CALL pair_cal( e, ep, TMeV, chem_e, l, outcome )

      TPK( jj, kk ) = therm3 * outcome 

    END DO
  END DO

  END SUBROUTINE TPKernel


  SUBROUTINE pair_cal( e, ep, TMeV, chem_e, l, outcome )
!--------------------------------------------------------------------
! Purpose:
!   To integrates the quantities
!              [1-Fe(Ee)]*[1-Febar(e + ep - Ee)]*Phil(e,ep,Ee)
!   from 0 to e + ep over Ee.
!   Ee     : electron energy
!   Phil   : compute by pair_Phil
!
! Input:
!   e      : neutrino energy [MeV]
!   ep     : antineutrino energy [MeV]
!   TMeV   : matter temperature [MeV]
!   chem_e : electron chemical potential [MeV]
!   l      : moment order index
!
! Output:
!  outcome : the value of pair annihilation kernal
!--------------------------------------------------------------------
  IMPLICIT NONE
    INTEGER,  INTENT(in)       :: l
    REAL(dp), INTENT(in)       :: e, ep, TMeV, chem_e
    REAL(dp), INTENT(out)      :: outcome

    INTEGER, PARAMETER         :: nquad = 30 ! same with Grey nquad
    INTEGER                    :: ii
    REAL(dp)                   :: phi, midFe, FEXP
    REAL(dp), DIMENSION(nquad) :: roots_Ee, weights

    CALL gaquad( nquad, roots_Ee, weights, 0.0_dp, e+ep )
    
    outcome = zero

    DO ii = 1, nquad

      CALL pair_midFe( e, ep, roots_Ee(ii), TMeV, chem_e, midFe )
      CALL pair_Phil(  e, ep, roots_Ee(ii), l, phi )

      outcome = outcome + weights(ii) * midFe * phi

    END DO


  END SUBROUTINE pair_cal


  SUBROUTINE pair_Phil &
             ( e, ep, Ee, l, outcome )
!--------------------------------------------------------------------
! Purpose:
!   To compute Phil(e,ep,Ee) needed by paircal
!--------------------------------------------------------------------
  IMPLICIT NONE

    REAL(dp), INTENT(in)  :: e, ep, Ee
    INTEGER,  INTENT(in)  :: l
    REAL(dp), INTENT(out) :: outcome

    REAL(dp)              :: Coeff1, Coeff2, JI, JII
    INTEGER               :: nEe

    Coeff1 = ( cv + ca ) * ( cv + ca )
    Coeff2 = ( cv - ca ) * ( cv - ca )

    CALL pair_J( e, ep, Ee, l, JI )
    CALL pair_J( ep, e, Ee, l, JII )
    
    outcome = Coeff1 * JI + Coeff2 * JII

  END SUBROUTINE pair_Phil


  SUBROUTINE pair_J( e, ep, Ee, l, J )

  IMPLICIT NONE
 
    INTEGER,  INTENT(in)  :: l
    REAL(dp), INTENT(in)  :: e, ep, Ee
    REAL(dp), INTENT(out) :: J

      IF ( e >= ep ) THEN
      
        CALL pr_w_gt_wp_wp_gt_e( e, ep, Ee, l, J )
        CALL pr_w_gt_e_e_gt_wp ( e, ep, Ee, l, J )
        CALL pr_w_gt_wp_w_lt_e ( e, ep, Ee, l, J )
        
      ELSE
      
        CALL pr_w_lt_wp_e_lt_w ( e, ep, Ee, l, J )
        CALL pr_w_lt_e_e_lt_wp ( e, ep, Ee, l, J )
        CALL pr_w_lt_wp_wp_lt_e( e, ep, Ee, l, J )
        
      END IF

      J = J / (e*ep)

  END SUBROUTINE pair_J


  SUBROUTINE pair_midFe( e, ep, Ee, TMeV, chem_e, outcome ) 

  IMPLICIT NONE
  REAL(dp), INTENT(in)  :: e, ep, Ee, TMeV, chem_e
  REAL(dp), INTENT(out) :: outcome

  REAL(dp) :: FEXP

  outcome = FEXP( (e+ep)/TMeV ) &
            / ( (FEXP( (Ee-chem_e)      / TMeV ) + 1.0_dp) &
              * (FEXP( (e+ep-Ee+chem_e) / TMeV ) + 1.0_dp) )

  END SUBROUTINE pair_midFe  


  SUBROUTINE pr_w_gt_wp_wp_gt_e( e, ep, Ee, l, J_a )
!-----------------------------------------------
! To compute a0( e, ep, Ee ) and a1( e, ep, Ee )
!-----------------------------------------------
  IMPLICIT NONE

    INTEGER,  INTENT(in)    :: l
    REAL(dp), INTENT(in)    :: e, ep, Ee
    REAL(dp), INTENT(inout) :: J_a

    REAL(dp) :: a0, a1
    REAL(dp) :: r4_3, r4_5, r4_15, r8_3, r16_35
    REAL(dp) :: eep, ep2, eep3, eep_2, Ee3, Ee4, Ee5


    eep    =  e   * ep
    ep2    =  ep  * ep
    eep3   =  eep * ep * ep
    eep_2  =  eep * eep
    Ee3    =  Ee  * Ee * Ee
    Ee4    =  Ee3 * Ee
    Ee5    =  Ee4 * Ee 

    r4_3   =  4.0_dp  / 3.0_dp
    r4_5   =  4.0_dp  / 5.0_dp
    r4_15  =  4.0_dp  / 15.0_dp
    r8_3   =  8.0_dp  / 3.0_dp
    r16_35 =  16.0_dp / 35.0_dp

    IF ( l == 0 ) THEN

      a0   =  (r4_15*Ee5 - r4_3*Ee4*ep + r8_3*Ee3*ep2) / eep 
      J_a  =   J_a + a0

    ELSE IF ( l == 1 ) THEN

      a1   =  ( r16_35*Ee3*Ee4 - r4_5*Ee3*Ee3*(e+3.0_dp*ep) &
               + r4_15*Ee5*ep*(13.0_dp*e+18.0_dp*ep) &
               - r4_3*Ee4*ep2*(4.0_dp*e+3.0_dp*ep) &
               + r8_3*Ee3*eep3) / eep_2 
      J_a  =  J_a + a1 

      END IF


  END SUBROUTINE pr_w_gt_wp_wp_gt_e


  SUBROUTINE pr_w_gt_wp_w_lt_e( e, ep, Ee, l, J_b )
!-----------------------------------------------
! To compute b0( e, ep, Ee ) and b1( e, ep, Ee )
!-----------------------------------------------

  IMPLICIT NONE

    INTEGER,  INTENT(in)    :: l
    REAL(dp), INTENT(in)    :: e, ep, Ee
    REAL(dp), INTENT(inout) :: J_b

    REAL(dp)  :: b0, b1, a0, a1
    REAL(dp)  :: r4_3, r4_5, r4_15, r8_3, r8_5, r8_7, r8_15, &
                 r12_5, r12_35, r16_5, r16_35, r28_5, r28_15
    REAL(dp)  :: e2, e3, e4, e5, ep2, ep3, ep4, ep5, Ee2, Ee3, &
                 Ee4, Ee5, eep, eep3, eep_2, eaep, eaep2, eaep3

    r4_3   = 4.0_dp  / 3.0_dp
    r4_5   = 4.0_dp  / 5.0_dp
    r4_15  = 4.0_dp  / 15.0_dp
    r8_3   = 8.0_dp  / 3.0_dp
    r8_5   = 8.0_dp  / 5.0_dp
    r8_7   = 8.0_dp  / 7.0_dp
    r8_15  = 8.0_dp  / 15.0_dp
    r12_5  = 12.0_dp / 5.0_dp
    r12_35 = 12.0_dp / 35.0_dp
    r16_5  = 16.0_dp / 5.0_dp
    r16_35 = 16.0_dp / 35.0_dp
    r28_5  = 28.0_dp / 5.0_dp
    r28_15 = 28.0_dp / 15.0_dp

    e2     =  e     * e
    e3     =  e2    * e
    e4     =  e3    * e
    e5     =  e4    * e
    ep2    =  ep    * ep
    ep3    =  ep2   * ep
    ep4    =  ep3   * ep
    ep5    =  ep4   * ep
    Ee2    =  Ee    * Ee
    Ee3    =  Ee2   * Ee
    Ee4    =  Ee3   * Ee
    Ee5    =  Ee4   * Ee
    eep    =  e     * ep
    eep3   =  eep   * ep * ep
    eep_2  =  eep   * eep
    eaep   =  e     + ep
    eaep2  =  eaep  * eaep
    eaep3  =  eaep2 * eaep

    IF ( l == 0 ) THEN

      a0   =  (r4_15*Ee5 - r4_3*Ee4*ep + r8_3*Ee3*ep2) / eep
      b0  = ( - a0 + r8_3*Ee2*(e3+ep3) &
                   - r4_3*Ee*eaep2*(ep2-2.0_dp*eep+3.0_dp*e2) &
                   + r4_15*eaep3*(ep2-3.0_dp*eep+6.0_dp*e2) ) &
            / eep 
      J_b = J_b + b0

    ELSE IF ( l == 1 ) THEN

      a1  = ( r16_35*Ee3*Ee4 - r4_5*Ee3*Ee3*(e+3.0_dp*ep) &
              + r4_15*Ee5*ep*(13.0_dp*e+18.0_dp*ep) &
              - r4_3*Ee4*ep2*(4.0_dp*e+3.0_dp*ep) &
              + r8_3*Ee3*eep3) / eep_2
      b1  = ( - a1 &
              - (r12_5*e5 + r4_3*e4*ep + r4_3*e*ep4 + r12_5*ep5) &
                 * Ee2 &
              + (r16_5*e3*e3 + r28_5*e5*ep + r8_3*e4*ep2 + &
                             + r28_15*e*ep5 + r8_5*ep3*ep3) &
                 * Ee &
              - (r8_7*e3*e4 + r16_5*e3*e3*ep + r16_5*e5*ep2 + r4_3*e4*ep3 &
                            + r8_15*e*ep3*ep3 + r12_35*ep3*ep4) &
             ) / eep_2 
      J_b = J_b + b1

      END IF


  END SUBROUTINE pr_w_gt_wp_w_lt_e


  SUBROUTINE pr_w_lt_e_e_lt_wp( e, ep, Ee, l, J_c )
!-----------------------------------------------
! To compute c0( e, ep, Ee )and c1( e, ep, Ee )
!-----------------------------------------------

  IMPLICIT NONE

    INTEGER,  INTENT(in)    :: l
    REAL(dp), INTENT(in)    :: e, ep, Ee
    REAL(dp), INTENT(inout) :: J_c

    REAL(dp)   :: c0, c1
    REAL(dp)   :: r4_3, r8_3, r8_5, r8_7, r12_5, r16_3, &
                  r16_5, r28_5
    REAL(dp)   :: e2, e3, ep2, ep3, Ee2, eep
   
    r4_3  = 4.0_dp  / 3.0_dp
    r8_3  = 8.0_dp  / 3.0_dp
    r8_5  = 8.0_dp  / 5.0_dp
    r8_7  = 8.0_dp  / 7.0_dp
    r12_5 = 12.0_dp / 5.0_dp
    r16_3 = 16.0_dp / 3.0_dp
    r16_5 = 16.0_dp / 5.0_dp
    r28_5 = 28.0_dp / 5.0_dp

    e2    = e   * e
    e3    = e2  * e
    ep2   = ep  * ep
    ep3   = ep2 * ep
    Ee2   = Ee  * Ee
    eep   = ep  * e


    IF ( l == 0 ) THEN

      c0  = ( r8_3 * ep2 + 4.0_dp * eep + r8_5 * e2 ) * e2 / ep &
            - ( r16_3 * e2 + 4 * e3 / ep ) * Ee &
            + r8_3 * Ee * Ee * e2 / ep
      J_c = J_c + c0

    ELSE IF ( l == 1 ) THEN

      c1  = - (r8_7*e3+r16_5*eep*e+r16_5*eep*ep+r4_3*ep3) * e2 /ep2 &
            + (r16_5*e2 + r28_5*eep + r8_3*ep2) * Ee * e2 / ep2 &
            - (r12_5*e + r4_3*ep ) * Ee2 * e2 / ep2  
      J_c = J_c + c1

      END IF

  END SUBROUTINE pr_w_lt_e_e_lt_wp


  SUBROUTINE pr_w_gt_e_e_gt_wp( e, ep, Ee, l, J_d )
!-----------------------------------------------
! To compute d0( e, ep, Ee )and d1( e, ep, Ee )
!-----------------------------------------------

  IMPLICIT NONE

    INTEGER,  INTENT(in)    :: l
    REAL(dp), INTENT(in)    :: e, ep, Ee
    REAL(dp), INTENT(inout) :: J_d

    REAL(dp)  :: d0, d1
    REAL(dp)  :: r4_3, r4_15, r8_3, r8_5, r8_15, r12_5, &
                 r12_35, r28_15
    REAL(dp)  :: ep2, ep3, ep4, e2, Ee2

    r4_3   = 4.0_dp  / 3.0_dp
    r4_15  = 4.0_dp  / 15.0_dp
    r8_3   = 8.0_dp  / 3.0_dp
    r8_5   = 8.0_dp  / 5.0_dp
    r8_15  = 8.0_dp  / 15.0_dp
    r12_5  = 12.0_dp / 5.0_dp
    r12_35 = 12.0_dp / 35.0_dp
    r28_15 = 28.0_dp / 15.0_dp

    ep2 = ep  * ep
    ep3 = ep2 * ep
    ep4 = ep3 * ep
    e2  = e   * e
    Ee2 = Ee  * Ee

    IF ( l == 0 ) THEN

      d0  = r4_15 * ep4 / e + r4_3 * ep3 * Ee / e &
           + r8_3 * ep2 * Ee * Ee / e
      J_d = J_d + d0

    ELSE IF ( l == 1 ) THEN

      d1  = - (r8_15*e + r12_35*ep ) * ep4 / e2 &
            + (r28_15*e + r8_5*ep ) * Ee * ep3 / e2 &
            - (r4_3*e + r12_5*ep ) * Ee2 * ep2 / e2
      J_d = J_d + d1

      END IF
    
  END SUBROUTINE pr_w_gt_e_e_gt_wp


  SUBROUTINE pr_w_lt_wp_e_lt_w( e, ep, Ee, l, J_a )
!-----------------------------------------------
! To compute a0( e, ep, Ee )and a1( e, ep, Ee )
!-----------------------------------------------

  IMPLICIT NONE

    INTEGER,  INTENT(in)    :: l
    REAL(dp), INTENT(in)    :: e, ep, Ee
    REAL(dp), INTENT(inout) :: J_a
   
    CALL pr_w_gt_wp_wp_gt_e( e, ep, Ee, l, J_a )

  END SUBROUTINE pr_w_lt_wp_e_lt_w


  SUBROUTINE pr_w_lt_wp_wp_lt_e( e, ep, Ee, l, J_b )
!-----------------------------------------------
! To compute b0( e, ep, Ee )and b1( e, ep, Ee )
!-----------------------------------------------

  IMPLICIT NONE
    INTEGER,  INTENT(in)    :: l
    REAL(dp), INTENT(in)    :: e, ep, Ee
    REAL(dp), INTENT(inout) :: J_b

    CALL pr_w_gt_wp_w_lt_e( e, ep, Ee, l, J_b )
    
  END SUBROUTINE pr_w_lt_wp_wp_lt_e
  
END MODULE B85_pair
