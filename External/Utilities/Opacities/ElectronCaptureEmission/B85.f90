MODULE B85
!-----------------------------------------------------------------------
!
!    File:         B85.90
!    Module:       B85
!    Type:         Module w/ Functions
!    Author:       R. Chu, Dept. Phys. & Astronomy
!                  U. Tennesee, Knoxville
!
!    Created:      3/31/16
!    WeakLib ver:  
!
!    Purpose:
!      Provides the function need considering the physics in Bruenn 85
!
!   
!
!    CONTAINS:
!      Function totalECapEm: gives opacity with np as approximation 
!
!    Modules used:
!      wlKindModule
!      wlExtPhysicalConstantsModule
!      ( function fexp is called )
!-----------------------------------------------------------------------
  USE wlKindModule, ONLY: dp
  USE wlExtPhysicalConstantsModule, ONLY: &
    kMeV, therm1, therm2, dmnp, me, mbG, mp, mn, cvel_inv, ergmev
 

  implicit none
   
  PUBLIC totalECapEm

CONTAINS
!========================Function=============================

  REAL(dp) FUNCTION &
    totalECapEm( energy, rho, T, Z, A, chem_e, chem_n, chem_p, xheavy, xn, xp )

    REAL(dp), INTENT(in) :: energy, rho, T, Z, A, chem_e, chem_n, chem_p, &
                            xheavy, xn, xp

    REAL(dp) :: TMeV, n, qpri, nhn, npz, etapn, jnucleon, jnuclear, midFe, &
                            midE, chem_v, inversefeq, mpG, mnG, feq, fexp, &
                            rop, ron, midFexpp, midFep, midEp, midCons
    REAL(dp) :: emitnp, absornp, emitni, absorni

    TMeV   = T * kmev                 ! kmev = 8.61733d-11 [MeV K^{-1}]
      N    = A - Z
    qpri   = chem_n - chem_p + 3.0_dp + dmnp! [MeV] 3MeV: energy of the 1f5/2 level
    chem_v = chem_e + chem_p - chem_n - dmnp! neutrino chemical potential
    
    inversefeq   = ( EXP( (energy - chem_v) / TMeV ) + 1.0_dp )   
     feq   = MAX( 1.0_dp / ( EXP( (energy - chem_v) / TMeV ) + 1.0_dp ),&
                  SQRT( TINY( 1.0_dp ) ) )   
   
     mpG   = mp * ergmev * cvel_inv * cvel_inv ! proton mass [g]
     mnG   = mn * ergmev * cvel_inv * cvel_inv ! neutron mass [g]

    if(n.le.34.0)               nhn = 6.0_dp
    if(n.gt.34.0.and.n.le.40.0) nhn = 40.0_dp - N
    if(n.gt.40.0)               nhn = 0.0_dp

    if(z.le.20.0)               npz = 0.0
    if(z.gt.20.0.and.z.le.28.0) npz = z - 20.0
    if(z.gt.28.0)               npz = 8.0
    
    
!    etapn = rho * ( xn - xp ) / ( mbG * ( FEXP( (chem_n-chem_p)/TMeV ) - 1.0_dp ) )
    rop = rho * xp / mpG                  ! Approxiation in the nondegenerate regime
    ron = rho * xn / mpG 
!----------------------------------------------------------------------------
! Stop the function if any of etapn/rop/ron is negative
!--------------------------------------------------------------------------- 
    IF ( ( rop < 0.0_dp ) .or. ( ron < 0.0_dp ) ) THEN
!    IF ( (etapn < 0.0_dp ) .or. ( rop < 0.0_dp ) .or. ( ron < 0.0_dp ) ) THEN
      WRITE(*,*)'xn - xp is ', xn - xp
      WRITE(*,*),'fexp term is ', FEXP( (chem_n-chem_p)/TMeV ) - 1.0_dp 
      WRITE(*,*)'chem_n - chem_p is ', chem_n - chem_p 
      STOP 
    END IF

!-----------------------------------------------------------------------------
!   j_nuclear(emitni) and chi_nuclear(absorni) 
!-----------------------------------------------------------------------------
    IF ( xheavy * npz * nhn == 0.0_dp ) THEN
      emitni = 0.0_dp
      absorni  = 0.0_dp
    ELSE
      midCons = therm2 * rho * xheavy * npz * nhn * midEp / (mbG * a)
      midEp = (energy+qpri)**2 * SQRT( 1.0_dp - ( me / (energy+qpri) )**2 )       
      midFexpp = FEXP( (energy+qpri-chem_e) / TMeV )
      midFep = 1.0_dp / ( midFexpp + 1.0_dp )

      emitni = midCons * midFep 
      absorni = midCons * FEXP( (chem_n + dmnp - chem_p - qpri) /TMeV ) &
               * ( 1.0_dp - midFep) 
    END IF
!-----------------------------------------------------------------------------
!   j_nucleon(emitnp) and chi_nucleon(absornp)
!-----------------------------------------------------------------------------
!------------------------------------------------------------------
!  Set emitnp + absornp = zero and return if both xn and xp are zero
!------------------------------------------------------------------
    IF ( xn == 0.0 .and. xp == 0.0 ) THEN
      totalECapEm = emitni + absorni
      WRITE(*,*) 'xn and xp = 0'
      RETURN
    END IF

    midFe = 1.0_dp / ( FEXP( (energy+dmnp-chem_e) / TMeV ) + 1.0_dp )
    midE = (energy+dmnp)**2 &
                 * SQRT( 1.0_dp - ( me / (energy+dmnp) )**2 )

    jnucleon = therm1 * rop * midE * midFe
    absornp  = therm1 * ron * midE * ( 1.0_dp - midFe )
    emitnp = jnucleon

    totalECapEm = ( emitni + absorni ) + ( emitnp + absornp )

    RETURN

  END FUNCTION totalECapEm

END MODULE B85
