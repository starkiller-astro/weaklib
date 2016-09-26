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
!
!      Function totalECapEm: 
!                   gives opacity with np as approximation 
!
!      Function totalElasticScatteringKernel: 
!                   gives isoenergetic scattering kernel with nMoment
!                   l = 0 or 1
!                            
!
!    Modules used:
!      wlKindModule
!      wlExtPhysicalConstantsModule
!      ( function fexp is called )
!-----------------------------------------------------------------------
  USE wlKindModule, ONLY: dp
  USE wlExtPhysicalConstantsModule, ONLY: &
    h, kMeV, therm1, therm2, dmnp, me, mbG, mp, mn, cvel_inv, ergmev,&
      cv_p, cv_n, ca_p, ca_n, gf, hbarc
  USE wlExtNumericalModule, ONLY: pi, half, twpi, zero 

  implicit none
   
  PUBLIC totalECapEm, &
         totalElasticScatteringKernel

CONTAINS

!========================Function=============================

  REAL(dp) FUNCTION totalECapEm &
    ( energy, rho, T, Z, A, chem_e, chem_n, chem_p, xheavy, xn, xp )

    REAL(dp), INTENT(in) :: energy, rho, T, Z, A, chem_e, chem_n, chem_p, &
                            xheavy, xn, xp

    REAL(dp) :: TMeV, n, qpri, nhn, npz, etapn, jnucleon, jnuclear, midFe, &
                            midE, chem_v, inversefeq, mpG, mnG, feq, fexp, &
                            rop, ron, midFexpp, midFep, midEp, midCons
    REAL(dp) :: emitnp, absornp, emitni, absorni

    TMeV   = T * kMeV                         ! kmev = 8.61733d-11 [MeV K^{-1}]
      N    = A - Z
    qpri   = chem_n - chem_p + 3.0_dp + dmnp  ! [MeV] 3 = energy of the 1f5/2 level 
    chem_v = chem_e + chem_p - chem_n - dmnp  ! neutrino chemical potential
    
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
     rop   = rho * xp / mpG                   ! Approxiation in the nondegenerate regime
     ron   = rho * xn / mnG 
!----------------------------------------------------------------------------
! Stop the function if any of etapn/rop/ron is negative
!--------------------------------------------------------------------------- 
    IF ( ( rop < 0.0_dp ) .or. ( ron < 0.0_dp ) ) THEN
!    IF ( (etapn < 0.0_dp ) .or. ( rop < 0.0_dp ) .or. ( ron < 0.0_dp ) ) THEN
      WRITE(*,*)'xn - xp is ', xn - xp
      WRITE(*,*)'fexp term is ', FEXP( (chem_n-chem_p)/TMeV ) - 1.0_dp 
      WRITE(*,*)'chem_n - chem_p is ', chem_n - chem_p 
      STOP 
    END IF

!-----------------------------------------------------------------------------
!   j_nuclear(emitni) and chi_nuclear(absorni) 
!-----------------------------------------------------------------------------
    IF ( xheavy * npz * nhn == 0.0_dp ) THEN
       emitni  = 0.0_dp
      absorni  = 0.0_dp
    ELSE
      midCons  = therm2 * rho * xheavy * npz * nhn * midEp / (mbG * a)
       midEp   = (energy+qpri)**2 * SQRT( 1.0_dp - ( me / (energy+qpri) )**2 )       
      midFexpp = FEXP( (energy+qpri-chem_e) / TMeV )
       midFep  = 1.0_dp / ( midFexpp + 1.0_dp )

       emitni  = midCons * midFep 
      absorni  = midCons * FEXP( (chem_n + dmnp - chem_p - qpri) /TMeV ) &
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

      midFe   = 1.0_dp / ( FEXP( (energy+dmnp-chem_e) / TMeV ) + 1.0_dp )
       midE   = (energy+dmnp)**2 &
                 * SQRT( 1.0_dp - ( me / (energy+dmnp) )**2 )

     jnucleon = therm1 * rop * midE * midFe
     absornp  = therm1 * ron * midE * ( 1.0_dp - midFe )
       emitnp = jnucleon

    totalECapEm = ( emitni + absorni ) + ( emitnp + absornp )

    RETURN
  END FUNCTION totalECapEm


  REAL(dp) FUNCTION totalElasticScatteringKernel&
     ( energy, rho, T, xh, A, Z, xn, xp, l )

!-----------------------------------------------------------------------
!   Input Variables
!-----------------------------------------------------------------------
    REAL(dp), INTENT(in) :: energy, rho, T, xh, A, Z, xn, xp
    INTEGER, INTENT(in)  :: l
!-----------------------------------------------------------------------
!   Physical Constants 
!-----------------------------------------------------------------------
    REAL(dp)             :: N, nucleiTP, & ! 'TP' for thermal parameter
                            nucleonTP, nucleiExp, Cv0, Cv1, etann, etapp

!-----------------------------------------------------------------------
!   Local Variables
!-----------------------------------------------------------------------   

    REAL(dp) :: ESNucleiKernel_0, ESNucleonKernel_0, &
                ESNucleiKernel_1, ESNucleonKernel_1

        N     = A - Z
       Cv0    = half * ( cv_p + cv_n) 
       Cv1    = cv_p - cv_n 
    
    nucleiExp = 4.0_dp * 4.8_dp * 10**(-6.0_dp) * &
                A**(2.0_dp/3.0_dp) * energy**2     
    nucleiExp = MAX( nucleiExp, SQRT( TINY( 1.0_dp ) ) )

    nucleiTP  = ( (twpi*gf)**2 / h ) * ( rho*xh/mbG ) * &
                A * ( Cv0 - ( (N-Z)*Cv1 )/(2.0_dp*A) )**2 * &
                EXP(- nucleiExp )

    CALL etaxx( rho, T, xn, xp, etann, etapp )

    nucleonTP = ( twpi * gf )**2 / h
  
    ESNucleiKernel_0  = (0.5_dp) * nucleiTP * &
                        ( EXP(nucleiExp) - EXP(-nucleiExp) ) &
                        / nucleiExp
 
    ESNucleiKernel_1 = (1.5_dp) * nucleiTP * &
                       ( ABS( ( nucleiExp - 1.0_dp) * EXP( nucleiExp ) &
                          + ( nucleiExp + 1.0_dp) * EXP( -nucleiExp ) ) ) &
                       / ( nucleiExp**2.0_dp )

    ESNucleonKernel_0 = (0.5_dp) * nucleonTP * &
                        ( etann * ( cv_n**2 + 3.0_dp * ca_n**2) + &
                          etapp * ( cv_p**2 + 3.0_dp * ca_p**2) )
 
    ESNucleonKernel_1 = (1.5_dp) * nucleonTP * &
                        ( etann * ( cv_n**2 - ca_n**2) + &
                          etapp * ( cv_p**2 - ca_p**2) )

    IF ( l == 0 ) THEN
    
     totalElasticScatteringKernel = ESNucleiKernel_0 + ESNucleonKernel_0

    ELSE IF ( l == 1) THEN

     totalElasticScatteringKernel = ESNucleiKernel_1 + ESNucleonKernel_1

    ELSE

     WRITE(*,*) "ERROR: Unable to provide Legendre Moment with &
                        l other than 0 and 1 "
    END IF

    RETURN

  END FUNCTION totalElasticScatteringKernel


  SUBROUTINE etaxx( rho, T, xn, xp, etann, etapp )
  

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

  REAL(dp), INTENT(in)    :: rho           ! density (g/cm3)
  REAL(dp), INTENT(in)    :: T             ! temperature [K]
  REAL(dp), INTENT(in)    :: xn            ! neutron mass fraction
  REAL(dp), INTENT(in)    :: xp            ! proton mass fraction

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

  REAL(dp), INTENT(out)   :: etann         ! neutron number corrected for blocking
  REAL(dp), INTENT(out)   :: etapp         ! proton number corrected for blocking

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

  REAL(dp)                :: nn            ! neutron number uncorrected for blocking (cm^{-3})
  REAL(dp)                :: np            ! proton number uncorrected for blocking (cm^{-3})
  REAL(dp)                :: d_n            ! neutron number uncorrected for blocking (fm^{-3})
  REAL(dp)                :: d_p            ! proton number uncorrected for blocking (fm^{-3})
  REAL(dp)                :: efn           ! degenerate expression
  REAL(dp)                :: efp           ! degenerate expression
  REAL(dp)                :: etanndgnt     ! nondegenerate expression
  REAL(dp)                :: etappdgnt     ! nondegenerate expression
  REAL(dp)                :: mpG, mnG
  REAL(dp), PARAMETER     :: tthird = 2.d0/3.d0
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  nn, np
!-----------------------------------------------------------------------

  
  mpG   = mp * ergmev * cvel_inv * cvel_inv ! proton mass [g]
  mnG   = mn * ergmev * cvel_inv * cvel_inv ! neutron mass [g]
  nn                 = xn * rho/mpG
  np                 = xp * rho/mnG

  IF ( nn <= zero  .or.  np <= zero ) THEN
    WRITE(*,*) "ERROR! nn or np less than zero."
  END IF

!-----------------------------------------------------------------------
!  etann, etanp (analytic approximation)
!-----------------------------------------------------------------------

  d_n          = nn * 1.d-39
  d_p          = np * 1.d-39
  efn          = ( hbarc**2/( 2.d+00 * mn ) ) &
                    * ( 3.d+00 * pi**2 * d_n )**tthird
  efp          = ( hbarc**2/( 2.d+00 * mp ) ) &
                    * ( 3.d+00 * pi**2 * d_p )**tthird
  etanndgnt    = 1.5d+00 * ( kMev * T/efn )
  etappdgnt    = 1.5d+00 * ( kMev * T/efp )
  etann        = nn * etanndgnt/DSQRT( 1.d+00 + etanndgnt**2 )
  etapp        = np * etappdgnt/DSQRT( 1.d+00 + etappdgnt**2 )

  END SUBROUTINE etaxx

END MODULE B85
