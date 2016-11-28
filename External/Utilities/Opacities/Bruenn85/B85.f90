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
    h, kMeV, therm1, therm2, dmnp, me, mbG, mp, mn, cvel_inv, cvel, ergmev,&
      cv_p, cv_n, ca_p, ca_n, gf, hbarc, cv, ca
  USE wlExtNumericalModule, ONLY: pi, half, twpi, zero 

  implicit none
   
  PUBLIC totalECapEm, &
         totalElasticScatteringKernel

CONTAINS

!========================Function=============================

  REAL(dp) FUNCTION totalECapEm &
    ( energy, rho, T, Z, A, chem_e, chem_n, chem_p, xheavy, xn, xp )
!------------------------------------------------------------------------------
! Purpose:
!   To compute the neutrino absorptivity.
!   (1) Absorptivity = emissivity + inverse of mean path
!------------------------------------------------------------------------------
  IMPLICIT NONE

    REAL(dp), INTENT(in) :: energy, rho, T, Z, A, chem_e, chem_n, chem_p, &
                            xheavy, xn, xp

    REAL(dp) :: TMeV, n, qpri, nhn, npz, etapn, jnucleon, jnuclear, midFe, &
                            midE, chem_v, mpG, mnG, fexp, &
                            rop, ron, midFexpp, midFep, midEp, midCons
    REAL(dp) :: emitnp, absornp, emitni, absorni

    TMeV   = T * kMeV                         ! kmev = 8.61733d-11 [MeV K^{-1}]
      N    = A - Z
    qpri   = chem_n - chem_p + 3.0_dp + dmnp  ! [MeV] 3 = energy of the 1f5/2 level 
    chem_v = chem_e + chem_p - chem_n - dmnp  ! neutrino chemical potential
    
   
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
       midEp   = (energy+qpri)**2 * SQRT( & 
                    MAX( 1.0_dp - ( me / (energy+qpri) )**2, 0.0_dp ) )     
      midFexpp = FEXP( (energy+qpri-chem_e) / TMeV )
       midFep  = 1.0_dp / ( midFexpp + 1.0_dp )
      midCons  = therm2 * rho * xheavy * npz * nhn * midEp / (mbG * A)

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

    IF ( ISNAN(totalECapEm) ) THEN
      WRITE(*,*) "totalECapEm is NAN! "
      STOP
    END IF

    RETURN
  END FUNCTION totalECapEm


  REAL(dp) FUNCTION totalElasticScatteringKernel&
     ( energy, rho, T, xh, A, Z, xn, xp, l )

!------------------------------------------------------------------------------
! Purpose:
!   To compute the zero and first legendre coefs for the neutrino-electron 
!   elastic scattering kernel. 
!------------------------------------------------------------------------------
!-----------------------------------------------------------------------
!   Input Variables
!-----------------------------------------------------------------------
    REAL(dp), INTENT(in) :: energy, rho, T, xh, A, Z, xn, xp
    INTEGER, INTENT(in)  :: l
!-----------------------------------------------------------------------
!   Physical Constants 
!-----------------------------------------------------------------------
    REAL(dp)             :: N, nucleiTP, & ! 'TP' for thermal parameter
                            nucleonTP, nucleiExp, Cv0, Cv1,&
                            etann, etapp, Npara

!-----------------------------------------------------------------------
!   Local Variables
!-----------------------------------------------------------------------   

    REAL(dp) :: ESNucleiKernel_0, ESNucleonKernel_0, &
                ESNucleiKernel_1, ESNucleonKernel_1, &
                tempC0, TempC1
    INTEGER  :: nquad = 20

        N     = A - Z
       Cv0    = half * ( cv_p + cv_n) 
       Cv1    = cv_p - cv_n
     Npara    = cvel_inv**4.0 * energy**2.0 / h**3.0
    
    nucleiExp = 4.0_dp * 4.8_dp * 10**(-6.0_dp) * &
                A**(2.0_dp/3.0_dp) * energy**2.0     
    nucleiExp = MAX( nucleiExp, SQRT( TINY( 1.0_dp ) ) )

    nucleiTP  = ( (twpi*gf)**2 / h ) * ( rho*xh/mbG ) * &
                A * ( Cv0 - ( (N-Z)*Cv1 )/(2.0_dp*A) )**2

    nucleonTP = ( twpi * gf )**2 / h

!------------------------------
!  scattering on nuclei
!------------------------------

    tempC0 = IntegralESNuclei( nucleiExp, 0, nquad )
    tempC1 = IntegralESNuclei( nucleiExp, 1, nquad )

    ESNucleiKernel_0 = nucleiTP * tempC0 / 2.0_dp

    ESNucleiKernel_1 = nucleiTP * tempC1 * 3.0_dp / 2.0_dp

    IF ( ISNAN(ESNucleiKernel_0) .or. ISNAN(ESNucleiKernel_1)  ) THEN
     WRITE(*,*) "ERROR AT B85.f90 MARK 1003 !"
     WRITE(*,*) "nucleiExp is ", nucleiExp
     WRITE(*,*) "ESNucleiKernel_0 ", ESNucleiKernel_0
     WRITE(*,*) "ESNucleiKernel_1 ", ESNucleiKernel_1
     STOP
    END IF

!--------------------------------------
!   Scattering on Nucleons
!-------------------------------------

    CALL etaxx( rho, T, xn, xp, etann, etapp )

    ESNucleonKernel_0 = (0.5_dp) * nucleonTP * &
                        ( etann * ( cv_n**2 + 3.0_dp * ca_n**2) + &
                          etapp * ( cv_p**2 + 3.0_dp * ca_p**2) )
 
    ESNucleonKernel_1 = (1.5_dp) * nucleonTP * &
                        ( etann * ( cv_n**2 - ca_n**2) + &
                          etapp * ( cv_p**2 - ca_p**2) )

    IF ( ISNAN(ESNucleonKernel_0) .or. ISNAN(ESNucleonKernel_1)  ) THEN
     WRITE(*,*) "ERROR AT B85.f90 MARK 1004 !"
     STOP
    END IF 

    IF ( l == 0 ) THEN
    
     totalElasticScatteringKernel = Npara * ( ESNucleiKernel_0 &
                                            + ESNucleonKernel_0 )

    ELSE IF ( l == 1) THEN

     totalElasticScatteringKernel = Npara * ( ESNucleiKernel_1 &
                                            + ESNucleonKernel_1 )

    ELSE

     WRITE(*,*) "ERROR: Unable to provide Legendre Moment with &
                        l other than 0 and 1 "
    END IF
    
    IF ( ISNAN(totalElasticScatteringKernel) ) THEN
      WRITE(*,*) "totalElasticScatteringKernel is NAN! "
      WRITE(*,*) "l is", l
      WRITE(*,*) "ESNucleiKernel_0 + ESNucleonKernel_0 ", ESNucleiKernel_0+ESNucleonKernel_0 
      WRITE(*,*) "ESNucleiKernel_1 + ESNucleonKernel_1 ", ESNucleiKernel_1+ESNucleonKernel_1 
      STOP
    END IF

    RETURN

  END FUNCTION totalElasticScatteringKernel

  
  REAL(dp) FUNCTION NESKern( energy1, energy2, omega, T, chem_e )  ! energy value or energy index?
!----------------------------------------------------------------------
! Purpose:
!    To compute the neutrino-electron scattering (OUT) kernel 
!    (1) e /= ep
!    R_out = cons * ( 1 / e / ep) * 
!                            ( beta1 * I1 + beta2 * I2 + beta3 * I3 )
!    (2) e == ep
!    new form of I1,2,3
!----------------------------------------------------------------------
  IMPLICIT NONE

    REAL(dp), INTENT(in) :: energy1, energy2, omega, T, chem_e      ! value of omega
    REAL(dp)             :: beta1, beta2, beta3
    REAL(dp)             :: cons, I1, I2, I3, sig0, delta, y0, A, B, C
    REAL(dp)             :: e, ep, tinv, eta, etap
    REAL(dp)             :: Gfun2, Gfun1, Gfun0

!    INTEGER, PARAMETER   :: jmax=8, & !number of angle bins
!                            nesp=20   !order of Guass integral
!    REAL(dp), DIMENSION(jmax, jmax, nesp ) :: omega
    
    tinv  = 1.0/T
     eta  = chem_e / T
      e   = energy1
      ep  = energy2
     etap = eta - tinv * (e-ep)
    beta1 = ( cv + ca )**2 
    beta2 = ( cv - ca )**2
    beta3 = ca**2 - cv**2
    sig0  = 1.764 * 10**(-44)           ! cm**(-2)
    cons  = half * pi * sig0 * cvel / ( twpi**3 * me**2 )

    delta = ( e*e + ep*ep - 2*e*ep*omega ) **(1/2)
      y0  = tinv * (- half * ( e - ep )  &
                    + half * delta &
                       * (1.0 + 2.0*(me*me)/( e*ep*(1-omega) ))**(1/2) )
      A   = e*e + ep*ep + e*ep*( 3.0 + omega )
      B   = e * ( 2.0*e*e + e*ep*(3.0-omega) - ep*ep * (1+3.0*omega) )
      C   = e*e * ( (e - ep*omega)**2 - half*ep**2 * (1-omega**2) &
                     - half*(1+omega)*me*me*delta*delta/((1-omega)*e*e)  )

    Gfun2 = FerInt( 2, etap-y0) - FerInt( 2, eta-y0) 
    Gfun1 = FerInt( 1, etap-y0) - FerInt( 1, eta-y0)
    Gfun0 = FerInt( 0, etap-y0) - FerInt( 0, eta-y0)

    IF( e .ne. ep ) THEN
       I1 = twpi * T *e*e * ep*ep* (1.0-omega)*(1.0-omega) * Fgamm((ep - e)/T) &
               * ( A*T*T*( Gfun2 + 2.0*y0*Gfun1 + y0*y0*Gfun0 ) &
                 + B*T*( Gfun1 +   y0*Gfun0 ) &
                 + C*Gfun0 ) &
             / delta**5
       I3 = twpi*T*e*ep*(1.0-omega)*me*me*Fgamm((ep-e)/T)*Gfun0/delta
       e  = -energy2
       ep = -energy1
       B  = e * ( 2.0*e*e + e*ep*(3.0-omega) - ep*ep * (1.0+3.0*omega) )
       C  = e*e * ( (e - ep*omega)**2 - half*ep**2 * (1.0-omega**2) &
                   - half*(1.0+omega)*me*me*delta*delta/((1.0-omega)*e*e)  )
       I2 = twpi * T *e*e * ep*ep* (1-omega)*(1-omega) * Fgamm((ep - e)/T) &
             * ( A*T*T*( Gfun2 + 2*y0*Gfun1 + y0*y0*Gfun0 ) &
                 + B*T*( Gfun1 +   y0*Gfun0 ) &
                 + C*Gfun0 ) &

    ELSE IF ( e .eq. ep ) THEN
       I1 = 1.0
       I2 = I1
       I3 = twpi*T*e*e*(1.0 - omega)*me*me*FerInt(-1, eta-y0)/delta
    END IF
    
    NESKern = cons * ( beta1 * I1 + beta2 * I2 + beta3 * I3 ) / ( e * ep )

  END FUNCTION NESKern

  FUNCTION Fgamm( x )
    REAL(dp), INTENT(in) :: x
    REAL(dp)             :: Fgamm, FEXP

    Fgamm = 1.0/( FEXP(x) - 1.0   )
    RETURN
  END FUNCTION Fgamm

  FUNCTION FerInt( n, eta )

    INTEGER,  INTENT(in) :: n
    REAL(dp), INTENT(in) :: eta
    REAL(dp)             :: FerInt
    FerInt = 0.0
    RETURN

  END FUNCTION FerInt
  

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


  FUNCTION IntegralESNuclei( a, l, nquad )

  REAL(dp), INTENT(in)       :: a
  INTEGER,  INTENT(in)       :: l, nquad

  REAL(dp), DIMENSION(nquad) :: roots, weights
  REAL(dp)                   :: IntegralESNuclei, buffer, func
  INTEGER                    :: ii

  CALL gaquad( nquad, roots, weights, -1.0_dp , 1.0_dp )

  buffer = 0.0

  IF ( l == 0 ) THEN
 
   DO ii = 1, nquad
     func = ( 1.0 + roots(ii) )* EXP( a * ( roots(ii) - 1.0 ) )
     buffer = buffer + weights(ii) * func 
   END DO

  ELSE IF ( l == 1 ) THEN

   DO ii = 1, nquad
     func = roots(ii) * ( roots(ii) + 1.0 ) * EXP ( a * ( roots(ii) - 1.0 ) )
     buffer = buffer + weights(ii) * func
   END DO

  ELSE 
    WRITE(*,*) "ERROR when calling IntegralESNuclei function."
  END IF

  IntegralESNuclei = buffer

  END FUNCTION IntegralESNuclei

END MODULE B85
