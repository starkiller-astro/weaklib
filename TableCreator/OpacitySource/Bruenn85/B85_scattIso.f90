MODULE B85_scattIso
!---------------------------------------------------------------------
!
!    File:         B85_scattIso.f90
!    Module:       B85_scattIso
!    Type:         Module w/ Funcstions and Routines
!    Author:       R. Chu, Dept. Phys. & Astronomy
!                  U. Tennesee, Knoxville
!
!    Created:      8/16/18
!
!    WeakLib ver:  weaklib/External/Utilities/Opacities/Bruenn85
!
!    Purpose:      The zero and first legendre coefs for the n-type 
!                  neutrino isoenergetic scattering on nucleon and 
!                  nuclei kernel considering the physics in Bruenn 85.
!
!    CONTAINS:     Function  TotalIsoScatteringKernel
!
!    Input arguments:
!                  energy  : neutrino energy [MeV]
!                  rho     : matter density [g cm^{-3}]
!                  T       : matter temperature [K]
!                  xh      : heavy nucleus mass fraction
!                  A       : heavy nucleus mass number
!                  Z       : heavy nucleus charge number
!                  xn      : free neutron mass fraction
!                  xp      : free proton mass fraction
!                  l       : order of the legendre coefs/moment
!
!    Output arguments:
!                  TotalIsoScatteringKernel
!
!    Modules used:
!                  wlKindModule
!                  wlExtPhysicalConstantsModule
!                  wlExtNumericalModule
!                  ( function fexp is called )
!
!---------------------------------------------------------------------

  USE wlKindModule, ONLY: dp
  USE wlExtPhysicalConstantsModule, ONLY: &
      h, kMeV, therm1, therm2, dmnp, me, mbG, mp, mn, cvel_inv, &
      cvel, ergmev, cv_p, cv_n, ca_p, ca_n, gf, hbarc, cv, ca
  USE wlExtNumericalModule, ONLY: &
      pi, half, twpi, zero, one 

  IMPLICIT NONE
   
  PUBLIC TotalIsoScatteringKernel

CONTAINS

  REAL(dp) FUNCTION TotalIsoScatteringKernel &
                      ( energy, rho, T, xh, A, Z, xn, xp, l )
!---------------------------------------------------------------------
! Purpose:
!   To compute the zero and first legendre coefs for the any-type 
!   neutrino isoenergetic scattering kernel. On both nucleon and 
!   nuclei.
! Ref: 
!   Mezzacappa & Bruenn(1993) AJ 405
!   Bruenn(1985) AJSS 58
!---------------------------------------------------------------------
!--------------------------------------------------------------
!   Input Variables
!--------------------------------------------------------------
  USE, INTRINSIC :: ieee_arithmetic, ONLY: IEEE_IS_NAN

  IMPLICIT NONE
    REAL(dp), INTENT(in) :: energy, rho, T, xh, A, Z, xn, xp
    INTEGER, INTENT(in)  :: l
!--------------------------------------------------------------
!   Physical Constants 
!--------------------------------------------------------------
    REAL(dp)             :: N, nucleiTP, & !'TP' for thermal parameter
                            nucleonTP, nucleiExp, Cv0, Cv1,&
                            etann, etapp, Npara
!---------------------------------------------------------------------
!   Local Variables
!---------------------------------------------------------------------   
    REAL(dp) :: ISNucleiKernel_0, ISNucleonKernel_0, &
                ISNucleiKernel_1, ISNucleonKernel_1, &
                tempC0, TempC1
    INTEGER  :: nquad = 20

        N       =   A - Z
       Cv0      =   half * ( cv_p + cv_n) 
       Cv1      =   cv_p - cv_n
     Npara      =   twpi * cvel_inv**4.0 * energy**2.0 / h**3.0

    nucleiExp   =   4.0_dp * 4.8_dp * 10**(-6.0_dp) * &
                    A**(2.0_dp/3.0_dp) * energy**2.0     
    nucleiExp   =   MAX( nucleiExp, SQRT( TINY( 1.0_dp ) ) )

    IF ( xh == zero ) THEN

      nucleiTP  =   zero

    ELSE
      
      nucleiTP  =   ( (twpi*gf)**2 / h ) * ( rho*xh/mbG ) * &
                    A * ( Cv0 - ( (N-Z)*Cv1 )/(2.0_dp*A) )**2
    END IF

    nucleonTP   =   ( twpi * gf )**2 / h

!------------------------------
!  scattering on nuclei
!------------------------------

    tempC0 = IntegralESNuclei( nucleiExp, 0, nquad )
    tempC1 = IntegralESNuclei( nucleiExp, 1, nquad )

    ISNucleiKernel_0 = nucleiTP * tempC0 / 2.0_dp

    ISNucleiKernel_1 = nucleiTP * tempC1 * 3.0_dp / 2.0_dp

    IF ( IEEE_IS_NAN(ISNucleiKernel_0) .or. IEEE_IS_NAN(ISNucleiKernel_1)  ) THEN
     WRITE(*,*) "ERROR AT B85.f90 MARK 1003 !"
     WRITE(*,*) "nucleiExp is ", nucleiExp
     WRITE(*,*) "ISNucleiKernel_0 ", ISNucleiKernel_0
     WRITE(*,*) "ISNucleiKernel_1 ", ISNucleiKernel_1
     STOP
    END IF

!--------------------------------------
!   Scattering on Nucleons
!-------------------------------------

    CALL etaxx( rho, T, xn, xp, etann, etapp )

    ISNucleonKernel_0 = (2.0_dp) * nucleonTP * &
                        ( etann * ( cv_n**2 + 3.0_dp * ca_n**2) + &
                          etapp * ( cv_p**2 + 3.0_dp * ca_p**2) ) 
 
    ISNucleonKernel_1 = (2.0_dp / 3.0_dp) * nucleonTP * &
                        ( etann * ( cv_n**2 - ca_n**2) + &
                          etapp * ( cv_p**2 - ca_p**2) )

    IF ( IEEE_IS_NAN(ISNucleonKernel_0) .or. IEEE_IS_NAN(ISNucleonKernel_1)  ) &
    THEN
      WRITE(*,*) "ERROR AT B85.f90 MARK 1004 !"
      WRITE(*,*) "etann is ", etann
      WRITE(*,*) "etapp is ", etapp
      WRITE(*,*) "ISNucleonKernel_0 ", ISNucleonKernel_0
      WRITE(*,*) "ISNucleonKernel_1 ", ISNucleonKernel_1
      STOP
    END IF 

    IF ( l == 0 ) THEN
    
     TotalIsoScatteringKernel = Npara * ( ISNucleiKernel_0 &
                                        + ISNucleonKernel_0 )

    ELSE IF ( l == 1) THEN

     TotalIsoScatteringKernel = Npara * ( ISNucleiKernel_1 &
                                        + ISNucleonKernel_1 )

    ELSE

      WRITE(*,*) "ERROR: Unable to provide Legendre Moment with &
                         l other than 0 and 1 "
    END IF
    
    IF ( IEEE_IS_NAN(TotalIsoScatteringKernel) ) THEN
      WRITE(*,*) "TotalIsoScatteringKernel is NAN! "
      WRITE(*,*) "l is", l
      WRITE(*,*) "ISNucleiKernel_0 + ISNucleonKernel_0 ", &
                  ISNucleiKernel_0 + ISNucleonKernel_0 
      WRITE(*,*) "ISNucleiKernel_1 + ISNucleonKernel_1 ", &
                  ISNucleiKernel_1 + ISNucleonKernel_1 
      STOP
    END IF

    RETURN

  END FUNCTION TotalIsoScatteringKernel


  FUNCTION IntegralESNuclei( a, l, nquad )

  IMPLICIT NONE
    REAL(dp), INTENT(in)       :: a
    INTEGER,  INTENT(in)       :: l, nquad
  
    REAL(dp), DIMENSION(nquad) :: roots, weights
    REAL(dp)                   :: IntegralESNuclei, buffer, func
    INTEGER                    :: ii
  
    CALL gaquad( nquad, roots, weights, -1.0_dp , 1.0_dp )
  
    buffer = 0.0
  
    IF ( l == 0 ) THEN
   
      DO ii = 1, nquad
        func   =  ( 1.0 + roots(ii) ) * EXP( a * ( roots(ii) - 1.0 ) )
        buffer =  buffer + weights(ii) * func 
      END DO
  
    ELSE IF ( l == 1 ) THEN
  
      DO ii = 1, nquad
        func   =  roots(ii) * ( roots(ii) + 1.0 ) &
                  * EXP( a * ( roots(ii) - 1.0 ) )
        buffer =  buffer + weights(ii) * func
      END DO
  
    ELSE 

       WRITE(*,*) "ERROR when calling IntegralESNuclei function."

    END IF
  
    IntegralESNuclei = buffer
  
  END FUNCTION IntegralESNuclei


  SUBROUTINE etaxx( rho, T, xn, xp, etann, etapp )
  
  IMPLICIT NONE
!---------------------------------------------------------------------
!        Input variables.
!---------------------------------------------------------------------
  REAL(dp), INTENT(in)    :: rho      ! density (g/cm3)
  REAL(dp), INTENT(in)    :: T        ! temperature [K]
  REAL(dp), INTENT(in)    :: xn       ! neutron mass fraction
  REAL(dp), INTENT(in)    :: xp       ! proton mass fraction

!---------------------------------------------------------------------
!        Output variables.
!---------------------------------------------------------------------
  REAL(dp), INTENT(out)   :: etann    
                              ! neutron number corrected for blocking
  REAL(dp), INTENT(out)   :: etapp 
                              ! proton number corrected for blocking

!---------------------------------------------------------------------
!        Local variables
!--------------------------------------------------------------------
  REAL(dp)                :: nn            
                   ! neutron number uncorrected for blocking (cm^{-3})
  REAL(dp)                :: np  
                   ! proton number uncorrected for blocking (cm^{-3})
  REAL(dp)                :: d_n   
                   ! neutron number uncorrected for blocking (fm^{-3})
  REAL(dp)                :: d_p     
                   ! proton number uncorrected for blocking (fm^{-3})
  REAL(dp)                :: efn            ! degenerate expression
  REAL(dp)                :: efp            ! degenerate expression
  REAL(dp)                :: etanndgnt      ! nondegenerate expression
  REAL(dp)                :: etappdgnt      ! nondegenerate expression
  REAL(dp)                :: mpG, mnG
  REAL(dp), PARAMETER     :: tthird = 2.d0/3.d0
  
  mpG  =  mp * ergmev * cvel_inv * cvel_inv ! proton mass [g]
  mnG  =  mn * ergmev * cvel_inv * cvel_inv ! neutron mass [g]
  nn   =  xn * rho/mpG
  np   =  xp * rho/mnG

  IF ( nn < zero  .or.  np < zero ) THEN
    WRITE(*,*) "ERROR! nn or np less than zero." 
  END IF

!---------------------------------------------------------------------
!  etann, etanp (analytic approximation)
!---------------------------------------------------------------------
  IF ( nn == zero ) THEN

     etann      =  0.0_dp

  ELSE

     d_n        =  nn * 1.d-39
     efn        =  ( hbarc**2/( 2.d+00 * mn ) ) &
                   * ( 3.d+00 * pi**2 * d_n )**tthird
     etanndgnt  =  1.5d+00 * ( kMev * T/efn )
     etann      =  nn * etanndgnt/DSQRT( 1.d+00 + etanndgnt**2 )

  END IF

  IF ( np == zero ) THEN

     etapp      =  0.0_dp

  ELSE

     d_p        =  np * 1.d-39
     efp        =  ( hbarc**2/( 2.d+00 * mp ) ) &
                   * ( 3.d+00 * pi**2 * d_p )**tthird
     etappdgnt  =  1.5d+00 * ( kMev * T/efp )
     etapp      =  np * etappdgnt/DSQRT( 1.d+00 + etappdgnt**2 )

  END IF

  END SUBROUTINE etaxx

END MODULE B85_scattIso
