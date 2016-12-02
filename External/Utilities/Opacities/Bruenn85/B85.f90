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

  
  SUBROUTINE NESKern( energygrid, omega, T, chem_e, nesktab )
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

    REAL(dp), DIMENSION(:), INTENT(in) :: energygrid, omega
    REAL(dp), INTENT(in) :: T, chem_e      
    REAL(dp), DIMENSION(:,:,:), INTENT(out) :: nesktab
    REAL(dp)             :: tsq, tinv, eta, FEXP, &
                            x1, x2, x2inv, &
                            FA0, FA1, &
                            y, tpiet, &
                            COMBO1, COMBO2, &
                            G0, G1, G2
    REAL(dp)             :: beta1, beta2, beta3, sig0, cons
    INTEGER              :: D_e, D_ome, i_e1, i_e2, i_ome
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: ediff, etap, fgamm, esq
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: delta, y0, I1, I2, I3,&
                                               A, B, C

       tsq = T*T
      tinv = 1.0/T
       eta = chem_e * tinv
     tpiet = twpi * T                      ! 2*pi*T
     beta1 = ( cv + ca )**2 
     beta2 = ( cv - ca )**2
     beta3 = ca**2 - cv**2
       D_e = SIZE( energygrid )
     D_ome = SIZE( omega )
     sig0  = 1.764 * 10**(-44)             ! cm**(-2)
     cons  = half * pi * sig0 * cvel / ( twpi**3 * me**2 )

    ALLOCATE(  ediff( D_e, D_e )        )  ! e1 - e2
    ALLOCATE(   etap( D_e, D_e )        )  ! eta prim
    ALLOCATE(   esq ( D_e, D_e )        )  ! e*e
    ALLOCATE(  fgamm( D_e, D_e )        )  ! gamma function
    ALLOCATE( delta ( D_ome, D_e, D_e ) )
    ALLOCATE(    y0 ( D_ome, D_e, D_e ) )
    ALLOCATE(     A ( D_ome, D_e, D_e ) )
    ALLOCATE(     B ( D_ome, D_e, D_e ) )
    ALLOCATE(     C ( D_ome, D_e, D_e ) )
    ALLOCATE(     I1( D_ome, D_e, D_e ) )
    ALLOCATE(     I2( D_ome, D_e, D_e ) )
    ALLOCATE(     I3( D_ome, D_e, D_e ) )

    DO i_e1 = 1, D_e
     DO i_e2 = 1, D_e

       ediff(i_e1,i_e2) = energygrid( i_e1 ) - energygrid( i_e2 )
        etap(i_e1,i_e2) = eta - ediff(i_e1,i_e2)*tinv
         esq(i_e1,i_e2) = energygrid( i_e1 ) * energygrid( i_e2 )
         
     IF ( i_e1 .ne. i_e2 ) THEN
       fgamm(i_e1,i_e2) = 1.0/(FEXP(-ediff(i_e1,i_e2)*tinv)-1.0)
     ELSE
     END IF     

     END DO ! i_e2
    END DO ! i_e1

    DO i_ome = 1, D_ome
     DO  i_e1 = 1, D_e
      DO  i_e2 = 1, D_e

       delta(i_ome,i_e1,i_e2) =                                          &
                             SQRT( esq(i_e1,i_e1) + esq(i_e2,i_e2)       &
                             - 2.0*esq(i_e1,i_e2)*omega(i_ome)  )
       y0(i_ome,i_e1,i_e2) =                                             &
                             tinv*(-0.5*ediff(i_e1,i_e2)                 &
                              +0.5*delta( i_ome, i_e1, i_e2 )            &
                                *sqrt(1.0+2.0*me*me                      &
                                /(esq(i_e1,i_e2)*(1.0-omega(i_ome)))     & 
                                     )                                   &
                                  )
       A(i_ome,i_e1,i_e2) = esq(i_e1,i_e1) + esq(i_e2,i_e2)              &
                          + esq(i_e1,i_e2)*(3.0 + omega(i_ome))

       B(i_ome,i_e1,i_e2) = energygrid(i_e1) *                           &
                            ( 2.0*esq(i_e1,i_e1)                         &
                            + esq(i_e1,i_e2)*(3.0 - omega(i_ome))        &
                            - esq(i_e2,i_e2)*(1.0 + 3.0*omega(i_ome) ) ) 
       C(i_ome,i_e1,i_e2) = esq(i_e1,i_e1) *                             &
                            ( (energygrid(i_e1) -                        &
                                  energygrid(i_e2)*omega(i_ome) )**2     &
                             - half*esq(i_e2,i_e2)*(1.0-omega(i_ome)**2) &
                             - half*( 1.0+omega(i_ome) )*me*me           &
                                         *(delta(i_ome,i_e1,i_e2)**2)    &
                               /( (1.0-omega(i_ome)) * esq(i_e1,i_e1) ) )  

      END DO ! i_e2
     END DO ! i_e1
    END DO ! i_ome  

!-------------------------------
!   energy_in == energy_out
!-------------------------------
   DO i_ome = 1, D_ome
    DO  i_e1 = 1, D_e
        
!       x1 = eta-y0
!       x2 = -FEXP(x1)
!       x2inv = 1.0/x2
!
!       FA0   = LOG(1.0-x2)
!
!       IF(x2.lt.-1.0) x2 = x2inv
!       y = (4.0*x2+1.0)/3.0
!
!      T0=1.0                                                                    
!      T1=y                                                                      
!      T2=2.0*y*T1-T0                                                            
!      T3=2.0*y*T2-T1                                                            
!      T4=2.0*y*T3-T2                                                            
!      T5=2.0*y*T4-T3                                                            
!      T6=2.0*y*T5-T4                                                            
!      T7=2.0*y*T6-T5                                                            
!      T8=2.0*y*T7-T6                                                            
!      T9=2.0*y*T8-T7                                                            
!      T10=2.0*y*T9-T8                                                           
!      T11=2.0*y*T10-T9                                                          
!      T12=2.0*y*T11-T10                                                         
!      T13=2.0*y*T12-T11                                                         
!      T14=2.0*y*T13-T12                                                         I
!      T15=2.0*y*T14-T13                                                         
!      T16=2.0*y*T15-T14                                                         
!      T17=2.0*y*T16-T15                                                         
!      T18=2.0*y*T17-T16                                                         
!      T19=2.0*y*T18-T17                                                         
!      T20=2.0*y*T19-T18                                                         
!      T21=2.0*y*T20-T19                                                         
!      T22=2.0*y*T21-T20                                                         
!      T23=2.0*y*T22-T21                                                         
!      T24=2.0*y*T23-T22                                                         
!      T25=2.0*y*T24-T23                                                         
!      T26=2.0*y*T25-T24 
!                                                                                
!      FA1= &
!           -x2*(0.5*CH10*T0+CH1(1)*T1+CH1(2)*T2+CH1(3)*T3+CH1(4)*T4 &
!                           +CH1(5)*T5+CH1(6)*T6+CH1(7)*T7+CH1(8)*T8 &
!                           +CH1(9)*T9+CH1(10)*T10+CH1(11)*T11       &
!                           +CH1(12)*T12+CH1(13)*T13+CH1(14)*T14     &
!                           +CH1(15)*T15+CH1(16)*T16+CH1(17)*T17     &
!                           +CH1(18)*T18+CH1(19)*T19+CH1(20)*T20     &
!                           +CH1(21)*T21+CH1(22)*T22+CH1(23)*T23     &
!                           +CH1(24)*T24+CH1(25)*T25+CH1(26)*T26     &
!             )                                                                  
!
!      if((1.0/x2inv).lt.-1.0) FA1= &
!                             -FA1+0.5*x1**2+pie**2/6.0
!
      I1(i_ome,i_e1,i_e1)= &
         ( tpiet*esq(i_e1,i_e1)*esq(i_e1,i_e1)                      &
              *(1.0-omega(i_ome))*(1.0-omega(i_ome))                &
              /delta(i_ome,i_e1,i_e1)**5  )                         &
              *(A(i_ome,i_e1,i_e1)*tsq                              &
                            *(2.0*FA1                               &
                             +2.0*y0(i_ome,i_e1,i_e1)*FA0           &
                             +y0(i_ome,i_e1,i_e1)**2                &
                                /(FEXP(-(eta                        &
                                  -y0(i_ome,i_e1,i_e1)))+1.0)       &
                             )                                      &
               +B(i_ome,i_e1,i_e1)*T                                &
                            *(FA0                                   &
                             +y0(i_ome,i_e1,i_e1)                   &
                                /(FEXP(-(eta                        &
                                  -y0(i_ome,i_e1,i_e1)))+1.0)       &
                             )                                      &
               +C(i_ome,i_e1,i_e1)/(FEXP(-(eta                      &
                                  -y0(i_ome,i_e1,i_e1)))+1.0)       &
               )      

      I2(i_ome,i_e1,i_e1) =  I1(i_ome,i_e1,i_e1)                    

      I3(i_ome,i_e1,i_e1) =                                         &
         tpiet*esq(i_e1,i_e1)                                       &
              *(1.0-omega(i_ome))                                   &
              *me*me/delta(i_ome,i_e1,i_e1)                         &
              /(FEXP(-(eta-y0(i_ome,i_e1,i_e1)))+1.0)             
                 
    END DO ! i_e1
   END DO ! i_ome

!-------------------------------
!   energy_in \= energy_out
!-------------------------------
    DO i_ome = 1, D_ome
     DO  i_e1 = 1, D_e
      DO  i_e2 = 1, D_e

      IF( i_e1.ne.i_e2 ) THEN      
!                                                                                
!      x3=eta-yo                                                              
!      x4=etap(k,m)-yo 
!
!      x1=x3
!      x2=-fexp(x1)                                                               
!      x2inv=1.0/x2 
!
!      FA00=log(1.0-x2)                                                           
!
!      if(x2.lt.-1.0) x2=x2inv 
!      y1=(4.0*x2+1.0)/3.0       
!
!      T0=1.0                                                                    
!      T1=y1                                                                     
!      T2=2.0*y1*T1-T0                                                           
!      T3=2.0*y1*T2-T1                                                           
!      T4=2.0*y1*T3-T2                                                           
!      T5=2.0*y1*T4-T3                                                           
!      T6=2.0*y1*T5-T4                                                           
!      T7=2.0*y1*T6-T5                                                           
!      T8=2.0*y1*T7-T6                                                           
!      T9=2.0*y1*T8-T7                                                           
!      T10=2.0*y1*T9-T8                                                          
!      T11=2.0*y1*T10-T9                                                         
!      T12=2.0*y1*T11-T10                                                        
!      T13=2.0*y1*T12-T11                                                        
!      T14=2.0*y1*T13-T12                                                        
!      T15=2.0*y1*T14-T13                                                        
!      T16=2.0*y1*T15-T14                                                        
!      T17=2.0*y1*T16-T15                                                        
!      T18=2.0*y1*T17-T16                                                        
!      T19=2.0*y1*T18-T17                                                        
!      T20=2.0*y1*T19-T18                                                        
!      T21=2.0*y1*T20-T19                                                        
!      T22=2.0*y1*T21-T20                                                        
!      T23=2.0*y1*T22-T21  
!
!      F2=-2.0*x2*(0.5*CH20*T0+CH2(1)*T1+CH2(2)*T2+CH2(3)*T3+CH2(4)*T4    & 
!                               +CH2(5)*T5+CH2(6)*T6+CH2(7)*T7+CH2(8)*T8  &
!                               +CH2(9)*T9+CH2(10)*T10+CH2(11)*T11        &
!                               +CH2(12)*T12+CH2(13)*T13+CH2(14)*T14      &
!                               +CH2(15)*T15+CH2(16)*T16+CH2(17)*T17      &
!                               +CH2(18)*T18+CH2(19)*T19+CH2(20)*T20      &
!                               +CH2(21)*T21+CH2(22)*T22+CH2(23)*T23      &
!                 )      
!
!      T24=2.0*y1*T23-T22                                                        
!      T25=2.0*y1*T24-T23                                                        
!      T26=2.0*y1*T25-T24                                                        
!                                                                                
!      F1=-x2*(0.5*CH10*T0+CH1(1)*T1+CH1(2)*T2+CH1(3)*T3+CH1(4)*T4        &
!                           +CH1(5)*T5+CH1(6)*T6+CH1(7)*T7+CH1(8)*T8      &
!                           +CH1(9)*T9+CH1(10)*T10+CH1(11)*T11            &
!                           +CH1(12)*T12+CH1(13)*T13+CH1(14)*T14          &
!                           +CH1(15)*T15+CH1(16)*T16+CH1(17)*T17          &
!                           +CH1(18)*T18+CH1(19)*T19+CH1(20)*T20          &
!                           +CH1(21)*T21+CH1(22)*T22+CH1(23)*T23          &
!                           +CH1(24)*T24+CH1(25)*T25+CH1(26)*T26          &
!             )                     
!
!      if((1.0/x2inv).lt.-1.0) then
!        F1=-F1+0.5*x1**2+pie**2/6.0
!        F2=F2+pie**2*x1/3.0+x1**3/3.0
!      endif
!
!      FA10=F1                                                                    
!      FA20=F2
!
!      x1=x4
!      x2=-fexp(x1)                                                               
!      x2inv=1.0/x2 
!
!      FA01=log(1.0-x2)                                                           
!  
!      if(x2.lt.-1.0) x2=x2inv 
!      y1=(4.0*x2+1.0)/3.0                                                       
!
!      T0=1.0                                                                    
!      T1=y1                  
!      T2=2.0*y1*T1-T0                                                           
!      T3=2.0*y1*T2-T1                                                           
!      T4=2.0*y1*T3-T2                                                           
!      T5=2.0*y1*T4-T3                                                           
!      T6=2.0*y1*T5-T4                                                           
!      T7=2.0*y1*T6-T5                                                           
!      T8=2.0*y1*T7-T6                                                           
!      T9=2.0*y1*T8-T7                                                           
!      T10=2.0*y1*T9-T8                                                          
!      T11=2.0*y1*T10-T9                                                         
!      T12=2.0*y1*T11-T10                                                        
!      T13=2.0*y1*T12-T11                                                        
!      T14=2.0*y1*T13-T12                                                        
!      T15=2.0*y1*T14-T13                                                        
!      T16=2.0*y1*T15-T14                                                        
!      T17=2.0*y1*T16-T15                                                        
!      T18=2.0*y1*T17-T16                                                        
!      T19=2.0*y1*T18-T17                                                        
!      T20=2.0*y1*T19-T18                                                        
!      T21=2.0*y1*T20-T19                                                        
!      T22=2.0*y1*T21-T20                                                        
!      T23=2.0*y1*T22-T21    
!
!      F2=-2.0*x2*(0.5*CH20*T0+CH2(1)*T1+CH2(2)*T2+CH2(3)*T3+CH2(4)*T4    &        
!                               +CH2(5)*T5+CH2(6)*T6+CH2(7)*T7+CH2(8)*T8  &
!                               +CH2(9)*T9+CH2(10)*T10+CH2(11)*T11        &
!                               +CH2(12)*T12+CH2(13)*T13+CH2(14)*T14      &
!                               +CH2(15)*T15+CH2(16)*T16+CH2(17)*T17      &
!                               +CH2(18)*T18+CH2(19)*T19+CH2(20)*T20      &
!                               +CH2(21)*T21+CH2(22)*T22+CH2(23)*T23      &
!                 )      
!
!      T24=2.0*y1*T23-T22                                                        
!      T25=2.0*y1*T24-T23                                                        
!      T26=2.0*y1*T25-T24                                                        
!                                                                                
!      F1=-x2*(0.5*CH10*T0+CH1(1)*T1+CH1(2)*T2+CH1(3)*T3+CH1(4)*T4        &
!                           +CH1(5)*T5+CH1(6)*T6+CH1(7)*T7+CH1(8)*T8      &
!                           +CH1(9)*T9+CH1(10)*T10+CH1(11)*T11            &
!                           +CH1(12)*T12+CH1(13)*T13+CH1(14)*T14          &
!                           +CH1(15)*T15+CH1(16)*T16+CH1(17)*T17          &
!                           +CH1(18)*T18+CH1(19)*T19+CH1(20)*T20          &
!                           +CH1(21)*T21+CH1(22)*T22+CH1(23)*T23          &
!                           +CH1(24)*T24+CH1(25)*T25+CH1(26)*T26          &
!             )                                                                  
!
!      if((1.0/x2inv).lt.-1.0) then
!        F1=-F1+0.5*x1**2+pie**2/6.0
!        F2=F2+pie**2*x1/3.0+x1**3/3.0
!      endif
!
!      FA11=F1                                                                    
!      FA21=F2
!
      G0=1.0!FA01-FA00
      G1=2.0!FA11-FA10
      G2=3.0!FA21-FA20
!                                                                                
      COMBO1=1.0!G2+2.0*y0*G1+y0*y0*G0                                              
      COMBO2=2.0!G1+y0*G0                                                           
!
      I1(i_ome,i_e1,i_e2) =                                              &
         tpiet* esq(i_e1,i_e1)*esq(i_e2,i_e2)                            &
              *(1.0-omega(i_ome))*(1.0-omega(i_ome))                     &
              /delta(i_ome,i_e1,i_e2)**5                                 &
              *fgamm(i_e1,i_e2)                                          &
              *(A(i_ome,i_e1,i_e2)*tsq*COMBO1                            &
               +B(i_ome,i_e1,i_e2)*T*COMBO2                              &
               +C(i_ome,i_e1,i_e2)*G0                                    &
               )                                                         
                                  
      I2(i_ome,i_e1,i_e2) =                                              &
         tpiet* esq(i_e1,i_e1)*esq(i_e2,i_e2)                            &
              *(1.0-omega(i_ome))*(1.0-omega(i_ome))                     &
              /delta(i_ome,i_e1,i_e2)**5                                 &
              *fgamm(i_e1,i_e2)                                          &
              *(A(i_ome,i_e1,i_e2)*tsq*COMBO1                            &
               -B(i_ome,i_e2,i_e1)*T*COMBO2                              &
               +C(i_ome,i_e2,i_e1)*G0                                    &
               )                                                         
                                                                                
      I3(i_ome,i_e1,i_e2) =                                              &
         tpiet* esq(i_e1,i_e2)*(1.0-omega(i_ome))                        &
             *me*me*fgamm(i_e1,i_e2)*G0/delta(i_ome,i_e1,i_e2)             

      ELSE
      END IF

      nesktab(i_ome,i_e1,i_e2) = cons * ( beta1 * I1(i_ome,i_e1,i_e2)    &
                                        + beta2 * I2(i_ome,i_e1,i_e2)    &
                                        + beta3 * I3(i_ome,i_e1,i_e2) )  &
                                      /  esq(i_e1,i_e2)
      END DO ! i_e2
     END DO ! i_e1
    END DO ! i_ome

   END SUBROUTINE NESKern

!  FUNCTION Fgamm( x )
!    REAL(dp), INTENT(in) :: x
!    REAL(dp)             :: Fgamm, FEXP
!
!    Fgamm = 1.0/( FEXP(x) - 1.0   )
!    RETURN
!  END FUNCTION Fgamm

!  FUNCTION FerInt( n, eta )
!
!    INTEGER,  INTENT(in) :: n
!    REAL(dp), INTENT(in) :: eta
!    REAL(dp)             :: FerInt
!    FerInt = 0.0
!    RETURN
!
!  END FUNCTION FerInt
  

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
