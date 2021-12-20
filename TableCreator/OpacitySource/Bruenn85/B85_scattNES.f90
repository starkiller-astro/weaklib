MODULE B85_scattNES
!---------------------------------------------------------------------
!
!    File:         B85_scattNES.f90
!    Module:       B85_scattNES
!    Type:         Module w/ Routines
!    Author:       R. Chu, Dept. Phys. & Astronomy
!                  U. Tennesee, Knoxville
!
!    Created:      8/16/18
!
!    WeakLib ver:  weaklib/External/Utilities/Opacities/Bruenn85
!
!    Purpose:      The zero and first legendre coefs for the n-type 
!                  neutrino electron scattering kernel considering the 
!                  physics in Bruenn 85.
!
!    CONTAINS:     Routine  TotalNESKernel
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
!                  NESK              
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

  PUBLIC TotalNESKernel

CONTAINS

  SUBROUTINE TotalNESKernel &
               ( energygrid, TMeV, chem_e, nquad, l, NESK, species )
!---------------------------------------------------------------------
! Purpose:
!   To compute the zero and first legendre coefs for the neutrino-
!   electron scattering kernel for electron-type neutrino and 
!   antineutrino 
! Ref: 
!   Mezzacappa & Bruenn(1993) AJ 410
!   Bruenn(1985) AJSS 58
!---------------------------------------------------------------------
  IMPLICIT NONE

  REAL(dp), DIMENSION(:), INTENT(in)      :: energygrid
  REAL(dp), INTENT(in)                    :: TMeV, chem_e
  INTEGER , INTENT(in)                    :: nquad, l, species
  REAL(dp), DIMENSION(:,:), INTENT(out)   :: NESK

  REAL(dp), DIMENSION(nquad)              :: roots, weights
  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: NESK_ome
  INTEGER                                 :: ii, jj, kk, nPointsE
  REAL(dp)                                :: outcome
  REAL(dp)                                :: UnitConvertConstant

  UnitConvertConstant  =  cvel_inv**4.0 / h**3.0
  nPointsE             =  SIZE(energygrid)

  ALLOCATE( NESK_ome( nquad, nPointsE, nPointsE ) )

  CALL gaquad( nquad, roots, weights, -1.0_dp , 1.0_dp )

  CALL NESKernelWithOmega &
         ( energygrid, roots, TMeV, chem_e, NESK_ome, species )

 
  IF ( l == 0 ) THEN

    DO jj = 1, nPointsE
      DO kk = 1, nPointsE

         outcome = zero
         DO ii = 1, nquad
           outcome = outcome + NESK_ome(ii,jj,kk) * weights(ii)     
         END DO
         NESK( jj, kk ) = UnitConvertConstant * outcome * 0.5_dp

      END DO
    END DO
  
  ELSE IF ( l == 1 ) THEN

    DO jj = 1, nPointsE
      DO kk = 1, nPointsE
         outcome = 0.0_dp
         DO ii = 1, nquad
           outcome = outcome &
                     + NESK_ome(ii,jj,kk) * weights(ii) * roots(ii)
         END DO
         NESK( jj, kk ) = UnitConvertConstant * outcome * 1.5_dp

      END DO
    END DO

  END IF

  END SUBROUTINE TotalNESKernel


  SUBROUTINE NESKernelWithOmega( energygrid, omega, TMeV, chem_e, nesktab, species )
!----------------------------------------------------------------------
! Purpose:
!    To compute the neutrino-electron scattering (OUT) kernel 
!    (1) e /= ep
!    R_out = cons * ( 1 / e / ep) * 
!                            ( beta1 * I1 + beta2 * I2 + beta3 * I3 )
!    (2) e == ep
!    new form of I1,2,3
!
!    species = 1 : electron-type neutrino
!    species = 2 : electron-type antineutrino
! Ref: 
!   Mezzacappa & Bruenn(1993) AJ 410
!   Bruenn(1985) AJSS 58
!----------------------------------------------------------------------
  IMPLICIT NONE

    REAL(dp), DIMENSION(:), INTENT(in) :: energygrid, omega
    INTEGER , INTENT(in)               :: species
    REAL(dp), INTENT(in) :: TMeV, chem_e       ! both unit = MeV
    REAL(dp), DIMENSION(:,:,:), INTENT(out) :: nesktab
    REAL(dp)             :: tsq, tinv, eta, FEXP, &
                            x1, x2, x2inv, &
                            FA_, FA0, FA1, &
                            y, tpiet, &
                            COMBO1, COMBO2, &
                            G0, G1, G2, buffer1, buffer2
    REAL(dp)             :: beta1, beta2, beta3, log10_sig0, log10_cons
    INTEGER              :: D_e, D_ome, i_e1, i_e2, i_ome
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: ediff, etap, fgamm, esq
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: delta, y0, I1, I2, I3,&
                                               A, B, C

       tsq = TMeV*TMeV
      tinv = 1.0/TMeV
       eta = chem_e * tinv
     tpiet = twpi * TMeV                      ! 2*pi*TMeV
     IF ( species == 1 ) THEN
       beta1 = ( cv + ca )**2 
       beta2 = ( cv - ca )**2
       beta3 = ca**2 - cv**2
     ELSE IF ( species == 2 ) THEN
       beta1 = ( cv - ca )**2
       beta2 = ( cv + ca )**2
       beta3 = ca**2 - cv**2
     END IF
     D_e = SIZE( energygrid )
     D_ome = SIZE( omega )
     log10_sig0  = LOG10(1.764) - 44  ! log10(sig0)            
     log10_cons  = LOG10(half * pi * cvel / ( twpi**3 * me**2 )) + log10_sig0

    ALLOCATE(  ediff( D_e, D_e )        )  ! e1 - e2 (ein -eout)
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
        
      CALL ComputeFValue_sameE( eta- y0(i_ome,i_e1,i_e1), FA_, FA0, FA1 )

      I1(i_ome,i_e1,i_e1)= &
         ( tpiet*esq(i_e1,i_e1)*esq(i_e1,i_e1)                      &
              *(1.0-omega(i_ome))*(1.0-omega(i_ome))                &
              /delta(i_ome,i_e1,i_e1)**5  )                         &
              *(A(i_ome,i_e1,i_e1)*tsq                              &
                            *(2.0*FA1                               &
                             +2.0*y0(i_ome,i_e1,i_e1)*FA0           &
                             +y0(i_ome,i_e1,i_e1)**2 *FA_           &
                             )                                      &
               +B(i_ome,i_e1,i_e1)*TMeV                                &
                            *(FA0                                   &
                             +y0(i_ome,i_e1,i_e1)*FA_               &
                             )                                      &
               +C(i_ome,i_e1,i_e1)*FA_                              &
               )      

      I2(i_ome,i_e1,i_e1) = I1(i_ome,i_e1,i_e1)                    

      I3(i_ome,i_e1,i_e1) =                                         &
         tpiet*esq(i_e1,i_e1)                                       &
              *(1.0-omega(i_ome))                                   &
              *me*me*FA_/delta(i_ome,i_e1,i_e1)            
    END DO ! i_e1
   END DO ! i_ome

!-------------------------------
!   energy_in \= energy_out
!-------------------------------
    DO i_ome = 1, D_ome
     DO  i_e1 = 1, D_e
      DO  i_e2 = 1, D_e

      IF( i_e1.ne.i_e2 ) THEN      

      CALL ComputeGValue_NESK( etap(i_e1,i_e2), &
                          eta, y0(i_ome,i_e1,i_e2), G0, G1, G2 )

      COMBO1=G2+2.0*y0(i_ome,i_e1,i_e2)*G1&
             +y0(i_ome,i_e1,i_e2)*y0(i_ome,i_e1,i_e2)*G0                                              
      COMBO2=G1+y0(i_ome,i_e1,i_e2)*G0                                                           

      I1(i_ome,i_e1,i_e2) =                                              &
         tpiet* esq(i_e1,i_e1)*esq(i_e2,i_e2)                            &
              *(1.0-omega(i_ome))*(1.0-omega(i_ome))                     &
              /delta(i_ome,i_e1,i_e2)**5                                 &
              *fgamm(i_e1,i_e2)                                          &
              *(A(i_ome,i_e1,i_e2)*tsq*COMBO1                            &
               +B(i_ome,i_e1,i_e2)*TMeV*COMBO2                              &
               +C(i_ome,i_e1,i_e2)*G0                                    &
               )                                                         
                                  
      I2(i_ome,i_e1,i_e2) =                                              &
         tpiet* esq(i_e1,i_e1)*esq(i_e2,i_e2)                            &
              *(1.0-omega(i_ome))*(1.0-omega(i_ome))                     &
              /delta(i_ome,i_e1,i_e2)**5                                 &
              *fgamm(i_e1,i_e2)                                          &
              *(A(i_ome,i_e1,i_e2)*tsq*COMBO1                            &
               -B(i_ome,i_e2,i_e1)*TMeV*COMBO2                              &
               +C(i_ome,i_e2,i_e1)*G0                                    &
               )                                                         
                                                                                
      I3(i_ome,i_e1,i_e2) =                                              &
         tpiet* esq(i_e1,i_e2)*(1.0-omega(i_ome))                        &
             *me*me*fgamm(i_e1,i_e2)*G0/delta(i_ome,i_e1,i_e2)             

      ELSE
      END IF

      buffer1                   = ( beta1 * I1(i_ome,i_e1,i_e2)    &
                                        + beta2 * I2(i_ome,i_e1,i_e2)    &
                                        + beta3 * I3(i_ome,i_e1,i_e2) )  &
                                      /  esq(i_e1,i_e2)
      buffer2 = LOG10(buffer1) + log10_cons 

      nesktab(i_ome,i_e1,i_e2) = 10**buffer2

      END DO ! i_e2
     END DO ! i_e1
    END DO ! i_ome

  END SUBROUTINE NESKernelWithOmega


  SUBROUTINE ComputeFValue_sameE( eta_y0, FA_, FA0, FA1) 

   IMPLICIT NONE

     REAL(dp), INTENT(in)    :: eta_y0
     REAL(dp), INTENT(out)   :: FA_, FA0, FA1
     REAL(dp)                 :: FEXP
     REAL(dp)                 :: x1, x2, x2inv
     REAL(dp)                 :: y, sumFA
     REAL(dp), DIMENSION(27)  :: Tarr, CH1arr
     INTEGER                  :: i_t

     FA_ = 1.0/( FEXP( -(eta_y0) ) + 1.0) 
  
     CH1arr(1)  = 1.935064300869969
     CH1arr(2)  = 0.166073032927855
     CH1arr(3)  = 0.024879329924228
     CH1arr(4)  = 0.004686361959447
     CH1arr(5)  = 0.001001627496164
     CH1arr(6)  = 0.000232002196094
     CH1arr(7)  = 0.000056817822718
     CH1arr(8)  = 0.000014496300557
     CH1arr(9)  = 0.000003816329463
     CH1arr(10) = 0.000001029904264
     CH1arr(11) = 0.000000283575385
     CH1arr(12) = 0.000000079387055
     CH1arr(13) = 0.000000022536705
     CH1arr(14) = 0.000000006474338
     CH1arr(15) = 0.000000001879117
     CH1arr(16) = 0.000000000550291
     CH1arr(17) = 0.000000000162421
     CH1arr(18) = 0.000000000048274
     CH1arr(19) = 0.000000000014437
     CH1arr(20) = 0.000000000004342
     CH1arr(21) = 0.000000000001312
     CH1arr(22) = 0.000000000000398
     CH1arr(23) = 0.000000000000121
     CH1arr(24) = 0.000000000000037
     CH1arr(25) = 0.000000000000011
     CH1arr(26) = 0.000000000000004
     CH1arr(27) = 0.000000000000001

     x1 = eta_y0
     x2 = -FEXP(x1)
     x2inv = 1.0/x2
      
     FA0 = LOG(1.0 - x2)

     IF(x2 .lt. -1.0) x2 = x2inv ! x1 .lt. 0
     y = (4.0*x2+1.0)/3.0

     Tarr(1) = 1.0  ! T0
     Tarr(2) = y    ! T1
     DO i_t = 3, SIZE(Tarr)
       Tarr(i_t) = 2.0*y*Tarr(i_t-1)-Tarr(i_t-2)
     END DO
    
     sumFA = 0.0
     DO i_t = 2, SIZE(Tarr)
       sumFA = CH1arr(i_t) * Tarr(i_t) + sumFA
     END DO
     FA1 = &
           -x2*( sumFA + 0.5*CH1arr(1)*Tarr(1) )

     IF( (1.0/x2inv) .lt. -1.0) FA1 = &  ! x2 .lt. -1
                               -FA1 + 0.5*x1**2 + pi*pi/6.0

   END SUBROUTINE ComputeFValue_sameE


   SUBROUTINE ComputeFValue_diffE&
              ( etap_y0, eta_y0, FA00, FA01, FA10, FA11, FA20, FA21)

   IMPLICIT NONE

     REAL(dp), INTENT(in)     :: etap_y0, eta_y0
     REAL(dp), INTENT(out)    :: FA00, FA01, FA10, FA11, FA20, FA21

     REAL(dp)                 :: x1, x2, x2inv, FEXP
     REAL(dp)                 :: y1, sumFA
     REAL(dp), DIMENSION(27)  :: Tarr, CH1arr
     REAL(dp), DIMENSION(24)  :: CH2arr
     INTEGER                  :: i_t

     CH1arr(1)  = 1.935064300869969
     CH1arr(2)  = 0.166073032927855
     CH1arr(3)  = 0.024879329924228
     CH1arr(4)  = 0.004686361959447
     CH1arr(5)  = 0.001001627496164
     CH1arr(6)  = 0.000232002196094
     CH1arr(7)  = 0.000056817822718
     CH1arr(8)  = 0.000014496300557
     CH1arr(9)  = 0.000003816329463
     CH1arr(10) = 0.000001029904264
     CH1arr(11) = 0.000000283575385
     CH1arr(12) = 0.000000079387055
     CH1arr(13) = 0.000000022536705
     CH1arr(14) = 0.000000006474338
     CH1arr(15) = 0.000000001879117
     CH1arr(16) = 0.000000000550291
     CH1arr(17) = 0.000000000162421
     CH1arr(18) = 0.000000000048274
     CH1arr(19) = 0.000000000014437
     CH1arr(20) = 0.000000000004342
     CH1arr(21) = 0.000000000001312
     CH1arr(22) = 0.000000000000398
     CH1arr(23) = 0.000000000000121
     CH1arr(24) = 0.000000000000037
     CH1arr(25) = 0.000000000000011
     CH1arr(26) = 0.000000000000004
     CH1arr(27) = 0.000000000000001

     CH2arr(1)  = 1.958417213383495
     CH2arr(2)  = 0.085188131486831
     CH2arr(3)  = 0.008559852220133
     CH2arr(4)  = 0.001211772144129
     CH2arr(5)  = 0.000207227685308
     CH2arr(6)  = 0.000039969586914 
     CH2arr(7)  = 0.000008380640658
     CH2arr(8)  = 0.000001868489447
     CH2arr(9)  = 0.000000436660867
     CH2arr(10) = 0.000000105917334
     CH2arr(11) = 0.000000026478920 
     CH2arr(12) = 0.000000006787000
     CH2arr(13) = 0.000000001776536 
     CH2arr(14) = 0.000000000473417
     CH2arr(15) = 0.000000000128121
     CH2arr(16) = 0.000000000035143
     CH2arr(17) = 0.000000000009754   
     CH2arr(18) = 0.000000000002736   
     CH2arr(19) = 0.000000000000775 
     CH2arr(20) = 0.000000000000221 
     CH2arr(21) = 0.000000000000064  
     CH2arr(22) = 0.000000000000018 
     CH2arr(23) = 0.000000000000005
     CH2arr(24) = 0.000000000000002

     x1 = eta_y0
     x2 = - FEXP(x1)
     x2inv = 1.0/x2

     FA00 = LOG(1.0 - x2)
     IF( x2 .lt. -1.0) x2 = x2inv  ! x1 .lt. 0
     y1 = (4.0*x2+1.0)/3.0

     Tarr(1) = 1.0
     Tarr(2) = y1
     DO i_t = 3, SIZE(Tarr)
       Tarr(i_t) = 2.0*y1*Tarr(i_t-1)-Tarr(i_t-2)
     END DO

     sumFA = 0.0
     DO i_t = 1, SIZE(CH2arr)
       sumFA = sumFA + CH2arr(i_t)*Tarr(i_t) 
     END DO

     FA20 = -2.0*x2*( sumFA - 0.5*CH2arr(1)*Tarr(1) )
    
     sumFA = 0.0 
     DO i_t = 1, SIZE(CH1arr)
       sumFA = sumFA + CH1Arr(i_t)*Tarr(i_t)
     END DO 

     FA10 = -x2*( sumFA - 0.5*CH1arr(1)*Tarr(1) )

     IF( (1.0/x2inv).lt.-1.0) THEN  ! x2 .lt. -1
       FA10 = -FA10 + 0.5*x1*x1+pi*pi/6.0
       FA20 =  FA20 + pi*pi*x1/3.0 + x1*x1*x1/3.0
     END IF

     x1 = etap_y0 
     x2 = - FEXP(x1)
     x2inv = 1.0/x2

     FA01 = LOG(1.0-x2)
     IF( x2 .lt. -1.0) x2 = x2inv
     y1 = (4.0*x2+1.0)/3.0

     Tarr(1) = 1.0
     Tarr(2) = y1
     DO i_t = 3, SIZE(Tarr)
       Tarr(i_t) = 2.0*y1*Tarr(i_t-1)-Tarr(i_t-2)
     END DO

     sumFA = 0.0
     DO i_t = 1, SIZE(CH2arr)
       sumFA = sumFA + CH2arr(i_t)*Tarr(i_t)
     END DO

     FA21 = -2.0*x2*( sumFA - 0.5*CH2arr(1)*Tarr(1) )

     sumFA = 0.0
     DO i_t = 1, SIZE(CH1arr)
       sumFA = sumFA + CH1Arr(i_t)*Tarr(i_t)
     END DO

     FA11 = -x2*( sumFA - 0.5*CH1arr(1)*Tarr(1) )

     IF( (1.0/x2inv).lt.-1.0) THEN
       FA11 = -FA11 + 0.5*x1*x1+pi*pi/6.0
       FA21 =  FA21 + pi*pi*x1/3.0 + x1*x1*x1/3.0
     END IF

   END SUBROUTINE ComputeFValue_diffE


  SUBROUTINE ComputeGValue_NESK( etap, eta, y0, G0, G1, G2 )

  IMPLICIT NONE

    REAL(dp), INTENT(in)      :: etap, eta, y0
    REAL(dp), INTENT(out)     :: G0, G1, G2
    REAL(dp)                 :: FA00, FA01, FA10, FA11, FA20, FA21

    CALL ComputeFValue_diffE( etap-y0, eta-y0, FA00, FA01, FA10, FA11, FA20, FA21)

    G0 = FA01 - FA00
    G1 = FA11 - FA10
    G2 = FA21 - FA20

  END SUBROUTINE ComputeGValue_NESK


END MODULE B85_scattNES
