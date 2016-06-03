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
!
!    Modules used:
!       
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
                            midE, chem_v, inversefeq, mpG, mnG, feq

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
    
    
!    etapn = rho * ( xn - xp ) / ( mbG * ( EXP( (chem_n-chem_p)/TMeV ) - 1.0_dp ) )
    etapn = rho * xp  / mpG                  ! Approxiation in the nondegenerate regime
    
    IF ( (etapn < 0.0_dp ) ) THEN
      WRITE(*,*)'etapn is negtive: ', etapn
      WRITE(*,*)'xn - xp is ', xn - xp
      WRITE(*,*),'exp term is ',  EXP( (chem_n-chem_p)/TMeV ) - 1.0_dp 
      WRITE(*,*)'chem_n - chem_p is ', chem_n - chem_p 
      STOP 
    END IF

    midFe = 1.0_dp / ( EXP( (energy+dmnp-chem_e) / TMeV ) + 1.0_dp )

    midE = (energy+dmnp)**2 &
                 * SQRT( 1.0_dp - ( me / (energy+dmnp) )**2 )

    jnucleon = therm1 * etapn * midE * midFe

    IF ( xheavy * npz * nhn == 0.0_dp ) THEN
      jnuclear = 0.0_dp
    ELSE
      jnuclear = therm2 * rho * xheavy * npz * nhn * (energy+qpri)**2 &
                 * SQRT( 1.0_dp - ( me / (energy+qpri) )**2 ) &
                 / ( mbG * a * ( EXP( (energy+qpri-chem_e) / TMeV ) + 1_dp ) )
    END IF

    totalECapEm = (jnucleon + jnuclear) / feq
 
    IF ( ISNAN(totalECapEm ) ) THEN    
      WRITE(*,*) 'totalECapEm is ', totalECapEm
      WRITE(*,*) 'LOG10(totalECapEm) is ', LOG10(totalECapEm)
      WRITE(*,*) 'jnucleon term is ', jnucleon
      WRITE(*,*) 'jnuclear term is ', jnuclear
      WRITE(*,*) 'therm2 * rho * xheavy * npz * nhn = ', &
                                          therm2 * rho * xheavy * npz * nhn
      WRITE(*,*) 'SQRT( 1.0_dp - ( me / (energy+qpri) )**2 ) = ', &
                                     SQRT( 1.0_dp - ( me / (energy+qpri) )**2 )
      WRITE(*,*) '( me / (energy+qpri) )**2 = ', ( me / (energy+qpri) )**2 
      WRITE(*,*) 'qpri = chem_n + dmnp +3.0_dp - chem_p is ', qpri
      WRITE(*,*) 'inversefeq is ', inversefeq
      WRITE(*,*) 'feq is ', feq
      WRITE(*,*) 'tiny (1.0) is', TINY( 1.0_dp )
      WRITE(*,*) ''
      STOP
    END IF
    RETURN 

  END FUNCTION totalECapEm

END MODULE B85
