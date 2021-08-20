SUBROUTINE scatical_weaklib &
           ( n, egrid, nez, rho, t, xn, xp, xhe, xh, ah, zh, cok )
!--------------------------------------------------------------------
!
!    Author:  R. Chu, Dept. Phys. & Astronomy
!             U. Tennesee, Knoxville
!
!    Created: 10/23/18
!
!    Purpose:
!      To calculate the zero and first legendre coefs for the 
!      n-type isoenergetic scattering functions.
!
!    Subprograms called:
!      etaxx, scatiicr
!
!    Input arguments:
!
!   n         : neutrino type (1, e-neutrino; 2, e-antineutrino)
!   egrid     : energy grid
!   nez       : size of energy grid
!   rho       : matter density [g cm^{-3}]
!   t         : matter temperature [K]
!   xn        : mass fraction of free neutrons
!   xp        : mass fraction of free protons
!   xhe       : mass fraction of helium nuclei
!   xh        : mass fraction of heavy nuclei
!   ah        : mean heavy nucleus mass number
!   zh        : mean heavy nucleus charge number
!
!    Output arguments:
!
!   cok       : neutrino-nucleus scattering kerneal
!  
!-------------------------------------------------------------------

USE wlKindModule, ONLY: dp
USE numerical_module, ONLY: zero, half, twothd, one, epsilon, pi
USE physcnst_module, ONLY: Gw, mp, hbar, cvel, cv, ga, rmu

IMPLICIT NONE

!--------------------------------------------------------------------
!        Input variables.
!--------------------------------------------------------------------

INTEGER,      INTENT(in)     :: n       ! neutrino flavor index
INTEGER,      INTENT(in)     :: nez     ! number of energy groups
REAL(dp), INTENT(in)     :: rho     ! density [g cm^{-3}]
REAL(dp), INTENT(in)     :: t       ! temperature [K]
REAL(dp), INTENT(in)     :: xn      ! neutron mass fraction
REAL(dp), INTENT(in)     :: xp      ! proton mass fraction
REAL(dp), INTENT(in)     :: xhe     ! helium mass fraction
REAL(dp), INTENT(in)     :: xh      ! heavy nucleus mass fraction
REAL(dp), INTENT(in)     :: ah      ! heavy nucleus mass number
REAL(dp), INTENT(in)     :: zh      ! heavy nucleus charge number
REAL(dp), DIMENSION(nez), INTENT(in) :: egrid 
                                        ! neutrino energy grid[MeV] 

!--------------------------------------------------------------------
!        Output variables.
!--------------------------------------------------------------------

REAL(dp), INTENT(out), DIMENSION(nez,2) :: cok 
                                  ! legendre coefs scattering kernel
                                  !  cok(:,1) - zeroth
                                  !  cok(:,2) - first
!--------------------------------------------------------------------
!        Local variables
!--------------------------------------------------------------------
REAL(dp), PARAMETER      :: g2 = ( Gw/mp**2 )**2 * hbar**5 * cvel**6 ! (C40)
REAL(dp), PARAMETER      :: cc = ( 2.d0 * pi )/( cvel * ( 2.d0 * pi * hbar * cvel )**3 ) * ( 4.d0 * pi ) * g2
                             ! 0.5 * 4*pi/(c*c(hc)**3) * 4*pi*g2

REAL(dp), PARAMETER      :: hvp =  one - cv     ! (C31)
                             ! proton vector coupling constant
REAL(dp), PARAMETER      :: hap =  0.5d0 * ga   ! (C32)
                             ! proton axial vector coupling constant
REAL(dp), PARAMETER      :: hvn = -0.5d0        ! (C33)
                             ! neutron vector coupling constant
REAL(dp), PARAMETER      :: han = -0.5d0 * ga   ! (C34)
                             ! neutron axial vector coupling constant
REAL(dp), PARAMETER      :: cv0 =  half * ( hvp + hvn ) 
                             ! coupling constant
REAL(dp), PARAMETER      :: cv1 =  hvp - hvn      
                             ! coupling constant

REAL(dp), PARAMETER      :: ap0 =  hvp**2 + 3.d0 * hap**2 
                             ! proton zero moment coupling constant
REAL(dp), PARAMETER      :: ap1 =  hvp**2 - hap**2   
                             ! proton first moment coupling constant
REAL(dp), PARAMETER      :: an0 =  hvn**2 + 3.d0 * han**2
                             ! neutron zero moment coupling constant
REAL(dp), PARAMETER      :: an1 =  hvn**2 - han**2    
                             ! neutron first moment coupling constant

REAL(dp)                 :: eeche        
                             ! helium opacity calculation
REAL(dp)                 :: eche          
                             ! helium opacity calculation
REAL(dp)                 :: eche2         
                             ! helium opacity calculation
REAL(dp)                 :: eche3         
                             ! helium opacity calculation
REAL(dp)                 :: b0he          
                             ! helium opacity calculation
REAL(dp)                 :: saghe         
                             ! helium opacity calculation
REAL(dp)                 :: sbghe         
                             ! helium opacity calculation

REAL(dp)                 :: eec           
                             ! heavy nucleus opacity calculation
REAL(dp)                 :: ec            
                             ! heavy nucleus opacity calculation
REAL(dp)                 :: ec2           
                             ! heavy nucleus opacity calculation
REAL(dp)                 :: ec3          
                             ! heavy nucleus opacity calculation
REAL(dp)                 :: b0h           
                             ! heavy nucleus opacity calculation
REAL(dp)                 :: sag           
                             ! heavy nucleus opacity calculation
REAL(dp)                 :: sbg           
                             ! heavy nucleus opacity calculation
REAL(dp)                 :: xnh          
                             ! heavy nucleus opacity calculation

REAL(dp)                 :: aisov         
                             ! vector coupling constant
REAL(dp), PARAMETER      :: aisosc = cv0  
                             ! vector coupling constant
REAL(dp), PARAMETER      :: a01 = cv0 * cv0 
                             ! vector coupling constant
REAL(dp)                 :: a02           
                             ! coupling constant
REAL(dp)                 :: e2            
                             ! neutrino energy squared
REAL(dp)                 :: etann         
                             ! neutron number corrected for blocking
REAL(dp)                 :: xnn           
                             ! neutron number corrected for blocking
REAL(dp)                 :: etapp         
                             ! proton number corrected for blocking
REAL(dp)                 :: xnp           ! proton number corrected for blocking
REAL(dp)                 :: xnhe          ! helium number
REAL(dp)                 :: xheaa         ! heavy nucleus number * a_he^2
REAL(dp)                 :: xhaa          ! heavy nucleus number * ah^2

REAL(dp)                 :: ciicr         ! ion-ion correlation correction

REAL(dp), DIMENSION(nez) :: xi_p_wm       ! weak magnetism correction for neutrino-proton scattering
REAL(dp), DIMENSION(nez) :: xi_n_wm       ! weak magnetism correction for neutrino-neutron scattering
REAL(dp), DIMENSION(nez) :: xib_p_wm      ! weak magnetism correction for antineutrino-proton scattering
REAL(dp), DIMENSION(nez) :: xib_n_wm      ! weak magnetism correction for antineutrino-neutron scattering

REAL(dp), DIMENSION(nez) :: rmdnns   ! mfp^-1 for neutrino-neutron scattering
REAL(dp), DIMENSION(nez) :: rmdnns0  ! mfp^-1 for neutrino-neutron scattering
REAL(dp), DIMENSION(nez) :: rmdnns1  ! mfp^-1 for neutrino-neutron scattering
REAL(dp), DIMENSION(nez) :: rmdnps   ! mfp^-1 for neutrino-proton scattering
REAL(dp), DIMENSION(nez) :: rmdnps0  ! mfp^-1 for neutrino-proton scattering
REAL(dp), DIMENSION(nez) :: rmdnps1  ! mfp^-1 for neutrino-proton scattering
REAL(dp), DIMENSION(nez) :: rmdnbns  ! mfp^-1 for antineutrino-neutron scattering
REAL(dp), DIMENSION(nez) :: rmdnbns0 ! mfp^-1 for antineutrino-neutron scattering
REAL(dp), DIMENSION(nez) :: rmdnbns1 ! mfp^-1 for antineutrino-neutron scattering
REAL(dp), DIMENSION(nez) :: rmdnbps  ! mfp^-1 for antineutrino-proton scattering
REAL(dp), DIMENSION(nez) :: rmdnbps0 ! mfp^-1 for antineutrino-proton scattering
REAL(dp), DIMENSION(nez) :: rmdnbps1 ! mfp^-1 for antineutrino-proton scattering
REAL(dp), DIMENSION(nez) :: rmdnhes  ! mfp^-1 for neutrino-helium scattering
REAL(dp), DIMENSION(nez) :: rmdnhes0 ! mfp^-1 for neutrino-helium scattering
REAL(dp), DIMENSION(nez) :: rmdnhes1 ! mfp^-1 for neutrino-helium scattering
REAL(dp), DIMENSION(nez) :: rmdnhs   ! mfp^-1 for neutrino-heavy nucleus scattering
REAL(dp), DIMENSION(nez) :: rmdnhs0  ! mfp^-1 for neutrino-heavy nucleus scattering
REAL(dp), DIMENSION(nez) :: rmdnhs1  ! mfp^-1 for neutrino-heavy nucleus scattering

INTEGER                           :: k             ! energy loop counter
INTEGER                           :: nezl          ! local energy loop extent (nez or nez+nezext)

!--------------------------------------------------------------------
!  Initialize  mfp^-1
!--------------------------------------------------------------------

rmdnps             = zero
rmdnns             = zero
rmdnbps            = zero
rmdnbns            = zero
rmdnhes            = zero
rmdnhs             = zero
rmdnps0            = zero
rmdnns0            = zero
rmdnbps0           = zero
rmdnbns0           = zero
rmdnhes0           = zero
rmdnhs0            = zero
rmdnps1            = zero
rmdnns1            = zero
rmdnbps1           = zero
rmdnbns1           = zero
rmdnhes1           = zero
rmdnhs1            = zero

!!-------------------------------------------------------------------
!!  Calculate the weak magnetism corrections
!!-------------------------------------------------------------------
!
!CALL nc_weak_mag( unuloc, xi_p_wm, xi_n_wm, xib_p_wm, xib_n_wm, nez )
!
!--------------------------------------------------------------------
!  Quatities needed to compute the scattering functions.
!--------------------------------------------------------------------

CALL etaxx( rho, t, xn, xp, etann, etapp )
xnn              = etann
xnp              = etapp
xnhe             = ( xhe/4.d0             ) * rho/rmu
xnh              = ( xh /( ah + epsilon ) ) * rho/rmu
xheaa            = xnhe * 16.d0
xhaa             = xnh * ah * ah
b0he             = 1.21d-5
b0h              = 4.80d-6 * ( ah )**twothd
 
!--------------------------------------------------------------------
!  Zero neutrino-nucleus scattering rate if ah and zh below 1.d-10
!--------------------------------------------------------------------

IF ( ah < 1.d-10  .or.  zh < 1.d-10 ) THEN
  a02            = zero
ELSE
  aisov          = half * cv1 * ( 2.d0 * zh - ah )/ah
  a02            = ( aisosc + aisov ) * ( aisosc + aisov )
END IF ! ah < 1.d-10  .or.  zh < 1.d-10

DO k = 1, nez

  e2             = egrid(k) * egrid(k) 

!--------------------------------------------------------------------
!  Quantities needed to compute the coherent scattering on helium.
!--------------------------------------------------------------------
 
  eche           = 4.d0 * b0he * e2
  eche2          = eche * eche
  eche3          = eche * eche2

  IF ( eche < 0.1d0 ) THEN
    saghe        = 1.d0 - twothd * eche
    sbghe        = ( 5.d0 - eche2 )/15.d0
  ELSE
    eeche        = zero
  IF ( eche < 15.d0 ) eeche = DEXP( -2.d0 * eche )
    saghe        = ( eche - 0.5d0 + 0.5d0 * eeche )/eche2
    sbghe        = ( eche2 - 1.5d0 * eche + 1.d0 - ( 0.5d0 * eche + 1.d0 ) * eeche )/ eche3
  END IF ! eche < 0.1

!--------------------------------------------------------------------
!  Quantities needed to compute the coherent scattering on heavy 
!  nuclei.
!--------------------------------------------------------------------

  ec             = 4.d0 * b0h * e2
  ec2            = ec * ec
  ec3            = ec * ec2

  IF ( ec < 0.1d0 ) THEN
    sag          = 1.d0 - twothd * ec
    sbg          = ( 5.d0 - ec2 )/15.d0
  ELSE
    eec          = zero
    IF ( ec < 15.d0 ) eec = DEXP( -2.d0 * ec )
    sag          = ( ec - 0.5d0 + 0.5d0 * eec )/ec2
    sbg          = ( ec2 - 1.5d0 * ec + 1.d0 - ( 0.5d0 * ec + 1.d0 ) * eec )/ec3
  END IF ! ec < 0.1d0

!--------------------------------------------------------------------
!  Inverse mean free paths for coherent scattering.
!--------------------------------------------------------------------

  rmdnps0(k)     = cc * e2 * xnp * ap0 !4pi/(cvel*(hc)**3)*e2*(C38p)/2
  rmdnns0(k)     = cc * e2 * xnn * an0 !4pi/(cvel*(hc)**3)*e2*(C38n)/2
  rmdnbps0(k)    = rmdnps0(k)
  rmdnbns0(k)    = rmdnns0(k)
  rmdnhes0(k)    = cc * e2 * xheaa * a01 * saghe
                                            ! 4pi/(cvel*(hc)**3)*e2*(C44he)/2
  rmdnhs0(k)     = cc * e2 * xhaa * a02 * sag
                                            ! 4pi/(cvel*(hc)**3)*e2*(C44a)/2

  rmdnps1(k)     = cc * e2 * xnp * ap1
                                            ! 4pi/(cvel*(hc)**3)*e2*(C39p)*3/2
  rmdnns1(k)     = cc * e2 * xnn * an1
                                            ! 4pi/(cvel*(hc)**3)*e2*(C39n)*3/2
  rmdnbps1(k)    = rmdnps1(k)
  rmdnbns1(k)    = rmdnns1(k)
  rmdnhes1(k)    = cc * e2 * xheaa * a01 * sbghe * 3.0d0
                                            ! 4pi/(cvel*(hc)**3)*e2*(C45he)/2
  rmdnhs1(k)     = cc * e2 * xhaa * a02 * sbg * 3.0d0
                                            ! 4pi/(cvel*(hc)**3)*e2*(C45a)/2

!!-------------------------------------------------------------------
!!  Ion-ion correlation correction for coherent scattering.
!!-------------------------------------------------------------------
!
!  CALL scatiicr( rho, t, egrid(k), xh, ah, zh, ciicr )
!
!  rmdnhs(k)      = rmdnhs(k)  * ciicr
!  rmdnhs0(k)     = rmdnhs0(k) * ciicr
!  rmdnhs1(k)     = rmdnhs1(k) * ciicr
!
!!-------------------------------------------------------------------
!!  Weak magnetism corrections for neutrino and antineutrino neutron and
!!   proton scattering.
!!-------------------------------------------------------------------
!
!  rmdnps0 (k)    = rmdnps0 (k) * xi_p_wm (k)
!  rmdnns0 (k)    = rmdnns0 (k) * xi_n_wm (k)
!  rmdnbps0(k)    = rmdnbps0(k) * xib_p_wm(k)
!  rmdnbns0(k)    = rmdnbns0(k) * xib_n_wm(k)
!  rmdnps1 (k)    = rmdnps1 (k) * xi_p_wm (k)
!  rmdnns1 (k)    = rmdnns1 (k) * xi_n_wm (k)
!  rmdnbps1(k)    = rmdnbps1(k) * xib_p_wm(k)
!  rmdnbns1(k)    = rmdnbns1(k) * xib_n_wm(k)
!  rmdnps  (k)    = rmdnps  (k) * xi_p_wm (k)
!  rmdnns  (k)    = rmdnns  (k) * xi_n_wm (k)
!  rmdnbps (k)    = rmdnbps (k) * xib_p_wm(k)
!  rmdnbns (k)    = rmdnbns (k) * xib_n_wm(k)

END DO ! k = 1, nezl

!--------------------------------------------------------------------
!  Coherent (net) scattering inverse mean free path
!--------------------------------------------------------------------

SELECT CASE (n)
  CASE (1) ! e-neutrino
     cok(:,1) = rmdnps0  + rmdnns0  + rmdnhes0 + rmdnhs0
     cok(:,2) = rmdnps1  + rmdnns1  + rmdnhes1 + rmdnhs1
  CASE (2) ! e-antineutrino
     cok(:,1) = rmdnbps0 + rmdnbns0 + rmdnhes0 + rmdnhs0
     cok(:,2) = rmdnbps1 + rmdnbns1 + rmdnhes1 + rmdnhs1
  CASE DEFAULT
     cok = zero
END SELECT

RETURN

END SUBROUTINE scatical_weaklib
