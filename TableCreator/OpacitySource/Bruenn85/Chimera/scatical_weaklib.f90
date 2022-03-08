SUBROUTINE scatical_weaklib &
           ( n, egrid, nez, rho, t, xn, xp, xhe, xh, ah, zh, cok )
!--------------------------------------------------------------------
!
!    Author:  R. Chu, Dept. Phys. & Astronomy
!             U. Tennesee, Knoxville
!
!    Created: 10/23/18
!    Edited : 12/15/19   --- Dropped the 2*pi from coeff.
!             03/08/22   --- Move legendre coeffs to creator driver
!
!    Purpose:
!      To calculate the zero and first legendre coefs for the 
!      n-type isoenergetic scattering functions without 1/2 and 3/2.
!      (1/2 and 3/2 are moved to creator driver)
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

USE kind_module, ONLY: double
USE numerical_module, ONLY: zero, half, twothd, one, epsilon, pi
USE physcnst_module, ONLY: Gw, mp, hbar, cvel, cv, ga, rmu

IMPLICIT NONE

!--------------------------------------------------------------------
!        Input variables.
!--------------------------------------------------------------------

INTEGER,      INTENT(in)     :: n       ! neutrino flavor index
INTEGER,      INTENT(in)     :: nez     ! number of energy groups
REAL(double), INTENT(in)     :: rho     ! density [g cm^{-3}]
REAL(double), INTENT(in)     :: t       ! temperature [K]
REAL(double), INTENT(in)     :: xn      ! neutron mass fraction
REAL(double), INTENT(in)     :: xp      ! proton mass fraction
REAL(double), INTENT(in)     :: xhe     ! helium mass fraction
REAL(double), INTENT(in)     :: xh      ! heavy nucleus mass fraction
REAL(double), INTENT(in)     :: ah      ! heavy nucleus mass number
REAL(double), INTENT(in)     :: zh      ! heavy nucleus charge number
REAL(double), DIMENSION(nez), INTENT(in) :: egrid 
                                        ! neutrino energy grid[MeV] 

!--------------------------------------------------------------------
!        Output variables.
!--------------------------------------------------------------------

REAL(double), INTENT(out), DIMENSION(nez,2) :: cok 
                                  ! legendre coefs scattering kernel
                                  !  cok(:,1) - zeroth
                                  !  cok(:,2) - first
!--------------------------------------------------------------------
!        Local variables
!--------------------------------------------------------------------
! Equation C40 Bruenn 85
REAL(double), PARAMETER      :: g2 = ( Gw/mp**2 )**2 * hbar**5 * cvel**6

! Compute 1/(c*(hc)**3) * 4*pi*g2 in (C38) (C39), where
! 1/(c*(hc)**3) from (A19), 4*pi from integral over phase space: int 4*pi*e**2 de
REAL(double), PARAMETER      :: cc = 1.0d0 /( cvel * ( 2.d0 * pi * hbar * cvel )**3 )&
                                * ( 4.d0 * pi ) * g2

! Equation C31 Bruenn 85, proton vector coupling constant
REAL(double), PARAMETER      :: hvp =  one - cv

! Equation C32 Bruenn 85, proton axial vector coupling constant
REAL(double), PARAMETER      :: hap =  0.5d0 * ga

! Equation C33 Bruenn 85, neutron vector coupling constant
REAL(double), PARAMETER      :: hvn = -0.5d0

! Equation C34 Bruenn 85, neutron axial vector coupling constant
REAL(double), PARAMETER      :: han = -0.5d0 * ga

! Coupling constant in Equation C43 Bruenn 85
REAL(double), PARAMETER      :: cv0 =  half * ( hvp + hvn )

! Coupling constant in Equation C43 Bruenn 85
REAL(double), PARAMETER      :: cv1 =  hvp - hvn

! Bracket in Equation C38 Bruenn 85, proton zero moment coupling constant
REAL(double), PARAMETER      :: ap0 =  hvp**2 + 3.d0 * hap**2 

! Bracket in Equation C39 Bruenn 85, proton first moment coupling constant
REAL(double), PARAMETER      :: ap1 =  hvp**2 - hap**2

! Bracket in Equation C38 Bruenn 85, neutron zero moment coupling constant
REAL(double), PARAMETER      :: an0 =  hvn**2 + 3.d0 * han**2

! Bracket in Equation C39 Bruenn 85, neutron first moment coupling constant
REAL(double), PARAMETER      :: an1 =  hvn**2 - han**2

REAL(double)                 :: eeche        
                             ! helium opacity calculation
REAL(double)                 :: eche          
                             ! helium opacity calculation
REAL(double)                 :: eche2         
                             ! helium opacity calculation
REAL(double)                 :: eche3         
                             ! helium opacity calculation
REAL(double)                 :: b0he          
                             ! helium opacity calculation
REAL(double)                 :: saghe         
                             ! helium opacity calculation
REAL(double)                 :: sbghe         
                             ! helium opacity calculation

REAL(double)                 :: eec           
                             ! heavy nucleus opacity calculation
REAL(double)                 :: ec            
                             ! heavy nucleus opacity calculation
REAL(double)                 :: ec2           
                             ! heavy nucleus opacity calculation
REAL(double)                 :: ec3          
                             ! heavy nucleus opacity calculation
REAL(double)                 :: b0h           
                             ! heavy nucleus opacity calculation
REAL(double)                 :: sag           
                             ! heavy nucleus opacity calculation
REAL(double)                 :: sbg           
                             ! heavy nucleus opacity calculation
REAL(double)                 :: xnh          
                             ! heavy nucleus opacity calculation

REAL(double)                 :: aisov         
                             ! vector coupling constant
REAL(double), PARAMETER      :: aisosc = cv0  
                             ! vector coupling constant
REAL(double), PARAMETER      :: a01 = cv0 * cv0 
                             ! vector coupling constant
REAL(double)                 :: a02           
                             ! coupling constant
REAL(double)                 :: e2            
                             ! neutrino energy squared
REAL(double)                 :: etann         
                             ! neutron number corrected for blocking
REAL(double)                 :: xnn           
                             ! neutron number corrected for blocking
REAL(double)                 :: etapp         
                             ! proton number corrected for blocking
REAL(double)                 :: xnp           ! proton number corrected for blocking
REAL(double)                 :: xnhe          ! helium number
REAL(double)                 :: xheaa         ! heavy nucleus number * a_he^2
REAL(double)                 :: xhaa          ! heavy nucleus number * ah^2

REAL(double)                 :: ciicr         ! ion-ion correlation correction

REAL(double), DIMENSION(nez) :: xi_p_wm       ! weak magnetism correction for neutrino-proton scattering
REAL(double), DIMENSION(nez) :: xi_n_wm       ! weak magnetism correction for neutrino-neutron scattering
REAL(double), DIMENSION(nez) :: xib_p_wm      ! weak magnetism correction for antineutrino-proton scattering
REAL(double), DIMENSION(nez) :: xib_n_wm      ! weak magnetism correction for antineutrino-neutron scattering

REAL(double), DIMENSION(nez) :: rmdnns   ! mfp^-1 for neutrino-neutron scattering
REAL(double), DIMENSION(nez) :: rmdnns0  ! mfp^-1 for neutrino-neutron scattering
REAL(double), DIMENSION(nez) :: rmdnns1  ! mfp^-1 for neutrino-neutron scattering
REAL(double), DIMENSION(nez) :: rmdnps   ! mfp^-1 for neutrino-proton scattering
REAL(double), DIMENSION(nez) :: rmdnps0  ! mfp^-1 for neutrino-proton scattering
REAL(double), DIMENSION(nez) :: rmdnps1  ! mfp^-1 for neutrino-proton scattering
REAL(double), DIMENSION(nez) :: rmdnbns  ! mfp^-1 for antineutrino-neutron scattering
REAL(double), DIMENSION(nez) :: rmdnbns0 ! mfp^-1 for antineutrino-neutron scattering
REAL(double), DIMENSION(nez) :: rmdnbns1 ! mfp^-1 for antineutrino-neutron scattering
REAL(double), DIMENSION(nez) :: rmdnbps  ! mfp^-1 for antineutrino-proton scattering
REAL(double), DIMENSION(nez) :: rmdnbps0 ! mfp^-1 for antineutrino-proton scattering
REAL(double), DIMENSION(nez) :: rmdnbps1 ! mfp^-1 for antineutrino-proton scattering
REAL(double), DIMENSION(nez) :: rmdnhes  ! mfp^-1 for neutrino-helium scattering
REAL(double), DIMENSION(nez) :: rmdnhes0 ! mfp^-1 for neutrino-helium scattering
REAL(double), DIMENSION(nez) :: rmdnhes1 ! mfp^-1 for neutrino-helium scattering
REAL(double), DIMENSION(nez) :: rmdnhs   ! mfp^-1 for neutrino-heavy nucleus scattering
REAL(double), DIMENSION(nez) :: rmdnhs0  ! mfp^-1 for neutrino-heavy nucleus scattering
REAL(double), DIMENSION(nez) :: rmdnhs1  ! mfp^-1 for neutrino-heavy nucleus scattering

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

! Compute eta in Equation C37 Bruenn 85
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

! Compute Equation C38 for proton with coeff. : (cvel*(hc)**3) * ( C38 )
  rmdnps0(k)     = cc * xnp * ap0

! Compute Equation C38 for neutron with coeff.: (cvel*(hc)**3) * ( C38 )
  rmdnns0(k)     = cc * xnn * an0

  rmdnbps0(k)    = rmdnps0(k)
  rmdnbns0(k)    = rmdnns0(k)

! Compute Equation C44 for He with coeff.:  (cvel*(hc)**3) * ( C44 )
  rmdnhes0(k)    = cc * xheaa * a01 * saghe

! Compute Equation C44 for heavy element with coeff.: (cvel*(hc)**3) * ( C44 )
  rmdnhs0(k)     = cc * xhaa * a02 * sag

! Compute Equation C39 for proton with coeff.:  (cvel*(hc)**3) * ( C39 )
  rmdnps1(k)     = cc * xnp * ap1 / 3.0d0

! Compute Equation C39 for neutron with coeff.: (cvel*(hc)**3) * ( C39 )
  rmdnns1(k)     = cc * xnn * an1 / 3.0d0

  rmdnbps1(k)    = rmdnps1(k)
  rmdnbns1(k)    = rmdnns1(k)

! Compute Equation C45 for He with coeff.: (cvel*(hc)**3) * ( C45 )
  rmdnhes1(k)    = cc * xheaa * a01 * sbghe

! Compute Equation C45 for heavy element with coeff.: (cvel*(hc)**3)*e2*( C45 )
  rmdnhs1(k)     = cc * xhaa * a02 * sbg

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

     ! Compute Equation C44 + Equation C38 Bruenn 85
     cok(:,1) = rmdnps0  + rmdnns0  + rmdnhes0 + rmdnhs0

     ! Compute Equation C45 + Equation C39 Bruenn 85
     cok(:,2) = rmdnps1  + rmdnns1  + rmdnhes1 + rmdnhs1

  CASE (2) ! e-antineutrino

     ! Compute Equation C44 + Equation C38 Bruenn 85
     cok(:,1) = rmdnbps0 + rmdnbns0 + rmdnhes0 + rmdnhs0

     ! Compute Equation C45 + Equation C39 Bruenn 85
     cok(:,2) = rmdnbps1 + rmdnbns1 + rmdnhes1 + rmdnhs1

  CASE DEFAULT
     cok = zero
END SELECT

RETURN

END SUBROUTINE scatical_weaklib
