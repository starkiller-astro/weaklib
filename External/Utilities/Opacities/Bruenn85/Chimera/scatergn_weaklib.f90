SUBROUTINE scatergn_weaklib &
           ( nez, egrid, tmev, eta, h0i, h0ii, h1i, h1ii )
!--------------------------------------------------------------------
!
!    Author:       R. Chu, Dept. Phys. & Astronomy
!                  U. Tennesee, Knoxville
!
!    Created:      12/04/2018
!
!    Purpose:
!      To compute the neutrino-electron scattering function elements
!      h0i/h1i and h1i/h1ii (zeroth and first order), which have the 
!      forms
!
!      h0i/h1i   = ( 2*pi/!(hc)**3 )*( g2/pi*w2*w'2 ) * houtLi(w,w')
!   
!      h1ii/h1ii = ( 2*pi/!(hc)**3 )*( g2/pi*w2*w'2 ) * houtLii(w,w') 
!
!--------------------------------------------------------------------
!
!    The Legendre moments of the neutrino-electron scattering functions
!    are given by
!
!      phiLout = ( 2*pi/!(hc)**3 )*( g2/pi*w2*w'2 ) 
!                  * cnes1(n)*houtLi(w,w') + cnes2(n)*houtLii(w,w')
!
!    where
!
!      houtLi(w,w') = Int{ de * Fe(e) * [1 - F(e + w - w')]
!                                     * hLo(e,w,w') }
!
!                   = (kT)**6 
!                   * Int{ d(e/kT) * Fe(e/kT) 
!                          * [1 - F(e/kT + w/kT - w'/kT)] 
!                          * hLo(e/kT,w/kT,w'/kT) }
!
!    and houtLii(w,w') is defined likewise.
! 
!    The integrations are performed in subroutine sctgldnv_weaklib.
!
!-------------------------------------------------------------------

USE kind_module, ONLY: double
USE numerical_module, ONLY: one, zero, pi
USE physcnst_module, ONLY: Gw, mp, hbar, cvel, kmev, cv, ca

IMPLICIT NONE

!-------------------------------------------------------------------
!        Input variables
!--------------------------------------------------------------------

INTEGER,      INTENT(in) :: nez    ! number of neutrino energy group
REAL(double), INTENT(in) :: tmev   ! temperature [MeV]
REAL(double), INTENT(in) :: eta    ! electron chemical potential/TMeV
REAL(double), DIMENSION(nez), INTENT(in) :: egrid   ! neutrino energy grid

!--------------------------------------------------------------------
!        Output variables
!-------------------------------------------------------------------

REAL(double), DIMENSION(nez,nez), INTENT(out) :: h0i
REAL(double), DIMENSION(nez,nez), INTENT(out) :: h0ii
REAL(double), DIMENSION(nez,nez), INTENT(out) :: h1i
REAL(double), DIMENSION(nez,nez), INTENT(out) :: h1ii

!--------------------------------------------------------------------
!        Local variables
!--------------------------------------------------------------------

INTEGER                     :: i             ! loop counter over points to build
INTEGER                     :: k             ! incomiong neutrino energy zone index
INTEGER                     :: kp            ! outgoing neutrino energy zone index

REAL(double), DIMENSION(nez) :: wk2   ! egrid**2
REAL(double)                :: enuin         ! incoming neutrino energy/tmev
REAL(double)                :: enuout        ! outgoing neutrino energy/tmev

REAL(double), PARAMETER     :: coc = 2.d0 * ( Gw/mp**2 )**2 * 1.d0/( 8.d0 * pi**3 * hbar * cvel )
REAL(double)                :: cxct          ! coc*t**6
REAL(double)                :: cxc           ! cxc/egrid**2
REAL(double)                 :: fexp          ! exponential function
REAL(double)                :: hout0i        ! zero moment of outgoing scattering function type i
REAL(double)                :: hout0ii       ! zero moment of outgoing scattering function type ii
REAL(double)                :: hout1i        ! first moment of outgoing scattering function type i
REAL(double)                :: hout1ii       ! first moment of outgoing scattering function type ii

REAL(double), DIMENSION(nez,nez) :: scate_0i  ! zero moment of electron scatering, type i
REAL(double), DIMENSION(nez,nez) :: scate_0ii ! zero moment of electron scatering, type ii
REAL(double), DIMENSION(nez,nez) :: scatp_0i  ! zero moment of positron scatering, type i
REAL(double), DIMENSION(nez,nez) :: scatp_0ii ! zero moment of positron scatering, type ii
REAL(double), DIMENSION(nez,nez) :: scate_1i  ! first moment of electron scatering, type i
REAL(double), DIMENSION(nez,nez) :: scate_1ii ! first moment of electron scatering, type ii
REAL(double), DIMENSION(nez,nez) :: scatp_1i  ! first moment of positron scatering, type i
REAL(double), DIMENSION(nez,nez) :: scatp_1ii ! first moment of positron scatering, type ii

EXTERNAL fexp
!--------------------------------------------------------------------
!--------------------------------------------------------------------

wk2(:)              = egrid(:) * egrid(:)
scate_0i (:,:)       = zero
scate_0ii(:,:)       = zero
scatp_0i (:,:)       = zero
scatp_0ii(:,:)       = zero
scate_1i (:,:)       = zero
scate_1ii(:,:)       = zero
scatp_1i (:,:)       = zero
scatp_1ii(:,:)       = zero

!--------------------------------------------------------------------
!  cxct has dimensions of [ energy / length ].
!--------------------------------------------------------------------

  cxct             = coc * (tmev)**6
!!!???????????

!--------------------------------------------------------------------
!  cxc has dimensions of 1 /[ energy length ]
!  cxc/wje2 has dimensions of 1 /[ energy**3 length ]
!--------------------------------------------------------------------

    DO k = 1,nez

      cxc            = cxct/wk2(k)
      enuin          = egrid(k)/tmev

!     DO kp = 1,nez !!! fix me
      DO kp = 1,k

        enuout       = egrid(kp)/tmev

        CALL sctlgndv_weaklib &
             ( enuin, enuout, eta, hout0i, hout0ii, hout1i, hout1ii )

        scate_0i (k,kp) = hout0i * cxc/wk2(kp) 
        scate_0ii(k,kp) = hout0ii* cxc/wk2(kp)
        scate_0i (kp,k) = scate_0i (k,kp) * fexp( (egrid(kp) - egrid(k) ) / TMeV )
        scate_0ii(kp,k) = scate_0ii(k,kp) * fexp( (egrid(kp) - egrid(k) ) / TMeV )

        scate_1i (k,kp) = hout1i * cxc/wk2(kp)
        scate_1ii(k,kp) = hout1ii* cxc/wk2(kp)
        scate_1i (kp,k) = scate_1i (k,kp) * fexp( (egrid(kp) - egrid(k) ) / TMeV )
        scate_1ii(kp,k) = scate_1ii(k,kp) * fexp( (egrid(kp) - egrid(k) ) / TMeV )


      END DO ! kp = 1,k
    END DO ! k = 1,nez

    h0i  = scate_0i
    h0ii = scate_0ii
    h1i  = scate_1i
    h1ii = scate_1ii

RETURN

END SUBROUTINE scatergn_weaklib
