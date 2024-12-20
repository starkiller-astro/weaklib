SUBROUTINE scatergn_weaklib &
           ( nez, egrid, tmev, eta, h0i, h0ii, h1i, h1ii, NPS )
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
!      h0i/h1i   = ( 1/!(hc)**3 )*( g2/pi*w2*w'2 ) * houtLi(w,w')
!   
!      h1ii/h1ii = ( 1/!(hc)**3 )*( g2/pi*w2*w'2 ) * houtLii(w,w') 
!
!--------------------------------------------------------------------
!
!    The Legendre moments of the neutrino-electron scattering functions
!    are given by
!
!      phiLout = ( 1/(hc)**3 )*( g2/pi*w2*w'2 ) 
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

USE wlKindModule, ONLY: dp
USE numerical_module, ONLY: one, zero, pi
USE physcnst_module, ONLY: Gw, mp, hbar, cvel, kmev, cv, ca
USE pair_module, ONLY: coef

IMPLICIT NONE

!-------------------------------------------------------------------
!        Input variables
!--------------------------------------------------------------------

INTEGER,      INTENT(in) :: nez    ! number of neutrino energy group
REAL(dp), INTENT(in) :: tmev   ! temperature [MeV]
REAL(dp), INTENT(in) :: eta    ! electron chemical potential/TMeV
REAL(dp), DIMENSION(nez), INTENT(in) :: egrid   ! neutrino energy grid
INTEGER,  INTENT(IN) :: NPS !add neutrino positron scattering

!--------------------------------------------------------------------
!        Output variables
!-------------------------------------------------------------------

REAL(dp), DIMENSION(nez,nez), INTENT(out) :: h0i
REAL(dp), DIMENSION(nez,nez), INTENT(out) :: h0ii
REAL(dp), DIMENSION(nez,nez), INTENT(out) :: h1i
REAL(dp), DIMENSION(nez,nez), INTENT(out) :: h1ii

!--------------------------------------------------------------------
!        Local variables
!--------------------------------------------------------------------

INTEGER                     :: i             ! loop counter over points to build
INTEGER                     :: k             ! incomiong neutrino energy zone index
INTEGER                     :: kp            ! outgoing neutrino energy zone index

REAL(dp), DIMENSION(nez) :: wk2   ! egrid**2
REAL(dp)                :: enuin         ! incoming neutrino energy/tmev
REAL(dp)                :: enuout        ! outgoing neutrino energy/tmev

REAL(dp)                :: cxct          ! coef*t**6
REAL(dp)                :: cxc           ! cxc/egrid**2
REAL(dp)                 :: fexp          ! exponential function
REAL(dp)                :: hout0i        ! zero moment of outgoing scattering function type i
REAL(dp)                :: hout0ii       ! zero moment of outgoing scattering function type ii
REAL(dp)                :: hout1i        ! first moment of outgoing scattering function type i
REAL(dp)                :: hout1ii       ! first moment of outgoing scattering function type ii

REAL(dp), DIMENSION(nez,nez) :: scate_0i  ! zero moment of electron scatering, type i
REAL(dp), DIMENSION(nez,nez) :: scate_0ii ! zero moment of electron scatering, type ii
REAL(dp), DIMENSION(nez,nez) :: scatp_0i  ! zero moment of positron scatering, type i
REAL(dp), DIMENSION(nez,nez) :: scatp_0ii ! zero moment of positron scatering, type ii
REAL(dp), DIMENSION(nez,nez) :: scate_1i  ! first moment of electron scatering, type i
REAL(dp), DIMENSION(nez,nez) :: scate_1ii ! first moment of electron scatering, type ii
REAL(dp), DIMENSION(nez,nez) :: scatp_1i  ! first moment of positron scatering, type i
REAL(dp), DIMENSION(nez,nez) :: scatp_1ii ! first moment of positron scatering, type ii

REAL(dp) :: eta_p !-eta

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

  cxct             = coef * (tmev)**6

!--------------------------------------------------------------------
!  cxc has dimensions of 1 /[ energy length ]
!  cxc/wje2 has dimensions of 1 /[ energy**3 length ]
!--------------------------------------------------------------------

    DO k = 1,nez

      cxc            = cxct/wk2(k)
      enuin          = egrid(k)/tmev

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

    IF (NPS > 0) THEN

      eta_p = -eta

      DO k = 1,nez

        cxc            = cxct/wk2(k)
        enuin          = egrid(k)/tmev

        DO kp = 1,k

          enuout       = egrid(kp)/tmev

          CALL sctlgndv_weaklib &
               ( enuin, enuout, eta_p, hout0i, hout0ii, hout1i, hout1ii )

          scatp_0i (k,kp) = hout0i * cxc/wk2(kp) 
          scatp_0ii(k,kp) = hout0ii* cxc/wk2(kp)
          scatp_0i (kp,k) = scatp_0i (k,kp) * fexp( (egrid(kp) - egrid(k) ) / TMeV )
          scatp_0ii(kp,k) = scatp_0ii(k,kp) * fexp( (egrid(kp) - egrid(k) ) / TMeV )

          scatp_1i (k,kp) = hout1i * cxc/wk2(kp)
          scatp_1ii(k,kp) = hout1ii* cxc/wk2(kp)
          scatp_1i (kp,k) = scatp_1i (k,kp) * fexp( (egrid(kp) - egrid(k) ) / TMeV )
          scatp_1ii(kp,k) = scatp_1ii(k,kp) * fexp( (egrid(kp) - egrid(k) ) / TMeV )


        END DO ! kp = 1,k
      END DO ! k = 1,nez

    ENDIF

    h0i  = scate_0i  + scatp_0ii
    h0ii = scate_0ii + scatp_0i
    h1i  = scate_1i  + scatp_1ii
    h1ii = scate_1ii + scatp_1i

RETURN

END SUBROUTINE scatergn_weaklib
