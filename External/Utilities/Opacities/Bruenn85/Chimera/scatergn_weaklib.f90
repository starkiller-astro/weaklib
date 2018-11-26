SUBROUTINE scatergn_weaklib( n, im, nez, egrid, tmev, cmpe, NES )
!-----------------------------------------------------------------------
!
!    Author:       R. Chu, Dept. Phys. & Astronomy
!                  U. Tennesee, Knoxville
!
!    Created:      11/2/18
!
!    Purpose:
!      To compute the zero and first legendre coefs for the n-type neutrino-electron 
!      scattering functions.
!
!-----------------------------------------------------------------------
!
!      The zero and first legendre coefs for the n-type neutrino scattering
!       functions are computed here These are included in the multi-group
!       diffusion equations, which have the form
!
!            pis1 = -yt*( dpsi0/dr - a1w*psi0 - c1w)
!
!            dpsi0/dt/! + d(r2*psi1)/dr/3./r2
!                 = xw + yw*psi0 + zw*dpsi0/dr
!
!       where
!
!            1/yt = zltot = jw + yaw - b1w
!            xw   = jw + c0w + b0w*c1w*yt
!            yw   = -jw - yaw + a0w + b0w*a1w*yt
!            zw   = -b0w*yt
!
!      Neutrino scattering contributes to the terms a0w,a1w, b0w, b1w, c0w, c1w as follows:
!
!            a0w  = -K   Int{w2'dw'[ sctin0(w,w')psi0(w') + sctot0(w,w')( 1 - psi0(w') ) ] }
!            b0w  = -K/3 Int{w2'dw'[ sctin1(w,w') - sctot1(w,w') ]psi1(w') }
!            c0w  =  K   Int{w2'dw'[ sctin0(w,w')psi0(w') ] }
!            a1w  = -K   Int{w2'dw'[ sctin1(w,w') - sctot1(w,w') ]psi1(w') }
!            b1w  = -K   Int{w2'dw'[ sctin0(w,w')psi0(w') + sctot0(w,w')( 1 - psi0(w') ) ] }
!            c1w  =  K   Int{w2'dw'[ sctin1(w,w')psi1(w') ] }
!
!       where
!
!            K    = 2*pi/! * 1/(hc)**3
!
!       for neutrino-electron scattering.
!
!      The Legendre moments of the neutrino-electron scattering functions
!       are given by
!
!          phiLin = ( 2*pi/!(hc)**3 )*( g2/pi*w2*w'2 ) * cnes1(n)*hinLi(w,w') + cnes2(n)*hinLii(w,w')
!
!       where!
!
!     hinLi(w,w') = Int{ de*Fe(e)*[1 - F(e + w - w')] * exp[(w - w')/kT]*hLi(e,w,w') }
!
!                 = (kT)**6 * Int{ d(e/kT) * Fe(e/kT)*[1 - F(e/kT + w/kT - w'/kT)]
!                  * exp[(w - w')/kT]*hLi(e/kT,w/kT,w'/kT) }
!
!      and hinLii(w,w') is defined likewise.
!
!         phiLout = ( 2*pi/!(hc)**3 )*( g2/pi*w2*w'2 ) * cnes1(n)*houtLi(w,w') + cnes2(n)*houtLii(w,w')
!
!      where
!
!    houtLi(w,w') = Int{ de*Fe(e)*[1 - F(e + w - w')]
!                  * hLo(e,w,w') }
!
!                 = (kT)**6 * Int{ d(e/kT) * Fe(e/kT)*[1 - F(e/kT + w/kT - w'/kT)] * hLo(e/kT,w/kT,w'/kT) }
!
!      and houtLii(w,w') is defined likewise.
!
!     The integrations are performed in subroutine sctgldnv.
!
!-----------------------------------------------------------------------
!
!    Subprograms called:
!
!     sctlgndv   : executes the computation of the neutrino-electron
!                  scattering rates
!     terminate  : terminates the run in the event of an error
!
!
!-----------------------------------------------------------------------
!        Units
!
!  Gw         : G_fermi * ( mpc^{2} )^{2} dimensionless
!  Gw/mp**2   : MeV^{-2}
!  coc        : MeV^{-5} cm^{-1}
!  cxct       : MeV cm^{-1}
!  cxc        : MeV^{-1} cm^{-1}
!  cxc/egrid**2: MeV^{-3} cm^{-1}
!
!-----------------------------------------------------------------------

USE kind_module, ONLY: double
USE numerical_module, ONLY: one, zero, pi
USE physcnst_module, ONLY: Gw, mp, hbar, cvel, kmev, cv, ca

USE vector_functions_module, ONLY: abflog_vec

IMPLICIT NONE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER,      INTENT(in) ::  n     ! neutrino flavor index 
INTEGER,      INTENT(in) :: im     ! moment order index 
INTEGER,      INTENT(in) :: nez    ! number of neutrino energy group
REAL(double), INTENT(in) :: cmpe   ! electron chemical potential
REAL(double), INTENT(in) :: tmev   ! temperature [MeV]
REAL(double), DIMENSION(nez), INTENT(in) :: egrid   ! neutrino energy grid

!--------------------------------------------------------------------
!        Output variables
!-------------------------------------------------------------------

REAL(double), DIMENSION(nez,nez), INTENT(out) :: NES 

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                     :: i             ! loop counter over points to build
INTEGER                     :: k             ! incomiong neutrino energy zone index
INTEGER                     :: kp            ! outgoing neutrino energy zone index

REAL(double)                :: eta           ! cmpe/tmev
REAL(double), DIMENSION(nez) :: wk2   ! egrid**2
REAL(double)                :: enuin         ! incoming neutrino energy/tmev
REAL(double)                :: enuout        ! outgoing neutrino energy/tmev

REAL(double), PARAMETER     :: coc = 2.d0 * ( Gw/mp**2 )**2 * 1.d0/( 8.d0 * pi**3 * hbar * cvel )
REAL(double)                :: cxct          ! coc*t**6
REAL(double)                :: cxc           ! cxc/egrid**2

REAL(double)                :: hin0i         ! zero moment of incoming scattering function type i
REAL(double)                :: hin0ii        ! zero moment of incoming scattering function type ii
REAL(double)                :: hout0i        ! zero moment of outgoing scattering function type i
REAL(double)                :: hout0ii       ! zero moment of outgoing scattering function type ii
REAL(double)                :: hin1i         ! first moment of incoming scattering function type i
REAL(double)                :: hin1ii        ! first moment of incoming scattering function type ii
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

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

wk2(:)              = egrid(:) * egrid(:)
scate_0i (:,:)       = zero
scate_0ii(:,:)       = zero
scatp_0i (:,:)       = zero
scatp_0ii(:,:)       = zero
scate_1i (:,:)       = zero
scate_1ii(:,:)       = zero
scatp_1i (:,:)       = zero
scatp_1ii(:,:)       = zero

!-----------------------------------------------------------------------
!  Compute equation of state quantities needed for the computation
!   of the neutrino scattering functions.
!
!  cxct has dimensions of [ energy / length ].
!-----------------------------------------------------------------------

  eta              = cmpe/tmev
  cxct             = coc * (tmev)**6

!SELECT CASE (n)
!  CASE (1)  ! electron-type neutrino

!-----------------------------------------------------------------------
!  Neutrino-electron down and iso-e scattering kernals
!
!  cxc has dimensions of 1 /[ energy length ]
!  cxc/wje2 has dimensions of 1 /[ energy**3 length ]
!-----------------------------------------------------------------------

    DO k = 1,nez

      cxc            = cxct/wk2(k)
      enuin          = egrid(k)/tmev

      DO kp = 1,k

        enuout       = egrid(kp)/tmev
        CALL sctlgndv( enuin, enuout, eta, hin0i, hin0ii, hout0i, hout0ii, &
&        hin1i, hin1ii, hout1i, hout1ii )

!-----------------------------------------------------------------------
!  Load scattering kernal arrays.
!
!  Down or isoenergetic scattering only is computed. Upscattering is
!   computed using detailed balance in scterate.f90. Thus:
!  scte0i(k,kp,icube) and scte0ii(k,kp,icube) 
!   are the out-scattering functions for a neutrino of incident in-beam
!   energy k to a final out-beam energy kp <= k.
!  scte0i(kp,k,icube) and scte0ii(kp,k,icube)
!   are the in-scattering function for a neutrino of incident out-beam
!   energy k to a final in-beam energy kp <= k.
!  Since the above in-beam and out-beam scatterings are the same scattering
!   and the corresponding scattering functions are therefore equal.
!  Consequently, scte0i(k,kp,icube) and scte0ii(k,kp,icube) for kp <= k 
!   represents "out" scattering, and for
!   kp > k represents "in" scattering. This difference is taken into
!   account in scterate.f90.
!-----------------------------------------------------------------------

        scate_0i (k,kp) = hout0i * cxc/wk2(kp)
        scate_0ii(k,kp) = hout0ii* cxc/wk2(kp)
        scate_0i (kp,k) = scate_0i (k,kp) 
        scate_0ii(kp,k) = scate_0ii(k,kp)
!------------?
        scate_1i (k,kp) = hout1i * cxc/wk2(kp)
        scate_1ii(k,kp) = hout1ii* cxc/wk2(kp)
        scate_1i (kp,k) = scate_1i (k,kp)
        scate_1ii(kp,k) = scate_1ii(k,kp)


      END DO ! kp = 1,k
    END DO ! k = 1,nez

  SELECT CASE (n)
    CASE (1)  
    SELECT CASE (im)
      CASE (0)
        NES = scate_0i *(cv+ca)**2 + scate_0ii * (ca-ca)**2
      CASE(1)
        NES = scate_1i *(cv+ca)**2 + scate_1ii * (ca-ca)**2
      CASE DEFAULT
        NES = zero
    END SELECT
    CASE (2)
    SELECT CASE (im)
      CASE (0)
        NES = scate_0i *(cv-ca)**2 + scate_0ii * (ca+ca)**2
      CASE(1)
        NES = scate_1i *(cv-ca)**2 + scate_1ii * (ca+ca)**2
      CASE DEFAULT
        NES = zero
    END SELECT
  END SELECT

!  CASE (2) ! electron-type antineutrino
!
!    eta              = -eta
!
!    DO k = 1,nez
!
!      cxc            = cxct/wk2(k)
!      enuin          = egrid(k)/tmev
!
!      DO kp = 1,k
!
!        enuout       = egrid(kp)/tmev
!        CALL sctlgndv( enuin, enuout, eta, hin0i, hin0ii, hout0i, hout0ii, &
!&        hin1i, hin1ii, hout1i, hout1ii )
!
!!-----------------------------------------------------------------------
!!  Load scattering kernal arrays.
!!
!!  Down or isoenergetic scattering only is computed. Upscattering is
!!   computed using detailed balance in scterate.f90. Thus:
!!  scte0i(k,kp,icube) and scte0ii(k,kp,icube) 
!!   are the out-scattering functions for a neutrino of incident in-beam
!!   energy k to a final out-beam energy kp <= k.
!!  scte0i(kp,k,icube) and scte0ii(kp,k,icube)
!!   are the in-scattering function for a neutrino of incident out-beam
!!   energy k to a final in-beam energy kp <= k.
!!  Since the above in-beam and out-beam scatterings are the same scattering
!!   and the corresponding scattering functions are therefore equal.
!!  Consequently, scte0i(k,kp,icube) and scte0ii(k,kp,icube) for kp <= k 
!!   represents "out" scattering, and for
!!   kp > k represents "in" scattering. This difference is taken into
!!   account in scterate.f90.
!!-----------------------------------------------------------------------
!
!        scatp_0i (k,kp) = hout0i * cxc/wk2(kp)
!        scatp_0ii(k,kp) = hout0ii* cxc/wk2(kp)
!        scatp_0i (kp,k) = scatp_0i (k,kp)
!        scatp_0ii(kp,k) = scatp_0ii(k,kp)
!!------------?
!        scatp_1i (k,kp) = hout1i * cxc/wk2(kp)
!        scatp_1ii(k,kp) = hout1ii* cxc/wk2(kp)
!        scatp_1i (kp,k) = scatp_1i (k,kp)
!        scatp_1ii(kp,k) = scatp_1ii(k,kp)
!
!      END DO ! kp = 1,k
!    END DO ! k = 1,nez
!
!    SELECT CASE (im)
!      CASE (0)
!        NES = scatp_0i *(cv+ca)**2 + scatp_0ii * (ca-ca)**2
!      CASE(1)
!        NES = scatp_1i *(cv+ca)**2 + scatp_1ii * (ca-ca)**2
!      CASE DEFAULT
!        NES = zero
!    END SELECT
!
!  CASE DEFAULT 
!
!    NES = zero
!
!END SELECT

RETURN
END SUBROUTINE scatergn_weaklib
