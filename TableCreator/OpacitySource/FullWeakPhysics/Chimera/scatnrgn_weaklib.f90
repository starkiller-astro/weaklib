SUBROUTINE scatnrgn_weaklib &
           ( nez, egrid_c, egrid_e, tmev, chem_n, chem_p, wm, g_strange, &
             phi0_nu_n, phi1_nu_n, phi0_nub_n, phi1_nub_n, &
             phi0_nu_p, phi1_nu_p, phi0_nub_p, phi1_nub_p )

!-----------------------------------------------------------------------
!
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/09/03
!
!    Modified:     C. Y. Cardall for weaklib
!
!    Purpose:
!      To set up the data for the computation of the neutrino-
!       nucleon elastic scattering functions, to call subroutine
!       scatncal for the computation of the scattering functions,
!       and to load the computed scattering functions in arrays
!       for interpolation.
!
!-----------------------------------------------------------------------
!
!      The zero and first legendre coefs for the n-type neutrino scattering
!       functions are computed here
!      These are included in the multi-group diffusion equations, which
!       have the form
!
!            pis1 = -yt*( dpsi0/dr - a1w*psi0 - c1w)
!
!            dpsi0/dt/! + d(r2*psi1)/dr/3./r2 = xw + yw*psi0 + zw*dpsi0/dr
!
!       where
!
!            1/yt = zltot = jw + yaw - b1w
!            xw   = jw + c0w + b0w*c1w*yt
!            yw   = -jw - yaw + a0w + b0w*a1w*yt
!            zw   = -b0w*yt
!
!      Neutrino scattering contributes to the terms a0w,a1w, b0w, b1w,
!       c0w, c1w as follows:
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
!       where
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
!    houtLi(w,w') = Int{ de*Fe(e)*[1 - F(e + w - w')] * hLo(e,w,w') }
!
!                 = (kT)**6 * Int{ d(e/kT)* Fe(e/kT)*[1 - F(e/kT + w/kT - w'/kT)]
!                  * hLo(e/kT,w/kT,w'/kT) }
!
!      and houtLii(w,w') is defined likewise.
!
!     The integrations are performed in subroutine sctgldnv.
!
!-----------------------------------------------------------------------
!
!    Subprograms called:
!  scatncal   : executes the computation of the neutrino-
!                nucleon elastic scattering rates
!
!    Input arguments:
!   nez       : size of energy grid
!   egrid_c   : neutrino energy grid, center values
!   egrid_e   : neutrino energy grid, edge values
!   tmev      : temperature, MeV
!   chem_n    : neutron chemical potential, Chimera
!   chem_p    : proton chemical potential, Chimera
!   wm        : include weak magnetism corrections
!   g_strange : strange quark contributions
!
!    Output arguments (without many-body corrections):
!   phi0_nu_n  : 0th Legendre coefficient for nu scattering from neutrons
!   phi1_nu_n  : 1st Legendre coefficient for nu scattering from neutrons
!   phi0_nub_n : 0th Legendre coefficient for nuBar scattering from neutrons
!   phi1_nub_n : 1st Legendre coefficient for nuBar scattering from neutrons
!   phi0_nu_p  : 0th Legendre coefficient for nu scattering from protons
!   phi1_nu_p  : 1st Legendre coefficient for nu scattering from protons
!   phi0_nub_p : 0th Legendre coefficient for nuBar scattering from protons
!   phi1_nub_p : 1st Legendre coefficient for nuBar scattering from protons
!
!-----------------------------------------------------------------------

USE wlKindModule, ONLY: dp
USE numerical_module, ONLY: half, one, epsilon
USE physcnst_module, ONLY: kmev, dmnp, mn, mp, ga, sin2W, kfm

IMPLICIT NONE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER,  INTENT(in)     :: nez    ! number of neutrino energy group
REAL(dp), DIMENSION(nez),   INTENT(in) :: egrid_c  ! neutrino energy grid
REAL(dp), DIMENSION(nez+1), INTENT(in) :: egrid_e  ! neutrino energy grid
REAL(dp), INTENT(in)     :: tmev   ! temperature [MeV]
REAL(dp), INTENT(in)     :: chem_n ! neutron chemical potential, Chimera
REAL(dp), INTENT(in)     :: chem_p ! proton chemical potential, Chimera
INTEGER,  INTENT(in)     :: wm
REAL(dp), INTENT(in)     :: g_strange

!--------------------------------------------------------------------
!        Output variables
!-------------------------------------------------------------------

REAL(dp), DIMENSION(nez,nez), INTENT(out) :: phi0_nu_n
REAL(dp), DIMENSION(nez,nez), INTENT(out) :: phi1_nu_n
REAL(dp), DIMENSION(nez,nez), INTENT(out) :: phi0_nub_n
REAL(dp), DIMENSION(nez,nez), INTENT(out) :: phi1_nub_n
REAL(dp), DIMENSION(nez,nez), INTENT(out) :: phi0_nu_p
REAL(dp), DIMENSION(nez,nez), INTENT(out) :: phi1_nu_p
REAL(dp), DIMENSION(nez,nez), INTENT(out) :: phi0_nub_p
REAL(dp), DIMENSION(nez,nez), INTENT(out) :: phi1_nub_p

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

LOGICAL, SAVE                :: first = .true.

! INTEGER                      :: i             ! loop counter over points to build
! INTEGER                      :: icube         ! index pointer to currently built cube

INTEGER                      :: k        ! incoming neutrino energy zone index
INTEGER                      :: kp       ! outgoing neutrino energy zone index

! REAL(dp), DIMENSION(igen)    :: rho          ! density [g cm^{-3}]
! REAL(dp), DIMENSION(igen)    :: temp         ! temperature [K]
! REAL(dp), DIMENSION(igen)    :: ye           ! electron fraction
! REAL(dp), DIMENSION(igen)    :: chem_n       ! neutron chemical potential
! REAL(dp), DIMENSION(igen)    :: chem_p       ! proton chemical potential
! REAL(dp), DIMENSION(igen)    :: Interpolants ! EoS table call output

! REAL(dp)                 :: tmev          ! temperature [MeV]
REAL(dp)                 :: sm            ! target particle rest mass
! REAL(dp)                 :: brydns        ! baryon density

REAL(dp)                 :: e_in    ! incoming neutrino zone-centered energy
REAL(dp)                 :: e_in_l  ! incoming neutrino lower zone-edged energy
REAL(dp)                 :: e_in_u  ! incoming neutrino upper zone-edged energy
REAL(dp)                 :: e_out   ! outgoing neutrino zone-centered energy
REAL(dp)                 :: e_out_l ! outgoing neutrino lower zone-edged energy
REAL(dp)                 :: e_out_u ! outgoing neutrino upper zone-edged energy

REAL(dp)                 :: sct0_n  ! zero momemt of the neutrino-neutron sct function
REAL(dp)                 :: sct1_n  ! first momemt of the neutrino-neutron sct function
REAL(dp)                 :: sct0_p  ! zero momemt of the neutrino-proton sct function
REAL(dp)                 :: sct1_p  ! first momemt of the neutrino-neutron sct function

REAL(dp), DIMENSION(nez)   :: xi_p_wm    ! weak magnetism correction for neutrino-proton scattering
REAL(dp), DIMENSION(nez)   :: xi_n_wm    ! weak magnetism correction for neutrino-neutron scattering
REAL(dp), DIMENSION(nez)   :: xib_p_wm   ! weak magnetism correction for antineutrino-proton scattering
REAL(dp), DIMENSION(nez)   :: xib_n_wm   ! weak magnetism correction for antineutrino-neutron scattering
! REAL(dp)                          :: S_tot      ! many-body modification to nu_N scattering

! REAL(dp), DIMENSION(nez+nezext)   :: unu_local  ! local copy of energy group centers
! REAL(dp), DIMENSION(nez+nezext+1) :: unub_local ! local copy of energy group edges
! REAL(dp), DIMENSION(nez+nezext,nez+nezext) :: sctmp 
! REAL(dp), DIMENSION(nez+nezext,nez+nezext) :: sctmpb

! CHARACTER(len=1)             :: flag          ! nuclear eos selection flag

! LOGICAL                      :: fail          ! did EoS fail to converge

! REAL(dp)                 :: press         ! pressure
! REAL(dp)                 :: energ         ! internal energy
! REAL(dp)                 :: entrop        ! entropy [kb/baryon]
REAL(dp)                 :: mu_n          ! free neutron chemical potential
REAL(dp)                 :: mu_p          ! free proton chemical potential
! REAL(dp)                 :: chem_e        ! electron chemical potential
! REAL(dp)                 :: xn_neut       ! free neutron fraction
! REAL(dp)                 :: xn_prot       ! free proton fraction
! REAL(dp)                 :: xn_heavy      ! free proton fraction
! REAL(dp)                 :: a_heavy       ! A for mean heavy nucleus
! REAL(dp)                 :: z_heavy       ! Z for mean heavy nucleus
! REAL(dp)                 :: be_heavy      ! Binding energy for mean heavy nucleus

REAL(dp), SAVE           :: cv_p_nscat    ! vector current neutrino-proton coupling constant
REAL(dp), SAVE           :: cv_n_nscat    ! vector current neutrino-neutron coupling constant
REAL(dp), SAVE           :: ca_p_nscat    ! axial vector current neutrino-proton coupling constant
REAL(dp), SAVE           :: ca_n_nscat    ! axial vector current neutrino-neutron coupling constant

REAL(dp)                 :: fexp          ! exponential function

EXTERNAL fexp

!-----------------------------------------------------------------------
!  Compute coupling constants
!-----------------------------------------------------------------------

IF (first) THEN
  cv_p_nscat       = 0.5d0 - 2.d0 * sin2W
  cv_n_nscat       = - 0.5d0
  ca_p_nscat       = (   ga - g_strange )/2.d0
  ca_n_nscat       = ( - ga - g_strange )/2.d0
  first            = .false.
END IF ! first

!-----------------------------------------------------------------------
!  Absolute chemical potentials
!-----------------------------------------------------------------------

  mu_p             = chem_p + dmnp + mp
  mu_n             = chem_n + dmnp + mn

!-----------------------------------------------------------------------
!  Weak magnetism corrections for neutrino and antineutrino neutron and
!   proton scattering.
!-----------------------------------------------------------------------

  !Make weak magnetism corrections non-operational initially
  !(these are multiplicative corrections)
  xi_p_wm  = 1.0d0
  xi_n_wm  = 1.0d0
  xib_p_wm = 1.0d0
  xib_n_wm = 1.0d0

  IF(wm == 1) THEN
    CALL nc_weak_mag_weaklib &
           ( egrid_c, xi_p_wm, xi_n_wm, xib_p_wm, xib_n_wm, nez )
  END IF

! !-----------------------------------------------------------------------
! !  Many body corrections for neutrino-nucleon neutral current
! !   scattering.
! !-----------------------------------------------------------------------

!   tmev             = kmev * temp(i)
!   brydns           = rho(i) * kfm

!   CALL nc_manybody_mod( brydns, tmev, ye(i), S_tot )

!-----------------------------------------------------------------------
!
!        \\\\\ COMPUTE SCATTERING KERNELS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Neutrino-nucleon downscattering kernels
!  Loop over incoming energy
!-----------------------------------------------------------------------

  DO k = nez,1,-1
    e_in           = egrid_c(k)
    e_in_l         = egrid_e(k)
    e_in_u         = egrid_e(k+1)

!-----------------------------------------------------------------------
!  Neutrino-nucleon downscattering kernels
!  Loop over outgoing energy
!-----------------------------------------------------------------------

    DO kp = 1,k
      e_out        = egrid_c(kp)
      e_out_l      = egrid_e(kp)
      e_out_u      = egrid_e(kp+1)

!-----------------------------------------------------------------------
!  Compute neutrino-neutron inelastic scattering kernels
!-----------------------------------------------------------------------

      sm           = mn
    
      CALL scatncal_weaklib &
             ( e_in, e_out, e_in_l, e_in_u, e_out_l, e_out_u, tmev, &
               mu_n, sm, sct0_n, sct1_n, cv_n_nscat, ca_n_nscat )

!-----------------------------------------------------------------------
!  Compute neutrino-proton inelastic scattering kernels
!-----------------------------------------------------------------------

      sm             = mp

      CALL scatncal_weaklib &
             ( e_in, e_out, e_in_l, e_in_u, e_out_l, e_out_u, tmev, &
               mu_p, sm, sct0_p, sct1_p, cv_p_nscat, ca_p_nscat )

!-----------------------------------------------------------------------
!  Load scattering kernal arrays.
!
!  Down or isoenergetic scattering only is computed. Upscattering is
!   computed using detailed balance in scatnrate. Thus:
!  sctn0(k,kp,icube) is the out-scattering function for a neutrino of
!   incident in-beam energy k to a final out-beam energy kp <= k.
!  sctn0(kp,k,icube) is the in-scattering function for a neutrino of
!   incident out-beam energy k to a final in-beam energy kp <= k. 
!  Since the above in-beam and out-beam scatterings are the same scattering
!   and the corresponding scattering functions are therefore equal.
!  Consequently, sctn0(k,kp,icube) for kp <= k represents "out"
!   scattering, and for kp > k represents "in" scattering. This 
!   difference is taken into account in scatnrate.
!  Likewise for the antineutrino scattering functions sctnb0(k,kp,icube)
!   and sctnb0(kp,k,icube).
!-----------------------------------------------------------------------

!       sctmp (k,kp) = DABS( sct0_n * xi_n_wm(k) * S_tot  + sct0_p * xi_p_wm(k) * S_tot )
!       sctmp (kp,k) = sctmp (k,kp)
!       sctmpb(k,kp) = DABS( sct0_n * xib_n_wm(k) * S_tot + sct0_p * xib_p_wm(k) * S_tot )
!       sctmpb(kp,k) = sctmpb(k,kp)

      !-- Follow scatergn_weaklib.f90 in applying detailed balance here

      !-- nu on n

      phi0_nu_n(k,kp) = sct0_n * xi_n_wm(k)
      phi1_nu_n(k,kp) = sct1_n * xi_n_wm(k)

      phi0_nu_n(kp,k) = phi0_nu_n(k,kp) &
                        * fexp ( (egrid_c(kp) - egrid_c(k)) / tmev )
      phi1_nu_n(kp,k) = phi1_nu_n(k,kp) &
                        * fexp ( (egrid_c(kp) - egrid_c(k)) / tmev )

      !-- nub on n

      phi0_nub_n(k,kp) = sct0_n * xib_n_wm(k)
      phi1_nub_n(k,kp) = sct1_n * xib_n_wm(k)

      phi0_nub_n(kp,k) = phi0_nub_n(k,kp) &
                        * fexp ( (egrid_c(kp) - egrid_c(k)) / tmev )
      phi1_nub_n(kp,k) = phi1_nub_n(k,kp) &
                        * fexp ( (egrid_c(kp) - egrid_c(k)) / tmev )

      !-- nu on p

      phi0_nu_p(k,kp) = sct0_p * xi_p_wm(k)
      phi1_nu_p(k,kp) = sct1_p * xi_p_wm(k)

      phi0_nu_p(kp,k) = phi0_nu_p(k,kp) &
                        * fexp ( (egrid_c(kp) - egrid_c(k)) / tmev )
      phi1_nu_p(kp,k) = phi1_nu_p(k,kp) &
                        * fexp ( (egrid_c(kp) - egrid_c(k)) / tmev )

      !-- nub on p

      phi0_nub_p(k,kp) = sct0_p * xib_p_wm(k)
      phi1_nub_p(k,kp) = sct1_p * xib_p_wm(k)

      phi0_nub_p(kp,k) = phi0_nub_p(k,kp) &
                        * fexp ( (egrid_c(kp) - egrid_c(k)) / tmev )
      phi1_nub_p(kp,k) = phi1_nub_p(k,kp) &
                        * fexp ( (egrid_c(kp) - egrid_c(k)) / tmev )

    END DO ! kp = 1,k
  END DO ! k = nez, 1, -1

RETURN
END SUBROUTINE scatnrgn_weaklib
