SUBROUTINE nu_N_absr_momts( enu_in, tmev, m_trgt_i, m_trgt_f, m_lep, &
& cmp_trgt_i, cmp_trgt_f, cmp_lep, ab_r0, ab_r1, e_out )
!-----------------------------------------------------------------------
!
!    File:         nu_N_absr_momts
!    Module:       nu_N_absr_momts
!    Type:         Subroutine
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         5/25/02
!
!    Purpose:
!      To compute moments of the neutrino emission and absorption
!       inverse mean free paths.
!
!    Input arguments:
!   enu_in     : absorbed or emitted neutrino [MeV]
!   tmev       : temperature [MeV]
!   m_trgt_i   : mass of the target particle [MeV]
!   m_trgt_f   : mass of the transformed target particle [MeV]
!   m_lep      : mass of the created lepton [MeV]
!   cmp_trgt_i : neutron chemical potential (including rest mass) [MeV]
!   cmp_trgt_f : proton chemical potential (including rest mass) [MeV]
!   cmp_lep    : electron chemical potential (including rest mass) [MeV]
!
!    Output arguments:
!   ab_r0      : zero angular moment of the inverse absorption mean free
!                 path (cm^{-1})
!   ab_r1      : first angular moment of the inverse absorption mean
!                 free path (cm^{-1})
!   e_out      : mean energy of the emitted lepton [MeV]
!
!    Subprograms called:
!      cc_difcs
!
!    Include files:
!  wlKindModule
!  numerical_module
!  physcnst_module
!
!  abem_module
!
!-----------------------------------------------------------------------

USE wlKindModule, ONLY: dp
USE numerical_module, ONLY: zero, half, one, epsilon
USE physcnst_module, ONLY: cvel, gv, ga

USE abem_module_weaklib, ONLY:  nleg_a, x_a, nleg_e, x_e, wt_e, wt_a

IMPLICIT none

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

REAL(dp), INTENT(in)    :: enu_in        ! incoming neutrino energy
REAL(dp), INTENT(in)    :: tmev          ! temperature [MeV]
REAL(dp), INTENT(in)    :: m_trgt_i      ! mass of the initial target particle [MeV]
REAL(dp), INTENT(in)    :: m_trgt_f      ! mass of the final target particle [MeV]
REAL(dp), INTENT(in)    :: m_lep         ! mass of the final lepton [MeV]
REAL(dp), INTENT(in)    :: cmp_trgt_i    ! chemical potential of the initial target particle [MeV]
REAL(dp), INTENT(in)    :: cmp_trgt_f    ! chemical potential of the transformed target particle [MeV]
REAL(dp), INTENT(in)    :: cmp_lep       ! chemcal potential of the secondary lepton [MeV]

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(dp), INTENT(out)   :: ab_r0         ! zero angular moment of the inverse absorption mean free path (cm^{-1})
REAL(dp), INTENT(out)   :: ab_r1         ! first angular moment of the inverse absorption mean free path (cm^{-1})
REAL(dp), INTENT(out)   :: e_out         ! mean energy of the emitted lepton

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

REAL(dp), PARAMETER     :: m_cm = 1.d-2  ! meters/cm
REAL(dp)                :: m_trgt_i2     ! m_trgt_i^{2}
REAL(dp)                :: m_trgt_f2     ! m_trgt_f^{2}
REAL(dp)                :: m_lep2        ! m_lep/cm^{2}
REAL(dp)                :: dm_trgt       ! m_trgt_i - m_trgt_f
REAL(dp)                :: dm_trgtp      ! - dm_trgt

INTEGER                     :: i_a           ! summation index of angular Gauss-Lagendre quadrature
REAL(dp)                :: xu_a          ! upper limit of lepton angular quadrature
REAL(dp)                :: xl_a          ! lower limit of lepton angular quadrature
REAL(dp)                :: mid_a         ! midpoint of lepton angular quadrature
REAL(dp)                :: width_a       ! half-width of lepton angular quadrature
REAL(dp)                :: c_a           ! scaled points of angular quadrature
REAL(dp)                :: costh         ! points of angular quadrature

REAL(dp)                :: ab_r0_e       ! zero moment of the neutrino absorption inverse mean free path per angle
REAL(dp)                :: ab_r1_e       ! first moment of the neutrino absorption inverse mean free path per angle
REAL(dp)                :: e_out_e       ! parameter to estimate energy quadrature widths per angle

REAL(dp), PARAMETER     :: mult = 10.d0  ! boundary multiplier
REAL(dp)                :: t_m           ! paremeter of the energy limits
REAL(dp)                :: prin          ! paremeter of the energy limits
REAL(dp)                :: radical       ! paremeter of the energy limits

INTEGER                     :: i_e           ! summation index of ouitgoing lepton energy Gauss-Lagendre quadrature
REAL(dp)                :: xu_e          ! upper limit of lepton energy quadrature
REAL(dp)                :: xl_e          ! lower limit of lepton energy quadrature
REAL(dp)                :: mid_e         ! midpoint of lepton energy quadrature
REAL(dp)                :: width_e       ! half-width of lepton energy quadrature
REAL(dp)                :: c_e           ! scaled points of energy quadrature
REAL(dp)                :: e_lep_f       ! energy quadrature points

REAL(dp)                :: absr_eomega   ! absorption cross section per unit energy and solid angle

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

m_trgt_i2          = m_trgt_i * m_trgt_i
m_trgt_f2          = m_trgt_f * m_trgt_f
m_lep2             = m_lep * m_lep
dm_trgt            = m_trgt_i - m_trgt_f
dm_trgtp           = - dm_trgt

!-----------------------------------------------------------------------
!  Initialize for dp (angle and energy) integration
!-----------------------------------------------------------------------

ab_r0              = zero
ab_r1              = zero
e_out              = zero

!-----------------------------------------------------------------------
!
!           \\\\\ INTEGRATION OVER FINAL LEPTON ANGLE /////
!
!-----------------------------------------------------------------------

xu_a               = one
xl_a               = - one
mid_a              = half * ( xu_a + xl_a )
width_a            = half * ( xu_a - xl_a)
t_m                = 2.d0 * mult * tmev/m_trgt_i
t_m                = DMIN1( t_m, 0.999d0 )


DO i_a = 1,nleg_a
  c_a              = x_a(i_a) * width_a
  costh            = mid_a + c_a

!-----------------------------------------------------------------------
!  Initialize for final lepton energy integration
!-----------------------------------------------------------------------

  ab_r0_e          = zero
  ab_r1_e          = zero
  e_out_e          = zero

!-----------------------------------------------------------------------
!
!          \\\\\ INTEGRATION OVER GFINAL LEPTON ENERGY /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Determine the energy quadrature limits
!-----------------------------------------------------------------------

  prin            = one - t_m * costh
  radical         = DSQRT( DMAX1( 2.d0 * t_m * ( one - costh ) - t_m * t_m * ( one - costh * costh ), epsilon ) )
  xu_e            = enu_in * ( prin + radical )/( one - t_m ) + dm_trgt
  xl_e            = enu_in * ( prin - radical )/( one - t_m ) + dm_trgt
  xl_e            = DMAX1( xl_e, zero   )

!      IF ( i_a == 2   .or.  i_a == 32  .or.  i_a == 64  .or.  i_a == 16  .or.  i_a == 48 ) &
!& WRITE (6,3005) t_m,prin,radical,xu_e,xl_e
! 3005 format (' t_m=',1pe10.3,' prin=',1pe10.3,' radical=',1pe10.3,' xu_e=',1pe10.3,' xl_e=',1pe10.3)

!-----------------------------------------------------------------------
!  Integrate over final lepton energy
!-----------------------------------------------------------------------

  mid_e            = half * ( xu_e + xl_e )
  width_e          = half * ( xu_e - xl_e )

  DO i_e = 1,nleg_e
    c_e            = x_e(i_e) * width_e
    e_lep_f        = mid_e + c_e
    IF ( e_lep_f <= zero ) CYCLE
    CALL cc_difcs( enu_in, e_lep_f, costh, tmev, m_trgt_i, dm_trgtp, &
&    cmp_trgt_i, cmp_trgt_f, cmp_lep, gv, ga, absr_eomega )

    absr_eomega    = DMAX1( absr_eomega, zero )
    ab_r0_e        = ab_r0_e + absr_eomega * wt_e(i_e)
    ab_r1_e        = ab_r1_e + absr_eomega * wt_e(i_e) * costh
    e_out_e        = e_out_e + absr_eomega * wt_e(i_e) * e_lep_f
  END DO ! i_e = 1,nleg_e
 
  ab_r0_e          = ab_r0_e * width_e
  ab_r1_e          = ab_r1_e * width_e
  e_out_e          = e_out_e * width_e

!-----------------------------------------------------------------------
!  Sum final result over angle and average the final electron energy
!-----------------------------------------------------------------------

  ab_r0            = ab_r0 + ab_r0_e * wt_a(i_a)
  ab_r1            = ab_r1 + ab_r1_e * wt_a(i_a)
  e_out            = e_out + e_out_e * wt_a(i_a) 

END DO ! i_a = 1,nleg_a

e_out              = e_out/( ab_r0 + epsilon )
ab_r0              = m_cm * ab_r0 * width_a
ab_r1              = m_cm * ab_r1 * width_a

RETURN
END SUBROUTINE nu_N_absr_momts
