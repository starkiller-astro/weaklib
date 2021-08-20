!-----------------------------------------------------------------------
!    Module:       nes_module
!    Author:       S. W. Bruenn
!    Date:         6/26/03
!-----------------------------------------------------------------------

MODULE nes_module

USE wlKindModule

!-----------------------------------------------------------------------
!  Neutrino energies
!-----------------------------------------------------------------------

REAL(dp)                :: w             ! e_in/kt
REAL(dp)                :: w2            ! w**2
REAL(dp)                :: w3            ! w**3
REAL(dp)                :: w4            ! w**4
REAL(dp)                :: w5            ! w**5
REAL(dp)                :: w6            ! w**6

REAL(dp)                :: wp            ! e_out/kt
REAL(dp)                :: wp2           ! wp**2
REAL(dp)                :: wp3           ! wp**3
REAL(dp)                :: wp4           ! wp**4
REAL(dp)                :: wp5           ! wp**5
REAL(dp)                :: wp6           ! wp**6

REAL(dp)                :: wwp           ! w*wp
REAL(dp)                :: w2pwp2        ! wwp**2

REAL(dp)                :: wp_w          ! wp - w
REAL(dp)                :: w_wp          ! w - wp
REAL(dp)                :: w_wp2         ! w_wp**2
REAL(dp)                :: w_wp3         ! w_wp**3
REAL(dp)                :: w_wp4         ! w_wp**4
REAL(dp)                :: w_wp5         ! w_wp**5

REAL(dp)                :: w_p_wp        ! w + wp
REAL(dp)                :: w_p_wp2       ! w_p_wp**2
REAL(dp)                :: w_p_wp3       ! w_p_wp**3

REAL(dp)                :: e             ! incident electron energy/kt
REAL(dp)                :: e2            ! e**2
REAL(dp)                :: e3            ! e**3
REAL(dp)                :: e4            ! e**4
REAL(dp)                :: e5            ! e**5
REAL(dp)                :: e6            ! e**6
REAL(dp)                :: e7            ! e**7
REAL(dp)                :: e_min         ! minimum incident electron energy/kt

!-----------------------------------------------------------------------
!  Constants
!-----------------------------------------------------------------------

REAL(dp)                :: r2_3         = 2.d+00/3.d+00
REAL(dp)                :: r2_15        = 2.d+00/15.d+00
REAL(dp)                :: r2_105       = 2.d+00/105.d+00
REAL(dp)                :: r4_3         = 4.d+00/3.d+00
REAL(dp)                :: r4_5         = 4.d+00/5.d+00
REAL(dp)                :: r4_15        = 4.d+00/15.d+00
REAL(dp)                :: r8_3         = 8.d+00/3.d+00
REAL(dp)                :: r8_5         = 8.d+00/5.d+00
REAL(dp)                :: r8_7         = 8.d+00/7.d+00
REAL(dp)                :: r8_15        = 8.d+00/15.d+00
REAL(dp)                :: r12_5        = 12.d+00/5.d+00
REAL(dp)                :: r12_35       = 12.d+00/35.d+00
REAL(dp)                :: r16_3        = 16.d+00/3.d+00
REAL(dp)                :: r16_5        = 16.d+00/5.d+00
REAL(dp)                :: r16_30       = 16.d+00/30.d+00
REAL(dp)                :: r16_35       = 16.d+00/35.d+00
REAL(dp)                :: r28_5        = 28.d+00/5.d+00
REAL(dp)                :: r28_15       = 28.d+00/15.d+00
REAL(dp)                :: r36_105      = 36.d+00/105.d+00


END module nes_module
