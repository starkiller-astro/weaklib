MODULE wlExtEOSWrapperModule
!-----------------------------------------------------------------------
!
!    File:         wlExtEOSWrapperModule.f90
!    Module:       wlExtEOSWrapperModule
!    Type:         Module w/ Subroutines
!    Author:       E. J. Lentz, Dept. Phys. & Astronomy
!                  U. Tennesee, Knoxville
!
!    Created:      7/30/13
!    WeakLib ver:  1/13/15
!
!    Purpose:
!      Wrappers to make direct calls equation of state subroutines
!       given the (rho,T,Ye) and type of EoS required
!
!    Units:
!      Pressure      :
!      Energy(int)   :
!      Entropy       : kb/baryon
!      Chem.Pot.elec : MeV 
!
!  CONTAINS:
!    get_elec_eos
!    get_full_eos
!
!    Modules used:
!  wlKindModule
!  wlExtNumericalModule
!  wlExtPhysicalConstantsModule
!-----------------------------------------------------------------------
!  NOTE: These routines are NOT thread safe! Do not call from within
!        shared memory parallel regions (i.e., OpenMP PARALLEL regions)
!        unless specific precautions are taken to prevent threads from
!        overwriting shared "backchannel" variables used to access and
!        within EoS routines
!-----------------------------------------------------------------------

USE wlKindModule, ONLY: dp
USE wlExtNumericalModule, ONLY: zero, half, one, pi
USE wlExtPhysicalConstantsModule, ONLY: dmnp, kmev, rmu, cm3fm3, ergmev, &
& asig, hbarc, mb
USE e_p_eos_module

PUBLIC wlGetElectronEOS
PUBLIC wlGetFullEOS

PRIVATE eos_bck_nse

REAL(dp), PARAMETER     :: kfm = cm3fm3/rmu    ! ( # nucleons/gram )( cm3/fm3 )
REAL(dp), PARAMETER     :: kp  = ergmev/cm3fm3 ! ( erg/cm3 ) / ( mev/fm3 )
REAL(dp), PARAMETER     :: ku  = ergmev/rmu    ! ( # nucleons/gram )( erg/mev )
REAL(dp), PARAMETER     :: UTOT0 = 8.9d0       ! change in the zero of energy [MeV]
REAL(dp), PARAMETER     :: me = 0.510998d+00   ! electron mass [MeV]

REAL(dp), DIMENSION(4)  :: inpvar_save         ! retained copy of input var for LS retries

INTEGER, PARAMETER :: n_ist = 10

REAL(dp), DIMENSION(n_ist), PARAMETER  :: drho0 = [ (0.04d0*i,  i=1,n_ist) ]!, 0.04d0, 0.06d0 0.08d0, 0.12d0, 0.16d0, 0.20d0 /)
REAL(dp), DIMENSION(n_ist), PARAMETER  :: dt0   = [ (0.04d0*i,  i=1,n_ist) ] !0.04d0, 0.08d0, 0.12d0, 0.16d0, 0.20d0 /)
REAL(dp), DIMENSION(n_ist), PARAMETER  :: dye0  = [ (0.01d0*i, i=1,n_ist) ] !, 0.02d0, 0.03d0, 0.04d0, 0.05d0 /)

PRIVATE

CONTAINS

SUBROUTINE wlGetElectronEOS( rho, temp, ye, press_e, entrop_e, energ_e, chem_e )
!-----------------------------------------------------------------------
!
!    Module:       wlGetElectronEOS
!    Type:         Module Subroutine
!    Author:       E. J. Lentz, Dept. Phys. & Astronomy
!                  U. Tennesee, Knoxville
!
!    Date:         7/30/13
!
!    Purpose:
!      Wrapper for electron EoS only
!
!    Input arguments:
!
!  rho            : Density [g/cm3]
!  temp           : Temperature [K]
!  Ye             : Electron fraction
!
!    Output arguments:
!  press_e        : electron/positron/photon pressure
!  entrop_e       : electron/positron/photon entropy
!  energ_e        : electron/positron/photon entropy
!  chem_e         : electron chemical potential
!                   (positron chemical potential = -chem_e)
!-----------------------------------------------------------------------

IMPLICIT none

!-----------------------------------------------------------------------
!  Input variables
!-----------------------------------------------------------------------

REAL(dp), INTENT(in)    :: rho         ! Density [g/cm3]
REAL(dp), INTENT(in)    :: temp        ! Temperature [K]
REAL(dp), INTENT(in)    :: Ye          ! Electron fraction

!-----------------------------------------------------------------------
!  Output variables
!-----------------------------------------------------------------------

REAL(dp), INTENT(out)   :: press_e     ! electron/positron/photon pressure
REAL(dp), INTENT(out)   :: entrop_e    ! electron/positron/photon entropy [kb/baryon]
REAL(dp), INTENT(out)   :: energ_e     ! electron/positron/photon energy
REAL(dp), INTENT(out)   :: chem_e      ! electron chemical potential [MeV]

!-----------------------------------------------------------------------
!  Local variables 
!-----------------------------------------------------------------------

REAL(dp)                :: pe          ! raw electron pressure [dynes cm^{-2}]
REAL(dp)                :: ee          ! raw electron energy [ergs cm^{-3}]
REAL(dp)                :: yeplus      ! positron fraction [unused]
REAL(dp)                :: rel         ! relativity parameter [unused]
REAL(dp)                :: tmev        ! temperature [MeV]
REAL(dp)                :: brydns      ! LS input baryon density [fm^-3]

!-----------------------------------------------------------------------
!  Convert to units for electron EoS
!-----------------------------------------------------------------------

brydns    = rho * kfm
tmev      = temp * kmev

!-----------------------------------------------------------------------
!  Compute electron equation of state
!-----------------------------------------------------------------------

CALL initialize_e_p_eos( )
CALL e_p_eos( brydns, tmev, ye, pe, ee, entrop_e, chem_e, yeplus, rel )

!-----------------------------------------------------------------------
!  Convert to units for Chimera to be returned
!-----------------------------------------------------------------------

press_e   = pe * kp
energ_e   = ee * ku

END SUBROUTINE wlGetElectronEOS

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

SUBROUTINE wlGetFullEOS( rho, temp, ye, flag, fail, press, energ,       &
& entrop, chem_n, chem_p, chem_e, xn_neut, xn_prot, xn_alpha, xn_heavy, &
& a_heavy, z_heavy, be_heavy, thermalenergy, gamma1, ii, jj, kk )
!-----------------------------------------------------------------------
!
!    Module:       wlGetFullEOS
!    Type:         Module Subroutine
!    Author:       E. J. Lentz, Dept. Phys. & Astronomy
!                  U. Tennesee, Knoxville
!
!    Date:         8/09/13
!
!    Purpose:
!      Wrapper for full nucleonic + electron EoS
!
!    Input arguments:
!
!  rho            : Density [g/cm3]
!  temp           : Temperature [K]
!  Ye             : Electron fraction
!  Flag           : eos flag
!
!    Output arguments:
!  fail           : failure flag, true if EoS did not converge
!  press          : electron/positron/photon pressure
!  energ          : electron/positron/photon energy
!  entrop         : electron/positron/photon entropy
!  chem_n         : free neutron cemical potential
!  chem_p         : free proton chemical potential
!  chem_e         : electron chemical potential
!                   (positron chemical potential = -chem_e)
!  xn_neut        : free neutron fraction
!  xn_prot        : free proton fraction
!  xn_alpha       : alpha fraction
!  xn_heavy       : heavy nucleus fraction
!  a_heavy        : A for mean heavy nucleus
!  z_heavy        : Z for mean heavy nucleus
!  be_heavy       : binding energy for mean heavy nucleus
!  thermalenergy  : internal energy minus rest mass and constant offsets
!  gamma1         : gamma from LS 
!-----------------------------------------------------------------------

USE el_eos_module, ONLY: EPRESS, EU, ES
USE eos_m4c_module, ONLY: PTOT, UTOT, STOT, XNUT, XPROT, XH, MUN,      &
& MUPROT, A, X, GAM_S, BUNUC
USE eos_bck_module, ONLY: dbck, tbck, yebck, dtran, xnbck,  &
& xpbck, xhbck, ptotbck, uhat, etot, stotbck, un, zabck, ahbck, theta, &
& xabck, b_hvy

IMPLICIT none

!-----------------------------------------------------------------------
!  Input variables
!-----------------------------------------------------------------------

CHARACTER(len=1), INTENT(inout) :: flag     ! nuclear eos selection flag

REAL(dp), INTENT(in)    :: rho         ! Density [g/cm3]
REAL(dp), INTENT(in)    :: temp        ! Temperature [K]
REAL(dp), INTENT(in)    :: Ye          ! Electron fraction

INTEGER, INTENT(in) :: ii,jj,kk

!-----------------------------------------------------------------------
!  Output variables
!-----------------------------------------------------------------------

LOGICAL, INTENT(out)    :: fail             ! did EoS fail to converge

REAL(dp), INTENT(out)   :: press            ! pressure
REAL(dp), INTENT(out)   :: energ            ! internal energy
REAL(dp), INTENT(out)   :: entrop           ! entropy [kb/baryon]
REAL(dp), INTENT(out)   :: chem_n           ! free neutron chemical potential
REAL(dp), INTENT(out)   :: chem_p           ! free proton chemical potential
REAL(dp), INTENT(out)   :: chem_e           ! electron chemical potential
REAL(dp), INTENT(out)   :: xn_neut          ! free neutron fraction
REAL(dp), INTENT(out)   :: xn_prot          ! free proton fraction
REAL(dp), INTENT(out)   :: xn_alpha         ! alpha fraction
REAL(dp), INTENT(out)   :: xn_heavy         ! heavy fraction
REAL(dp), INTENT(out)   :: a_heavy          ! A for mean heavy nucleus
REAL(dp), INTENT(out)   :: z_heavy          ! Z for mean heavy nucleus
REAL(dp), INTENT(out)   :: be_heavy         ! Binding energy for mean heavy nucleus
REAL(dp), INTENT(out)   :: thermalenergy    ! Internal energy without rest mass and constant offsets
REAL(dp), INTENT(out)   :: gamma1           ! Gamma from LS 

!-----------------------------------------------------------------------
!  Local variables 
!-----------------------------------------------------------------------

INTEGER                 :: iflag            ! LS input type flag
INTEGER                 :: eosflg           ! LS output EoS type flag
INTEGER                 :: forflg           ! LS EoS forcing flag
INTEGER                 :: sf               ! eos success flag
INTEGER                 :: ieos             ! eos counter
INTEGER                 :: ist              ! do index for staddling rho in a failed LS EoS attempt

REAL(dp), DIMENSION(4)  :: inpvar           ! input variables for LS EoS

REAL(dp)                :: pe               ! raw electron pressure [dynes cm^{-2}]
REAL(dp)                :: ee               ! raw electron energy [ergs cm^{-3}]
REAL(dp)                :: entrop_e         ! entropy of electron/photon gas
REAL(dp)                :: yeplus           ! positron fraction [unused]
REAL(dp)                :: rel              ! relativity parameter [unused]
REAL(dp)                :: tmev             ! temperature [MeV]
REAL(dp)                :: brydns           ! LS input baryon density [fm^-3]

REAL(dp)                :: xprev            ! LS input protron fraction
REAL(dp)                :: pprev            ! previous value of proton density
REAL(dp)                :: t_old            ! old value of temperature

REAL(dp)                :: rho_p            ! rho + drho new trial density in the case of a failed LS convergence
REAL(dp)                :: rho_m            ! rho - drho new trial density in the case of a failed LS convergence
REAL(dp)                :: t_p              ! t + dt new trial temperature in the case of a failed LS convergence
REAL(dp)                :: t_m              ! t - dt new trial temperature in the case of a failed LS convergence
REAL(dp)                :: ye_p             ! ye + dye new trial electron fraction in the case of a failed LS convergence
REAL(dp)                :: ye_m             ! ye - dye new trial electron fraction in the case of a failed LS convergence
REAL(dp)                :: drho             ! fraction of rho added and subtracted
REAL(dp)                :: dt               ! fraction of temperature added and subtracted
REAL(dp)                :: dye              ! electron fraction added and subtracted

REAL(dp)                :: press_p          ! pressure for rho + drho
REAL(dp)                :: energ_p          ! internal energy for rho + drho
REAL(dp)                :: entrop_p         ! entropy [kb/baryon] for rho + drho
REAL(dp)                :: chem_n_p         ! free neutron chemical potential for rho + drho
REAL(dp)                :: chem_p_p         ! free proton chemical potential for rho + drho
REAL(dp)                :: xn_neut_p        ! free neutron fraction for rho + drho
REAL(dp)                :: xn_prot_p        ! free proton fraction for rho + drho
REAL(dp)                :: xn_heavy_p       ! free proton fraction for rho + drho
REAL(dp)                :: xn_alpha_p
REAL(dp)                :: a_heavy_p        ! A for mean heavy nucleus for rho + drho
REAL(dp)                :: z_heavy_p        ! Z for mean heavy nucleus for rho + drho
REAL(dp)                :: be_heavy_p       ! Binding energy for mean heavy nucleus for rho + drho
REAL(dp)                :: thermalenergy_p
REAL(dp)                :: gamma1_p         ! thermodynamic index gamma for rho + drho

REAL(dp)                :: press_m          ! pressure for rho - drho
REAL(dp)                :: energ_m          ! internal energy for rho - drho
REAL(dp)                :: entrop_m         ! entropy [kb/baryon] for rho - drho
REAL(dp)                :: chem_n_m         ! free neutron chemical potential for rho - drho
REAL(dp)                :: chem_p_m         ! free proton chemical potential for rho - drho
REAL(dp)                :: xn_neut_m        ! free neutron fraction for rho - drho
REAL(dp)                :: xn_prot_m        ! free proton fraction for rho - drho
REAL(dp)                :: xn_heavy_m       ! free proton fraction for rho - drho
REAL(dp)                :: xn_alpha_m
REAL(dp)                :: a_heavy_m        ! A for mean heavy nucleus for rho - drho
REAL(dp)                :: z_heavy_m        ! Z for mean heavy nucleus for rho - drho
REAL(dp)                :: be_heavy_m       ! Binding energy for mean heavy nucleus for rho - drho
REAL(dp)                :: thermalenergy_m
REAL(dp)                :: gamma1_m         ! thermodynamic index gamma for rho - drho

REAL(dp), PARAMETER     :: c0 = ( mb/( 2.d0 * pi * hbarc**2 ) )**1.5d0
REAL(dp)                :: therm            ! thermal parameter in the chemical potential
REAL(dp), PARAMETER     :: x_min = 1.d-30   ! minimum mass fraction fraction

!-----------------------------------------------------------------------
!  Convert to units for EoS and other intializations
!-----------------------------------------------------------------------

fail      = .false.
press     = zero
energ     = zero
entrop    = zero
chem_n    = zero
chem_p    = zero
chem_e    = zero
xn_neut   = zero
xn_prot   = zero
xn_alpha  = zero
xn_heavy  = zero
a_heavy   = zero
z_heavy   = zero
be_heavy  = zero
thermalenergy  = zero
gamma1    = zero
brydns    = rho * kfm
tmev      = temp * kmev
t_old     = tmev

iflag     = 1
forflg    = 0
sf        = 1

!-----------------------------------------------------------------------
!  Compute electron equation of state
!-----------------------------------------------------------------------

CALL initialize_e_p_eos( )
CALL e_p_eos( brydns, tmev, ye, pe, ee, entrop_e, chem_e, yeplus, rel )

SELECT CASE ( flag )

  CASE ('L') 

!-----------------------------------------------------------------------
! Lattimer-Swesty EoS
!-----------------------------------------------------------------------
pprev     = ye * brydns
inpvar(1) = tmev
inpvar(2) = 0.155d0
inpvar(3) = -15.d0
inpvar(4) = -10.d0

!write(*,'(es12.3,f5.2,f5.2)') rho, tmev, ye 
!write(*,*) brydns, tmev, ye


CALL inveos( inpvar, t_old, ye, brydns, iflag, eosflg, forflg, sf,      &
&   xprev, pprev )
!-----------------------------------------------------------------------
!
!   ||||| If sf /= 0 (LS EOS has converged), proceed normally. |||||
!
!-----------------------------------------------------------------------

IF ( sf == 1 ) THEN

  inpvar_save(1) = pprev
  inpvar_save(2) = inpvar(2)
  inpvar_save(3) = inpvar(3)
  inpvar_save(4) = inpvar(4)

!-----------------------------------------------------------------------
!  Convert to units for Chimera to be returned
!-----------------------------------------------------------------------
  press         = ( PTOT - EPRESS + pe ) * kp
  energ         = ( UTOT + UTOT0 - EU + ee ) * ku
  entrop        = ( STOT - ES ) + entrop_e
  chem_n        = ( MUN - dmnp )
  chem_p        =   MUPROT
  xn_neut       =   XNUT
  xn_prot       =   XPROT
  xn_heavy      =   XH
  xn_alpha      =   MAX( zero, one - xn_neut - xn_prot -xn_heavy ) 
  a_heavy       =   A
  z_heavy       =   X * A
  be_heavy      =   BUNUC
  thermalenergy =   ku * ( UTOT - EU + ee + dmnp * XPROT + 7.075 * xn_alpha &
                -   BUNUC + 1.5d0 * tmev * XH/A - Ye * me  )
  gamma1        =   GAM_S

!-----------------------------------------------------------------------
!  Calculate XPROT and XNUT if ye < 0.03
!-----------------------------------------------------------------------

  IF ( ye < 0.03d0 ) THEN
    XPROT         = ye
    XNUT          = one - ye
  END IF

!-----------------------------------------------------------------------
!
!  ||||| If sf /= 0 (LS EOS has failed to converge), call the |||||
!  |||||   Lattimer-Swesty EOS again with straddling the      |||||
!  |||||         density and average the results.             |||||
!  |||||   Successively widen the straddle until convergence  |||||
!  |||||           until convergence is attained.             |||||
!
!-----------------------------------------------------------------------

ELSE  ! sf /= 1
  DO ist = 1, n_ist
!write(*,*) 'LS call failed, trying loop'
!write(*,'(i1,es12.3,f5.2,f5.2)') ist, rho, tmev, ye 
  

!-----------------------------------------------------------------------
!  Try converging with rhod --> rhod + drho(ist) and td + dt(ist)
!   or ye --> ye + dye
!  Determine the variable increment
!-----------------------------------------------------------------------

    IF ( temp < 1.1d+10  .and.  rho < 2.d+12  .and.  ye > 0.154d0  .and.  ye < 0.172d0 ) THEN
      dye           = dye0(ist)
      drho          = zero
      dt            = zero
    ELSE
      dye           = zero
      drho          = drho0(ist)
      dt            = dt0(ist)
    END IF ! d < 1.1d+10  .and.  rhod < 2.d+12  .and.  ye > 0.154d0  .and.  ye < 0.172d0

    rho_p           = ( 1.d0 + drho ) * rho
    t_p             = ( 1.d0 + dt ) * temp
    ye_p            = ye + dye
    brydns          = rho_p * kfm
    tmev            = t_p * kmev

    pprev     = ye_p * brydns !ye * brydns
    inpvar(1) = tmev
    inpvar(2) = inpvar_save(2) !0.155d0
    inpvar(3) = inpvar_save(3) !-15.d0
    inpvar(4) = inpvar_save(4) !-10.d0
  
!    pprev = ye * brydns
!    inpvar(1) = tmev
!    inpvar(2) = 0.155d0
!    inpvar(3) = -15.d0
!    inpvar(4) = -10.d0

    CALL inveos( inpvar, t_old, ye_p, brydns, iflag, eosflg, forflg, sf,  &
&     xprev, pprev )

!-----------------------------------------------------------------------
!  If LS EoS fails to converge, try straddling with a greater density
!   width
!-----------------------------------------------------------------------

    IF ( sf /= 1 )  THEN
      CYCLE
    END IF ! sf /= 1

    inpvar_save(1) = pprev
    inpvar_save(2) = inpvar(2)
    inpvar_save(3) = inpvar(3)
    inpvar_save(4) = inpvar(4)

!-----------------------------------------------------------------------
!  Calculate XPROT and XNUT if ye < 0.03
!-----------------------------------------------------------------------

    IF ( ye_p < 0.03d0 ) THEN
      XPROT   = ye_p
      XNUT    = one - ye_p
    END IF

!-----------------------------------------------------------------------
!  Convert to units for Chimera to be returned
!-----------------------------------------------------------------------
  press_p         = ( PTOT - EPRESS + pe ) * kp
  energ_p         = ( UTOT + UTOT0 - EU + ee ) * ku
  entrop_p        = ( STOT - ES ) + entrop_e
  chem_n_p        = ( MUN - dmnp )
  chem_p_p        =   MUPROT
  xn_neut_p       =   XNUT
  xn_prot_p       =   XPROT
  xn_heavy_p      =   XH
  xn_alpha_p      =   MAX( zero, one - xn_neut_p - xn_prot_p -xn_heavy_p ) 
  a_heavy_p       =   A
  z_heavy_p       =   X * A
  be_heavy_p      =   BUNUC
  thermalenergy_p =   ku * ( UTOT - EU + ee + dmnp * XPROT + 7.075 * xn_alpha_p &
                  -   BUNUC + 1.5d0 * tmev * XH/A - Ye_p * me  )
  gamma1_p        =   GAM_S

!-----------------------------------------------------------------------
!  Now try converging with rho --> rho - drho(ist)
!-----------------------------------------------------------------------

    rho_m           = ( 1.d0 - drho ) * rho
    t_m             = ( 1.d0 - dt ) * temp
    ye_m            = ye - dye
    brydns          = rho_m * kfm
    tmev            = t_m * kmev

    pprev     = ye_m * brydns !ye * brydns
    inpvar(1) = tmev
    inpvar(2) = inpvar_save(2) !0.155d0
    inpvar(3) = inpvar_save(3) !-15.d0
    inpvar(4) = inpvar_save(4) !-10.d0

!    pprev     = ye * brydns
!    inpvar(1) = tmev
!    inpvar(2) = 0.155d0
!    inpvar(3) = -15.d0
!    inpvar(4) = -10.d0

    CALL inveos( inpvar, t_old, ye_m, brydns, iflag, eosflg, forflg, sf,  &
&     xprev, pprev )

!-----------------------------------------------------------------------
!  If LS EoS fails to converge, try straddling with a greater density
!   width
!-----------------------------------------------------------------------

    IF ( sf /= 1 )  THEN
      CYCLE
    END IF ! sf /= 1

    inpvar_save(1) = pprev
    inpvar_save(2) = inpvar(2)
    inpvar_save(3) = inpvar(3)
    inpvar_save(4) = inpvar(4)

!-----------------------------------------------------------------------
!  Calculate XPROT and XNUT if ye < 0.03
!-----------------------------------------------------------------------

    IF ( ye_m < 0.03d0 ) THEN
      XPROT   = ye_m
      XNUT    = one - ye_m
    END IF

!-----------------------------------------------------------------------
!  Convert to units for Chimera to be returned
!-----------------------------------------------------------------------
  press_m         = ( PTOT - EPRESS + pe ) * kp
  energ_m         = ( UTOT + UTOT0 - EU + ee ) * ku
  entrop_m        = ( STOT - ES ) + entrop_e
  chem_n_m        = ( MUN - dmnp )
  chem_p_m        =   MUPROT
  xn_neut_m       =   XNUT
  xn_prot_m       =   XPROT
  xn_heavy_m      =   XH
  xn_alpha_m      =   MAX( zero, one - xn_neut_m - xn_prot_m -xn_heavy_m ) 
  a_heavy_m       =   A
  z_heavy_m       =   X * A
  be_heavy_m      =   BUNUC
  thermalenergy_m =   ku * ( UTOT - EU + ee + dmnp * XPROT + 7.075 * xn_alpha_m &
                  -   BUNUC + 1.5d0 * tmev * XH/A - Ye_m * me  )
  gamma1_m        =   GAM_S

    press         = 0.5d0 * ( press_p    + press_m    )
    energ         = 0.5d0 * ( energ_p    + energ_m    )
    entrop        = 0.5d0 * ( entrop_p   + entrop_m   )
    chem_n        = 0.5d0 * ( chem_n_p   + chem_n_m   )
    chem_p        = 0.5d0 * ( chem_p_p   + chem_p_m   )
    xn_neut       = 0.5d0 * ( xn_neut_p  + xn_neut_m  )
    xn_prot       = 0.5d0 * ( xn_prot_p  + xn_prot_m  )
    xn_heavy      = 0.5d0 * ( xn_heavy_p + xn_heavy_m )
    xn_alpha      = 0.5d0 * ( xn_alpha_p + xn_alpha_m )
    a_heavy       = 0.5d0 * ( a_heavy_p  + a_heavy_m  )
    z_heavy       = 0.5d0 * ( z_heavy_p  + z_heavy_m  )
    be_heavy      = 0.5d0 * ( be_heavy_p + be_heavy_m )
    thermalenergy = 0.5d0 * (thermalenergy_p + thermalenergy_m)
    gamma1        = 0.5d0 * ( gamma1_p   + gamma1_m    )


    EXIT

  END DO ! ist = 1, 5

  IF ( sf == 0 ) THEN
    fail = .true.
write(*,*) 'Straddling calls to LS EOS failed at: ', tmev, brydns, ye
write(*,*) rho, tmev, ye
    RETURN
  END IF

END IF ! sf == 1

!-----------------------------------------------------------------------
!  Adjust xn_neut and xn_prot if ye > 0.51 to cover regions out of
!   range of the LS EOS.
!-----------------------------------------------------------------------

IF ( ye > 0.50d0 ) THEN
  xn_alpha       = DMAX1( one - xn_neut - xn_prot - xn_heavy, zero )
  XNUT           = xn_neut
  XPROT          = xn_prot
  IF ( xn_neut + xn_prot + 0.5d0 * xn_alpha + xn_heavy * ( z_heavy/a_heavy ) > ye ) THEN
    xn_prot      = DMAX1( ye - 0.5d0 * xn_alpha - xn_heavy * z_heavy/a_heavy, x_min )
    xn_neut      = DMAX1( XNUT + XPROT - xn_prot, x_min )
    therm        = brydns/( c0 * tmev * SQRT(tmev) )
    chem_n       = tmev * DLOG( half * therm * xn_neut ) - dmnp
    chem_p       = tmev * DLOG( half * therm * xn_prot ) - dmnp
  END IF ! xn_neut + xn_prot + 0.5d0 * xn_alpha + xn_heavy * ( z_heavy/a_heavy ) > ye
END IF ! ye > 0.51d0

RETURN

!-----------------------------------------------------------------------
!  BCK EoS (NSE only!!!)
!-----------------------------------------------------------------------

  CASE ('B')


! --- eosdtgen_x ---

!-----------------------------------------------------------------------
!  Call the nuclear equation of state
!-----------------------------------------------------------------------

  dbck     = brydns
  tbck     = tmev
  yebck    = ye
  un       = zero
  uhat     = zero
  theta    = zero
  zabck    = zero
  xabck    = zero
  dtran    = zero

  CALL eos_bck_nse

  IF ( xnbck <= zero  .or.  xpbck <= zero  .or.  stotbck <= zero ) THEN
    dtran       = dbck
    CALL eos_bck_nse
  END IF ! xnbck < 0 or xpbck < 0 or stot < 0

!-----------------------------------------------------------------------
!  Convert thermodynamic quantities from units of mev and fm to cgs
!   units.
!-----------------------------------------------------------------------

  press    = ( ptotbck + pe ) * kp
  energ    = ( etot - yebck * dmnp + ee ) * ku
  entrop   = stotbck + entrop_e
  chem_n   = un
  chem_p   = un - uhat
  IF ( dbck >= dtran ) then
    xnbck  = one - ye
    xpbck  = ye
    xhbck  = zero
  END IF ! dbck >= dtran
  xn_neut  = DMAX1( xnbck , zero )
  xn_prot  = DMAX1( xpbck , zero )
  xn_heavy = DMAX1( xhbck , zero )
  xn_alpha = MAX( zero, one - xn_neut - xn_prot -xn_heavy ) 
  a_heavy  = ahbck
  z_heavy  = zabck * ahbck
  be_heavy = b_hvy
  thermalenergy = ku * ( etot + ee - UTOT0 + 7.075 * xn_alpha - xhbck * b_hvy - yebck * me )

!-----------------------------------------------------------------------
!  Stop if unimplemented EoS is requested
!-----------------------------------------------------------------------

  CASE DEFAULT

  WRITE(*,*) "Bad EoS flag in get_full_eos ", flag
  STOP "Bad EoS flag in get_full_eos"

END SELECT

END SUBROUTINE wlGetFullEOS

!-----------------------------------------------------------------------
!  PRIVATE copy of bck "eos" driver routine for NSE computations only
!-----------------------------------------------------------------------

SUBROUTINE eos_bck_nse

USE eos_bck_module, ONLY: dbck, tbck, yebck, erad, prad, enu, pnu,      &
& sneu, dtran, etot, eh, ed, ptotbck, ph, pd, stotbck, sh, sd, dedt, un

IMPLICIT none

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

REAL(dp)           :: at4           ! photon energy (MeV fm^{-3})
REAL(dp)           :: srad          ! entropy (nucleon^{-1})
REAL(dp)           :: dtest         ! estimate of transition density to nuclear matter

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!*************************************************** radiation *********

at4             = asig * tbck * tbck * tbck * tbck
erad            = at4/dbck
prad            = at4/3.d0
srad            = ( erad + prad/dbck )/tbck
enu             = erad
pnu             = prad
sneu            = srad

!**************************************************** nucleons *********

!........Matter is in NSE

dtest           = .16d0 * ( 1.d0 - 3.d0 * ( .5d0 - yebck )**2 )
IF ( dtran <= 0.0d0 ) dtran = dtest
dtran           = DMIN1( dtran, dtest )

!........Transition to nuclear matter

IF ( dbck >= dtran ) THEN
  CALL eosnm
ELSE
  CALL saha
END IF

!************************************************ get totals ***********
!-----------------------------------------------------------------------
!        Energy offset
!
!  The subroutine  eosnm, saha, and eosnuc_x return the energy including
!   the binding energy.
!  The quantity - yebck * dmnp is added so that the binding energy is
!   relative to free neutrons rather than free neutrons and protons
!  An energy offset 8.9d0 MeV/nucleon is added to prevent the internal
!   energy from becoming negative
!  Subtracting dmnp from un also subtracts it from up since the latter
!   is given by un - uhat. This normalizes these chemical potentials to
!   those given by the LS EOS.
!-----------------------------------------------------------------------

etot            = eh + erad + ed + UTOT0
ptotbck         = ph + prad + pd
stotbck         = sh + srad + sd
un              = un - dmnp

!***************************************** approx for energy inversion *

dedt            = sh + erad*4.d0/tbck + ed/tbck

RETURN
END SUBROUTINE eos_bck_nse

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

END MODULE wlExtEOSWrapperModule
