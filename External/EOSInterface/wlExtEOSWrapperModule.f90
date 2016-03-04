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

USE wlKindModule, ONLY: double => dp
USE wlExtNumericalModule, ONLY: zero, one
USE wlExtPhysicalConstantsModule, ONLY: dmnp, kmev, rmu, cm3fm3, ergmev, &
& asig
USE e_p_eos_module

PUBLIC wlGetElectronEOS
PUBLIC wlGetFullEOS

PRIVATE eos_bck_nse

REAL(double), PARAMETER     :: kfm = cm3fm3/rmu    ! ( # nucleons/gram )( cm3/fm3 )
REAL(double), PARAMETER     :: kp  = ergmev/cm3fm3 ! ( erg/cm3 ) / ( mev/fm3 )
REAL(double), PARAMETER     :: ku  = ergmev/rmu    ! ( # nucleons/gram )( erg/mev )
REAL(double), PARAMETER     :: UTOT0 = 8.9d0       ! change in the zero of energy [MeV]
REAL(double), PARAMETER     :: me = 0.510998d+00   ! electron mass [MeV]

REAL(double), DIMENSION(4)  :: inpvar_save         ! retained copy of input var for LS retries

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

REAL(double), INTENT(in)    :: rho         ! Density [g/cm3]
REAL(double), INTENT(in)    :: temp        ! Temperature [K]
REAL(double), INTENT(in)    :: Ye          ! Electron fraction

!-----------------------------------------------------------------------
!  Output variables
!-----------------------------------------------------------------------

REAL(double), INTENT(out)   :: press_e     ! electron/positron/photon pressure
REAL(double), INTENT(out)   :: entrop_e    ! electron/positron/photon entropy [kb/baryon]
REAL(double), INTENT(out)   :: energ_e     ! electron/positron/photon energy
REAL(double), INTENT(out)   :: chem_e      ! electron chemical potential [MeV]

!-----------------------------------------------------------------------
!  Local variables 
!-----------------------------------------------------------------------

REAL(double)                :: pe          ! raw electron pressure [dynes cm^{-2}]
REAL(double)                :: ee          ! raw electron energy [ergs cm^{-3}]
REAL(double)                :: yeplus      ! positron fraction [unused]
REAL(double)                :: rel         ! relativity parameter [unused]
REAL(double)                :: tmev        ! temperature [MeV]
REAL(double)                :: brydns      ! LS input baryon density [fm^-3]

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

REAL(double), INTENT(in)     :: rho         ! Density [g/cm3]
REAL(double), INTENT(in)     :: temp        ! Temperature [K]
REAL(double), INTENT(in)     :: Ye          ! Electron fraction

INTEGER, INTENT(in) :: ii,jj,kk

!-----------------------------------------------------------------------
!  Output variables
!-----------------------------------------------------------------------

LOGICAL, INTENT(out)        :: fail        ! did EoS fail to converge

REAL(double), INTENT(out)   :: press       ! pressure
REAL(double), INTENT(out)   :: energ       ! internal energy
REAL(double), INTENT(out)   :: entrop      ! entropy [kb/baryon]
REAL(double), INTENT(out)   :: chem_n      ! free neutron chemical potential
REAL(double), INTENT(out)   :: chem_p      ! free proton chemical potential
REAL(double), INTENT(out)   :: chem_e      ! electron chemical potential
REAL(double), INTENT(out)   :: xn_neut     ! free neutron fraction
REAL(double), INTENT(out)   :: xn_prot     ! free proton fraction
REAL(double), INTENT(out)   :: xn_alpha    ! alpha fraction
REAL(double), INTENT(out)   :: xn_heavy    ! heavy fraction
REAL(double), INTENT(out)   :: a_heavy     ! A for mean heavy nucleus
REAL(double), INTENT(out)   :: z_heavy     ! Z for mean heavy nucleus
REAL(double), INTENT(out)   :: be_heavy    ! Binding energy for mean heavy nucleus
REAL(double), INTENT(out)   :: thermalenergy ! Internal energy without rest mass and constant offsets
REAL(double), INTENT(out)   :: gamma1      ! Gamma from LS 

!-----------------------------------------------------------------------
!  Local variables 
!-----------------------------------------------------------------------

INTEGER                     :: iflag       ! LS input type flag
INTEGER                     :: eosflg      ! LS output EoS type flag
INTEGER                     :: forflg      ! LS EoS forcing flag
INTEGER                     :: sf          ! eos success flag
INTEGER                     :: ieos        ! eos counter

REAL(double), DIMENSION(4)  :: inpvar      ! input variables for LS EoS

REAL(double)                :: pe          ! raw electron pressure [dynes cm^{-2}]
REAL(double)                :: ee          ! raw electron energy [ergs cm^{-3}]
REAL(double)                :: entrop_e    ! entropy of electron/photon gas
REAL(double)                :: yeplus      ! positron fraction [unused]
REAL(double)                :: rel         ! relativity parameter [unused]
REAL(double)                :: tmev        ! temperature [MeV]
REAL(double)                :: brydns      ! LS input baryon density [fm^-3]

REAL(double)                :: xprev       ! LS input protron fraction
REAL(double)                :: pprev       ! previous value of proton density
REAL(double)                :: t_old       ! old value of temperature

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

  CALL inveos( inpvar, t_old, ye, brydns, iflag, eosflg, forflg, sf,    &
&   xprev, pprev )

!-----------------------------------------------------------------------
!
!   ||||| If sf /= 0 (LS EOS has converged), proceed normally. |||||
!
!-----------------------------------------------------------------------

  IF ( sf /= 0 ) THEN

    inpvar_save(1) = pprev
    inpvar_save(2) = inpvar(2)
    inpvar_save(3) = inpvar(3)
    inpvar_save(4) = inpvar(4)

!-----------------------------------------------------------------------
!  Calculate XPROT and XNUT if ye < 0.03
!-----------------------------------------------------------------------

    IF ( ye < 0.03d0 ) THEN
      XPROT         = ye
      XNUT          = one - ye
    END IF

!-----------------------------------------------------------------------
!
!  ||||| If sf = 0 (LS EOS has failed to converge), call the  |||||
!  |||||   Lattimer-Swesty EOS again with inpvar_save from    |||||
!  |||||         previous call and hope it works!             |||||
!
!-----------------------------------------------------------------------

  ELSE  ! sf = 0
    pprev     = DMIN1( ye * brydns, inpvar_save(1) )
    inpvar(1) = tmev
    inpvar(2) = inpvar_save(2)
    inpvar(3) = inpvar_save(3)
    inpvar(4) = inpvar_save(4)

    CALL inveos( inpvar, t_old, ye, brydns, iflag, eosflg, forflg, sf,  &
&     xprev, pprev )

!-----------------------------------------------------------------------
!  Calculate XPROT and XNUT if ye < 0.03
!-----------------------------------------------------------------------

    IF ( ye < 0.03d0 ) THEN
      XPROT         = ye
      XNUT          = one - ye
    END IF

  END IF ! sf /= 0

!-----------------------------------------------------------------------
!  Failure mode if retry failed
!-----------------------------------------------------------------------

  IF ( sf == 0 ) THEN
    fail = .true.
    RETURN
  END IF

!-----------------------------------------------------------------------
!  Convert to units for Chimera to be returned
!-----------------------------------------------------------------------

  press    = ( PTOT - EPRESS + pe ) * kp
  energ    = ( UTOT + UTOT0 - EU + ee ) * ku
  entrop   = ( STOT - ES ) + entrop_e
  chem_n   = ( MUN - dmnp )
  chem_p   =   MUPROT
  xn_neut  =   XNUT
  xn_prot  =   XPROT
  xn_heavy =   XH
  xn_alpha =   MAX( zero, one - xn_neut - xn_prot -xn_heavy ) 
  a_heavy  =   A
  z_heavy  =   X * A
  be_heavy =   BUNUC
  thermalenergy = ku * ( UTOT - EU + ee + dmnp * XPROT + 7.075 * xn_alpha - BUNUC + 1.5d0 * tmev * XH/A - Ye * me  )
  gamma1   =   GAM_S

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

REAL(double)           :: at4           ! photon energy (MeV fm^{-3})
REAL(double)           :: srad          ! entropy (nucleon^{-1})
REAL(double)           :: dtest         ! estimate of transition density to nuclear matter

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
