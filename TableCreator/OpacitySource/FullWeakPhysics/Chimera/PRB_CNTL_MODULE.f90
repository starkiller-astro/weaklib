!-----------------------------------------------------------------------
!    Module:       prb_cntl_module
!    Author:       S. W. Bruenn
!    Date:         8/16/02
!-----------------------------------------------------------------------

MODULE prb_cntl_module

USE kind_module

!-----------------------------------------------------------------------
!
!             \\\\\ ZERO TRANSVERSE VELOCITY OPTION /////
!
!-----------------------------------------------------------------------
!  v_trans_0      : zero transverse velocity above shock switch
!
!     v_trans_0 = ye : transverse velocities 0 above shock
!     v_trans_0 = no : transverse velocities 0 above shock computed
!-----------------------------------------------------------------------

CHARACTER(LEN=2)                           :: v_trans_0

!-----------------------------------------------------------------------
!  Hydrodynamic Controls
!-----------------------------------------------------------------------
!  ihydro   : Hydrodynamics switch.
!
!     ihydro = 0   : hydrodynamics bypassed.
!     ihydro = 1   : hydrodynamics included.
!
!  irelhy   : Newtonian - GR hydrodynamics toggle.
!
!     irelhy = 0   : Newtonian hydrodynamics.
!     irelhy = 1   : PN hydrodynamics.
!     irelhy = 2   : general relativistic hydrodynamics.
!
!  ilapsehy : Newtonian - PN - GR lapse toggle.
!
!     ilapsehy = 0 : lapse = 1 in the hydrodynamics and transport.
!     ilapsehy = 1 : lapse computed to PN approximation in the /; b4..040
!      hydrodynamics if irelhy = 1, and in the transport if ireltrns = 1.
!     ilapsehy = 2 : lapse computed with full GR in the 
!      hydrodynamics if irelhy = 1, and in the transport if ireltrns = 1.
!
!-----------------------------------------------------------------------

INTEGER                                    :: ihydro
INTEGER                                    :: irelhy
INTEGER                                    :: ilapsehy

!-----------------------------------------------------------------------
!  Hydrodynamic Mode Control
!-----------------------------------------------------------------------
!
!     e_mode = 'e  ' : evolve total energy (internal + gravitational + kinetic)
!     e_mode = 'e-g' : evolve partial energy (internal + kinetic)
!-----------------------------------------------------------------------

CHARACTER(LEN=3)                           :: e_mode

!-----------------------------------------------------------------------
!  Gravitation Options
!-----------------------------------------------------------------------
!  i_grav     : gravitation key
!
!     i_grav = 1 : spherical symmetric gravity
!     i_grav = 2 : nonspherical components computed from by a Newtonian
!                   Poisson solver
!
!-----------------------------------------------------------------------

INTEGER                                    :: i_grav

!-----------------------------------------------------------------------
!
!                     \\\\\ MOVING GRID OPTION /////
!
!-----------------------------------------------------------------------
!  m_grid_1 = la : Radial grid above the inner-outer grid boundary
!                   is Lagrangian if lagr and eul = no. (Control)
!  m_grid_1 = no : Radial grid above the inner-outer grid boundary
!                   is Eulerian if lagr and eul = no. (Control)
!
!  m_grid_2       This option is not currently functional
!
!-----------------------------------------------------------------------

CHARACTER(LEN=2)                           :: m_grid_1
CHARACTER(LEN=2)                           :: m_grid_2

!-----------------------------------------------------------------------
!  Transport controls
!-----------------------------------------------------------------------
!  inutrn   : Transport switch.
!
!     inutrn = 0 : neutrino transport bypassed.
!     inutrn = 1 : neutrino transport included.
!
!  ireltrns : Newtonian - GR transport toggle.
!
!     ireltrns = 0 : Newtonian neutrino transport.
!     ireltrns = 1 : general relativistic neutrino transport.
!
!  idiff    : Transport boundary conditions toggle.
!
!     idiff = 0 : dc(i,k,n) computed normally at each iteration.
!     idiff = 1 : dc(i,k,n) is the average of the current and preceding
!      timestep value.
!     idiff = 2 : dc(i,k,n) = 0 for i = inumax, dc(i,k,n) computed normally
!      at each iteration for 
!      i /= inumax.
!     idiff = 3 : dc(i,k,n) = 0 for i = imax, dc(i,k,n) as in idiff = 1
!      for i /= imax
!     idiff = 4 : dc(i,k,n) = 0 for all i.
!     idiff = 5 : dc(i,k,n) is computed only on the first iteration
!     idiff = 6 : dc(i,k,n) = 0 for i = imax, computed only on the first
!      iteration for i /= imax
!
!  inumin   : inner radial zone for which neutrino transport is computed.
!  inumax   : outer radial zone for which neutrino transport is computed.
!
!  itnu(n)  : transport temperature change switch.
!
!     itnu(n) = 0 : temperature change due to n-type neutrinos bypassed.
!     itnu(n) = 1 : temperature change due to n-type neutrinos included.
!
!  jexect   : transport temperature change for zone jexect switch.
!
!     jexect = 0 : computation of the temperature is normal.
!     jexect = j : if itnu = 0, computation of temperature is normal for
!      j = jexect.
!
!  iyenu    : transport electron fraction change switch.
!
!     iyenu = 0 : electron fraction change due to n-type neutrinos bypassed.
!     iyenu = 1 : electron fraction change due to n-type neutrinos included.
!
!  jexecy   : transport electron fraction change for zone jexect switch.
!
!     jexecy = 0 : computation of the temperature is normal.
!     jexecy = j : if iyenu = 0, computation of temperature is normal for
!      j = jexecy.
!-----------------------------------------------------------------------

INTEGER                                    :: inutrn
INTEGER                                    :: ireltrns
INTEGER                                    :: idiff
INTEGER                                    :: inumin
INTEGER                                    :: inumax
INTEGER, ALLOCATABLE, DIMENSION(:)         :: itnu
INTEGER                                    :: jexect
INTEGER                                    :: iyenu
INTEGER                                    :: jexecy

!-----------------------------------------------------------------------
!
!         \\\\\ LAGRANGIAN - EULERIAN TRANSPORT SWITCHES /////
!
!-----------------------------------------------------------------------
!     lag_tr  = ye:
!     eul_tr  = ye:
!  Lagrangian transport when moving grid option on, otherwise Eulerian
!   transport. (Control)
!
!     lag_tr  = ye:
!     eul_tr  = no:
!  Lagrangian transport always. (Control)
!
!     lag_tr  = no:
!     eul_tr  = ye:
!  Eulerian transport always. (Control)
!-----------------------------------------------------------------------

CHARACTER(LEN=2)                           :: lag_tr
CHARACTER(LEN=2)                           :: eul_tr

!-----------------------------------------------------------------------
!
!        \\\\\ TEMPERATURE AND ELECTRON FRACTION UPDATES /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  t_chg   = ye: compute matter temperature change due to source and
!                transport
!  t_chg   = no: temperature change to source and transport bypassed
!
!  ye_chg  = ye: compute electron fraction change due to source and
!                transport
!  ye_chg  = no: electron fraction change to source and transport
!                bypassed
!-----------------------------------------------------------------------

CHARACTER(LEN=2)                           :: t_chg
CHARACTER(LEN=2)                           :: ye_chg

!-----------------------------------------------------------------------
!  Absorption and emission controls
!-----------------------------------------------------------------------
!  iaefnp   : emission and absorption of e-neutrinos and e-antineutrinos
!   on free neutrons and protons switch.
!
!     iaefnp = 0: emission and absorption of e-neutrinos and
!      e-antineutrinos on free neutrons and protons turned off.
!     iaefnp = 1: emission and absorption of e-neutrinos and
!      e-antineutrinos on free neutrons and protons turned on.
!
!  i_aeps   : emission and aborption of e-neutrinos and e-antineutrinos
!   on free neutrons and protons phase space switch.
!
!     i_aeps = 0: emission and absorption of e-neutrinos and
!      e-antineutrinos on free neutrons and protons are computed
!      assuming no recoil, or thermal motions, and approximate
!      nucleon phase space blocking.
!     iaefnp = 1: emission and absorption of e-neutrinos and
!      e-antineutrinos on free neutrons and protons are computed
!      taking into account recoil, thermal motions, and nucleon
!      phase space blocking.
!
!  rhoaefnp : density above which emission and absorption of e-neutrinos
!   and e-antineutrinos on free neutrons and protons is turned off.
!
!  iaence   : emission and absorption of e-neutrinos on nuclei as given
!   by the FFN prescription  switch
!
!     iaence = 0: emission and absorption of e-neutrinos on nuclei as
!      given by the FFN prescription turned off.
!     iaence = 1: emission and absorption of e-neutrinos on nuclei as
!      given by the FFN prescription turned on.
!
!  edmpe    : difference between the mean excited state energy of the
!   daughter nucleus and that of the parent nucleus [MeV].
!
!  iaenca   : emission and absorption of e-antineutrinos on nuclei as
!   given by the FFN prescription  switch
!
!     iaenca = 0: emission and absorption of e-antineutrinos on nuclei as
!      given by the FFN prescription turned off.
!     iaenca = 1: emission and absorption of e-antineutrinos on nuclei as
!      given by the FFN prescription turned on.
!
!  edmpa    : difference between the mean excited state energy of the
!   daughter nucleus and that of  the parent nucleus [MeV].
!
!  iaenct: emission and absorption of e-neutrinos on nuclei as given by
!   the NSE-folded tabular data for Hix et al (2003) 
!
!     iaenct  = 0: emission and absorption of e-neutrinos on nuclei as
!      given by the tabular data is turned off.
!     iaenct  = 1: emission and absorption of e-neutrinos on nuclei as
!      given by the tabular data is turned on.
!
!  roaenct: density above which emission and absorption of e-neutrinos on
!   nuclei as given by Hix et al (2003) is turned off.
!-----------------------------------------------------------------------

INTEGER                                    :: iaefnp
INTEGER                                    :: i_aeps
INTEGER                                    :: iaence
INTEGER                                    :: iaenca
INTEGER                                    :: iaenct

REAL(double)                               :: rhoaefnp
REAL(double)                               :: edmpe
REAL(double)                               :: edmpa
REAL(double)                               :: roaenct

!-----------------------------------------------------------------------
!  Weak magnetism correction for charged current interactions
!-----------------------------------------------------------------------
!  icc_wm : charged current weak magnetism switch.
!
!     icc_wm = 0: weak magnetism correction for charged current
!      interactions turned off.
!     icc_wm = 1: weak magnetism correction for charged current
!      interactions turned on.
!-----------------------------------------------------------------------

INTEGER                                    :: icc_wm

!-----------------------------------------------------------------------
!  General scattering switch
!-----------------------------------------------------------------------
!  iscat : general scattering switch.
!
!     iscat = 0 : all scattering processes turned off.
!     iscat = 1 : any scattering process is allowed to proceed if its
!      key is turned on.                                           
!-----------------------------------------------------------------------

INTEGER                                    :: iscat

!-----------------------------------------------------------------------
!  Neutrino nucleon scattering controls
!-----------------------------------------------------------------------
!  in     : neutrino-neutron scattering switch
!
!     in = 0 : neutrino-neutron isoenergetic scattering turned off.
!     in = 1 : neutrino-neutron isoenergetic scattering turned on.
!
!  ip     : neutrino-proton scattering switch
!
!     ip = 0: neutrino-proton isoenergetic scattering turned off.
!     ip = 1: neutrino-proton isoenergetic scattering turned on.
!
!  ietann : final state blocking switch
!
!     ietann = 0 : final state blocking due to neutron or proton final
!      state occupancy in neutrino-neutron and neutrino-proton
!      isoenergetic scattering turned off.
!     ietann = 1 : final state blocking due to neutron or proton final
!      state occupancy in neutrino-neutron and neutrino-proton
!      isoenergetic scattering turned on.
!-----------------------------------------------------------------------

INTEGER                                    :: in
INTEGER                                    :: ip
INTEGER                                    :: ietann

!-----------------------------------------------------------------------
!  Neutrino nuclei scattering controls
!-----------------------------------------------------------------------
!  ihe    : neutrino-helium isoenergetic scattering switch
!
!     ihe = 0: neutrino-helium isoenergetic scattering turned off.
!     ihe = 1: neutrino-helium isoenergetic scattering turned on.
!
!  iheavy : neutrino-heavy nucleus isoenergetic scattering switch
!
!     iheavy = 0: neutrino-heavy nucleus isoenergetic scattering turned off.
!     iheavy = 1: neutrino-heavy nucleus isoenergetic scattering turned on.
!
!  iicor  : ion-ion correlation switch.
!
!     iicor = 0: neutrino-heavy nucleus isoenergetic scattering
!      correction due to ion-ion correlation effects turned off.
!     iicor = 1: neutrino-heavy nucleus isoenergetic scattering
!      correction due to ion-ion correlation effects turned on.
!-----------------------------------------------------------------------

INTEGER                                    :: ihe
INTEGER                                    :: iheavy
INTEGER                                    :: iicor

!-----------------------------------------------------------------------
!  Weak magnetism correction for neutral current interactions
!-----------------------------------------------------------------------
!  inc_wm : neutral current weak magnetism switch.
!
!     inc_wm = 0: weak magnetism correction for neutral current
!      interactions turned off.
!     inc_wm = 1: weak magnetism correction for neutral current
!      interactions turned on.
!-----------------------------------------------------------------------

INTEGER                                    :: inc_wm

!-----------------------------------------------------------------------
!  Neutrino-electron scattering controls
!-----------------------------------------------------------------------
!  nes      : neutrino-electron scattering switch.
!
!     nes = 0 : neutrino-electron scattering turned off.
!     nes = 1 : neutrino-electron scattering turned on.
!
!  rhonesmn : density below which n-neutrino-electron scattering is
!   turned off.
!
!  rhonesmx : density above which n-neutrino-electron scattering is
!   turned off.
!-----------------------------------------------------------------------

INTEGER                                    :: nes

REAL(double)                               :: rhonesmn
REAL(double)                               :: rhonesmx

!-----------------------------------------------------------------------
!  Neutrino-positron scattering controls
!-----------------------------------------------------------------------
!  nps      : neutrino-positron scattering switch.
!
!     nes = 0 : neutrino-positron scattering turned off.
!     nes = 1 : neutrino-positron scattering turned on.
!
!  rhonpsmn : density below which n-neutrino-positron scattering is
!              turned off.
!
!  rhonpsmx : density above which n-neutrino-positron scattering is
!              turned off.
!-----------------------------------------------------------------------

INTEGER                                    :: nps

REAL(double)                               :: rhonpsmn
REAL(double)                               :: rhonpsmx

!-----------------------------------------------------------------------
!  Neutrino-nucleon elastic scattering controls
!-----------------------------------------------------------------------
!  isctn      : neutrino-nucleon inelastic scattering switch.
!
!     isctn = 0: neutrino-nucleon inelastic scattering turned off.
!     isctn = 1: neutrino-nucleon inelastic scattering turned on.
!                Reddy et al.
!     isctn = 2: neutrino-nucleon inelastic scattering turned on.
!                Hannested & Raffelt
!
!  rhosctnemn : density below which n-neutrino-neutrino-nucleon,
!                elastic scattering is turned off.
!
!  rhosctnemx : density above which n-neutrino-neutrino-nucleon,
!                elastic scattering is turned off.
!-----------------------------------------------------------------------

INTEGER                                    :: isctn

REAL(double)                               :: rhosctnmn
REAL(double)                               :: rhosctnmx

!-----------------------------------------------------------------------
!  Neutrino-nucleus inelastic scattering controls
!-----------------------------------------------------------------------
!  isctnA      : neutrino-nuclei inelastic scattering switch.
!
!     isctnA = 0 : neutrino-nucleus inelastic scattering turned off.
!     isctnA = 1 : neutrino-nucleus inelastic scattering turned on.
!
!  rhosctnAemn : density below which e-neutrino-neutrino-nucleus,
!   e-antineutrino-neutrino-nucleus  inelastic scattering is turned off.
!
!  rhosctnAemx : density above which e-neutrino-neutrino-nucleus,
!   e-antineutrino-neutrino-nucleus inelastic scattering is turned off.
!
!  rhosctnAtmn : density below which x-neutrino-neutrino-nucleus,
!   x-antineutrino-neutrino-nucleus  inelastic scattering is turned off.
!
!  rhosctnAtmx : density above which x-neutrino-neutrino-nucleus,
!   x-antineutrino-neutrino-nucleus scattering is turned off.
!-----------------------------------------------------------------------

INTEGER                                    :: isctnA

REAL(double)                               :: rhosctnAemn
REAL(double)                               :: rhosctnAemx
REAL(double)                               :: rhosctnAtmn
REAL(double)                               :: rhosctnAtmx

!-----------------------------------------------------------------------
!  Electron-positron pair annihilation controls
!-----------------------------------------------------------------------
!  ipair      : electron-positron pair annihilation switch.
!
!     ipair = 0 : electron-positron pair annihilation process turned off.
!     ipair = 1 : electron-positron pair annihilation process turned on.
!
!  rhopairemn : density below which e-neutrino - e-antineutrino pair
!   production by electron-positron pair annihilation is turned off.
!
!  rhopairemx : density above which e-neutrino - e-antineutrino pair
!   production by electron-positron pair annihilation is turned off.
!
!  rhopairtmn : density below which x-neutrino - x-antineutrino pair
!   production by electron-positron pair annihilation is turned off..
!
!  rhopairtmx : density above which x-neutrino - x-antineutrino pair
!   production by electron-positron pair annihilation is turned off.
!-----------------------------------------------------------------------

INTEGER                                    :: ipair

REAL(double)                               :: rhopairemn
REAL(double)                               :: rhopairemx
REAL(double)                               :: rhopairtmn
REAL(double)                               :: rhopairtmx

!-----------------------------------------------------------------------
!  Bremsstrahlung pair annihilation controls
!-----------------------------------------------------------------------
!  ibrem      : bremsstrahlung pair annihilation switch.
!
!     ibrem = 0 : bremsstrahlung pair annihilation process turned off.
!     ibrem = 1 : bremsstrahlung pair annihilation process turned on.
!
!  rhobrememn : density below which e-neutrino - e-antineutrino pair
!   production by bremsstrahlung pair annihilation is turned off.
!
!  rhobrememx : density above which e-neutrino - e-antineutrino pair
!   production by bremsstrahlung pair annihilation is turned off.
!-----------------------------------------------------------------------

INTEGER                                    :: ibrem

REAL(double)                               :: rhobrememn
REAL(double)                               :: rhobrememx

!-----------------------------------------------------------------------
!  NSE flashing and deflashing controls
!-----------------------------------------------------------------------
!  tnse  : temperature at which material is flashed to nse [K]
!
!  tdnse : temperature at which material is deflashed from nse [K]
!
!  deflash_pre_bounce : deflashing permitted before bounce
!-----------------------------------------------------------------------

REAL(double)                               :: tnse
REAL(double)                               :: tdnse

LOGICAL                                    :: deflash_pre_bounce

!-----------------------------------------------------------------------
!  Time window averaging of the GR pseudo-potential
!-----------------------------------------------------------------------
!  t_wndw_grpot  : time window for averaging the GR pseudo-potential
!   [micro-seconds]
!
!  t_apply_gr_ave: time after bounce at which time-window averaging
!   is to be applied [milli-seconds]
!-----------------------------------------------------------------------

REAL(double)                               :: t_wndw_grpot
REAL(double)                               :: t_apply_gr_ave

!-----------------------------------------------------------------------
!  ave_pole: pole average width
!-----------------------------------------------------------------------

INTEGER                                    :: ave_pole

END module prb_cntl_module
