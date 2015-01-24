!-----------------------------------------------------------------------
!    Module:       eos_bck_module
!    Author:       S. W. Bruenn
!    Date:         8/16/02
!-----------------------------------------------------------------------

MODULE eos_bck_module

USE wlKindModule, ONLY: dp

!-----------------------------------------------------------------------
!  Quantities used in the BCK EOS
!-----------------------------------------------------------------------
!  bad    : auxilary variable.
!
!  fmbad  : auxilary variable.
!
!  itsav  : the number of iterations done in SAHA.
!
!  af     : auxilary variable.
!
!  ahbck  : A[heavy nucleus].
!
!  b_hvy  : mean binding energy per particle of the of the heavies
!            excluding n, p, and He).
!
!  b      : mean binding energy per particle.
!
!  con    : auxilary variable.
!
!  dbck   : density (baryons/fm**3).
!
!  d0     : saturation density of dense phase, or density of nuclear matter
!   (baryons/fm**3).
!
!  d00    : zero value of saturation density of dense phase, or density of
!   Fe56 (baryons/fm**3).
!
!  dbdlu  : auxilary variable.
!
!  dedt - dE / dT  [approximate!].
!
!  delbet : ( ue - uneu - uhat ) [affinity] : not functional.
!
!  dtran  : transition density to nuclear matter.
!
!  ed     : energy in drip (including heavy translation).
!
!  ee     : energy in electrons.
!
!  egy0   : zero value of energy of Fe56 [MeV].
!
!  eh     : energy in heavy nuclei (or nuclear matter) w/o tranlation.
!
!  enu    : energy in neutrinos ; not calculated.
!
!  erad   : photon energy.
!
!  etot   : total energy ; no neutrinos ; includes radiation.
!
!  fa     : auxilary variable.
!
!  fm     : auxilary variable.
!
!  pd     : drip pressure (including heavy translation).
!
!  pe     : electron pressure.
!
!  ph     : heavy pressure (no translation).
!
!  phi    : auxilary variable.
!
!  pnu    : neutrino pressure ; not calculated.
!
!  prad   : photon pressure
!
!  ptotbck: total pressure (includes radiation).
!
!  rel    : 3 * pressure / kinetic energy (from 1(rel) to 2(nonrel)).
!
!  sd     : entropy in drip (includes heavy translation).
!
!  se     : entropy in electrons.
!
!  sh     : entropy in heavy nuclei (intrinsic only).
!
!  size0  : zero value of size of Fe56.
!
!  sneu   : neutrino entropy (not calculated).
!
!  stotbck: total entropy ; includes radiation.
!
!  tbck   : temperature [MeV].
!
!  therm  : auxilary variable.
!
!  theta  : density of dense phase / saturation density of dense phase.
!
!  theta0 : zero value of Fe56 density / saturation density of Fe56.
!
!  tsi    : auxilary variable.
!
!  upack  : packing fraction of dense phase.
!
!  u1pack : 1 - upack.
!
!  ue     : electron chemical potential (mu_e).
!
!  uhat   :  mu_n - mu_p.
!
!  uhtra  :  auxiliary variable
!
!  un     :  neutron chemical potential (mu_n).
!
!  uneu   :  neutrino chemical potential (mu_nu) ; not calculated.
!
!  xabck  :  alpha particle mass fraction ; = 0 for nuclear matter.
!
!  xhbck  :  heavy nuclei mass fraction ; = 1 for nuclear matter.
!
!  xmstar :  m* / m for heavy nuclei or nuclear matter.
!
!  xnbck  :  free neutron mass fraction ; = 0 for nuclear matter.
!
!  xpbck  :  free proton mass fraction  ; = 0 for nuclear matter.
!
!  ya     :  auxilary variable.
!
!  yeFe   :  zero value of Y_e^- - Y_e^+ for Fe56.
!
!  yebck  :  Y_e^- - Y_e^+.
!
!  yeplus :  Y_e^+.
!
!  yh     :  auxilary variable.
!
!  ynbck  :  auxilary variable.
!
!  ynu    :  Y_nu ; not used.
!
!  ypbck  :  auxilary variable.
!
!  zabck  :  Z / A for heavy nucleus ; = ye for nuclear matter. 
!
!-----------------------------------------------------------------------

LOGICAL                               :: bad
LOGICAL                               :: fmbad

INTEGER                               :: itsav

INTEGER, PARAMETER                    :: nnc_nse_bck = 17

REAL(dp)                              :: af
REAL(dp)                              :: ahbck
REAL(dp)                              :: b
REAL(dp)                              :: b_hvy
REAL(dp)                              :: con
REAL(dp)                              :: dbck
REAL(dp)                              :: d0
REAL(dp)                              :: d00
REAL(dp)                              :: dbdlu
REAL(dp)                              :: dedt
REAL(dp)                              :: delbet
REAL(dp)                              :: dtran  = 0.0d0
REAL(dp)                              :: ed
REAL(dp)                              :: ee
REAL(dp)                              :: egy0
REAL(dp)                              :: eh
REAL(dp)                              :: enu
REAL(dp)                              :: erad
REAL(dp)                              :: etot
REAL(dp)                              :: fa
REAL(dp)                              :: fm
REAL(dp)                              :: pd
REAL(dp)                              :: pe
REAL(dp)                              :: ph
REAL(dp)                              :: phi
REAL(dp)                              :: pnu
REAL(dp)                              :: prad
REAL(dp)                              :: ptotbck
REAL(dp)                              :: rel
REAL(dp)                              :: sd
REAL(dp)                              :: se
REAL(dp)                              :: sh
REAL(dp)                              :: size0
REAL(dp)                              :: sneu
REAL(dp)                              :: stotbck
REAL(dp)                              :: tbck
REAL(dp)                              :: therm
REAL(dp)                              :: theta
REAL(dp)                              :: theta0
REAL(dp)                              :: tsi
REAL(dp)                              :: upack
REAL(dp)                              :: u1pack
REAL(dp)                              :: ue
REAL(dp)                              :: uhat
REAL(dp)                              :: uhtra
REAL(dp)                              :: un
REAL(dp)                              :: uneu
REAL(dp)                              :: xabck
REAL(dp)                              :: xhbck
REAL(dp)                              :: xmstar
REAL(dp)                              :: xnbck
REAL(dp)                              :: xpbck
REAL(dp)                              :: ya
REAL(dp)                              :: yeFe
REAL(dp)                              :: yebck
REAL(dp)                              :: yeplus
REAL(dp)                              :: yh
REAL(dp)                              :: ynbck
REAL(dp)                              :: ynu
REAL(dp)                              :: ypbck
REAL(dp)                              :: zabck

!-----------------------------------------------------------------------
!  eos_gamma:  value of gamma for the gamma-law equation of state
!-----------------------------------------------------------------------

REAL(dp)                              :: eos_gamma

!-----------------------------------------------------------------------
!  ye_nnc_nse_min:  lower bound of ye for use of the 17 nuclear specie
!                    system in NSE for nse = 1 and for 'B' EOS
!-----------------------------------------------------------------------

REAL(dp), PARAMETER                   :: ye_nnc_nse_min = 26.d0/56.d0

END module eos_bck_module
