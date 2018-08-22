!-----------------------------------------------------------------------
!    Module:       edit_module
!    Author:       S. W. Bruenn
!    Date:         8/16/02
!-----------------------------------------------------------------------

MODULE edit_module

USE kind_module

!-----------------------------------------------------------------------
!
!               \\\\\ GENERAL EDIT PARAMETERS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  data_path     : the path directing output to the data directories.
!
!  log_path      : path to the simulation log
!
!  reset_path    : path to write the restart key file, reset.d
!
!  nread         : the unit number for file 'reset.d', where the initial
!   model or the information directing input to the restart file is stored.
!
!  nrrst         : the unit number of the directory that subroutine readst
!   is directed to look for the initial model or the information directing
!   input to the restart file. Typically nnrst = nread.
!
!  iprint        : the parameter used in the call to an edit subroutine;
!
!     0 : print to a print file is bypassed.
!     1 : print to a print file is implemented.
!
!  nprint        : a parameter used in the call to an edit subroutine
!   denoting the unit number of the file to which the print file is to be
!   sent. It is usually redefined before a call but nominally carries the
!   value denoting the unit number of the file 'superdump' to which any
!   diagnostics are written during the course of a run.
!
!  nlog          : the unit number of the file to which the simulation log
!   is directed
!
!  prnttest      : a logical parameter used in the call to an edit
!   subroutine. Calls to any edit subroutine can be indexed so that
!   particular sections of the subroutine are activated only after a
!   prescribed number of calls. Edit subroutines are routinely called
!   after each hydro cycle, and this permits each section of an edit
!   subroutine to print only after a prescribed number of cycles have
!   elapsed since it last printed.
!
!     T : call to the edit subroutine is indexed.
!     F : call to the edit subroutine is not indexed.
!
!  intprt        : the number of cycles between the closing of an old print
!   file and the opening of a new print file. The print files are named
!   modelxxx.d, where xxx is the right justified consecutive number of
!   the print file, e.g., the fifth print file is named 'model005.d'.
!
!  nprt(i_ray)   : a print counter giving the number of cycles since the
!   last closing of a print file.
!
!  iprtbgn       : a print parameter regulating the print file created
!   when the calculation is initiated (cycle_flag = 0) or restarted
!   (cycle_flag ne 0);
!
!     iprtbgn = 0, cycle_flag = 0 : abbreviated edit
!     iprtbgn = 0, cycle_flag > 0 : configuration edit
!     iprtbgn = 1, cycle_flag = 0 : abbreviated edit
!     iprtbgn = 1, cycle_flag > 0 : abbreviated edit
!     iprtbgn = 2, cycle_flag = 0 : full edit
!     iprtbgn = 2, cycle_flag > 0 : configuration edit
!     iprtbgn = 3, cycle_flag = 0 : full edit
!     iprtbgn = 3, cycle_flag > 0 : abreviated edit
!     iprtbgn = 4, cycle_flag = 0 : full edit
!     iprtbgn = 4, cycle_flag > 0 : full edit
!
!  iflprt        : a print parameter regulating the print file created when the
!   calculation is terminated (ncycle = ncymax)
!
!     iflprt = 0:  no printout at termination
!     iflprt = 1:  configuration printout at termination
!     iflprt = 2:  full printout except editn at termination
!     iflprt = 3:  full printout except editng at termination
!     iflprt = 4:  full printout at termination
!
!  nmodel(i_ray) : the number of the current print file.
!
!  head : the character string bearing the name of the calculation. It is
!   printed at the beginning of each edit subfile.
!
!  rhoprint(i)   : the ith central density at which a print file is created.
!   The print file is named modeldxx.d, where xx is the value of i right
!   justified.
!
!  itprnt : elapsed time - post-bounce time flag
!
!  tprint : the ith elapsed time [ms] (itprnt = 1) or post-bounce
!   time [ms] (itprnt = 2) at which a print file is created. The print
!   file is named
!       modeleltxx_yy.d    (itprnt = 1)
!   or
!       modelpbtxx_yy.d    (itprnt = 1)
!      
!   where xx is the processor number and yy value of ith post-bounce
!   time in ms right justified.
!
!  iplot         : a parameter used in the call to an edit subroutine -
!
!     0 : print to a plot file is bypassed.
!     1 : print to a plot file is implemented.
!
!  nplot         : a parameter used in the call to an edit subroutine
!   denoting the unit number of the file to which the plot file is to be
!   sent.
!
!  intplf        : the number of cycles between the printing to a plot
!    file. 
!
!  npltf         : the number of cycles since the last printing to a
!   plot file.
!
!  nplotc        : the unit number for downloading model configuration
!   plot data.
!
!  nplote        : the unit number for downloading electron neutrino
!   plot data.
!
!  nplota        : the unit number for downloading electron antineutrino
!   plot data.
!
!  nplott        : the unit number for downloading muon and tau neutrino
!   plot data.
!
!  i_edit        : edit parameter,
!     O - edits are executed on the basis of the counter values
!     1 - brief model edit
!     2 - full model edit
!     3 - full model edit plus differential neutrino edit
!     4 - full model edit plus differential plus integrated neutrino edit
!
!  it_edit       : key for editing at selected postbounce time intervals.
!
!  dt_edit       : postbounce time intervals for editing (ms).
!-----------------------------------------------------------------------

LOGICAL                                             :: prnttest

CHARACTER(len = 128)                                :: data_path
CHARACTER(len = 128)                                :: log_path
CHARACTER(len = 128)                                :: reset_path
CHARACTER(len = 128)                                :: head

INTEGER                                             :: nread
INTEGER                                             :: nrrst
INTEGER                                             :: iprint
INTEGER                                             :: nprint
INTEGER                                             :: nlog
INTEGER                                             :: intprt
INTEGER, ALLOCATABLE, DIMENSION(:)                  :: nprt
INTEGER                                             :: iprtbgn
INTEGER                                             :: iflprt
INTEGER, ALLOCATABLE, DIMENSION(:)                  :: nmodel
INTEGER                                             :: iplot
INTEGER                                             :: nplot
INTEGER                                             :: intplf
INTEGER                                             :: npltf
INTEGER                                             :: nplotc
INTEGER                                             :: nplote
INTEGER                                             :: nplota
INTEGER                                             :: nplott
INTEGER                                             :: i_edit
INTEGER                                             :: it_edit
INTEGER                                             :: itprnt

REAL(KIND=double), DIMENSION(10)                    :: rhoprint
REAL(KIND=double), DIMENSION(50)                    :: tprint
REAL(KIND=double)                                   :: dt_edit

!-----------------------------------------------------------------------
!  Parameters for configuration edits
!-----------------------------------------------------------------------
!  intedc(i) : the number of cycles between the implementation of
!   subsection i of subroutine editc. 
!
!  nedc(i)   : the number of cycles since the last implementation of
!   subsection i of subroutine editc.
!
!  idxedc(i) : a parameter used in subroutine editc. Subsection i of
!   subroutine editc will print data for every idxedc(i)'th radial zone,
!   starting with the outermost zone. The data corresponding to the
!   innermost radial zone will also be printed. Setting idxedc(i) = 1 will
!   cause data for all radial zones to be printed.
!-----------------------------------------------------------------------

INTEGER, DIMENSION(20)                              :: intedc
INTEGER, DIMENSION(20)                              :: nedc
INTEGER, DIMENSION(20)                              :: idxedc

!-----------------------------------------------------------------------
!  Parameters for kinetic, internal, and gravitational energy edits
!-----------------------------------------------------------------------
!  intede(i) : the number of cycles between the implementation of
!   subsection i of subroutine edit_e. 
!
!  nede(i)   : the number of cycles since the last implementation of
!   subsection i of subroutine edit_e.
!
!  idxede(i) : a parameter used in subroutine edit_e. Subsection i of
!   subroutine editc will print data for every idxede(i)'th radial zone,
!   starting with the outermost zone. The data corresponding to the
!   innermost radial zone will also be printed. Setting idxede(i) = 1 will
!   cause data for all radial zones to be printed.
!-----------------------------------------------------------------------

INTEGER, DIMENSION(20)                              :: intede
INTEGER, DIMENSION(20)                              :: nede
INTEGER, DIMENSION(20)                              :: idxede

!-----------------------------------------------------------------------
!  Parameters for editmi edits
!-----------------------------------------------------------------------
!  intdmi(i) : the number of cycles between the implementation of
!   subsection i of subroutine editmi.
!
!  nedmi(i)  : the number of cycles since the last implementation of
!   subsection i of subroutine editmi.
!
!  idxemi(i) : a parameter used in subroutine editmi. Subsection i of
!   subroutine editmi will print
!   data for every idxemi(i)'th radial zone, starting with the outermost
!   zone. The data corresponding to the innermost radial zone will also
!   be printed. Setting idxemi(i) = 1 will cause data for all radial
!   zones to be printed.
!-----------------------------------------------------------------------

INTEGER, DIMENSION(20)                              :: intdmi
INTEGER, DIMENSION(20)                              :: nedmi
INTEGER, DIMENSION(20)                              :: idxemi

!-----------------------------------------------------------------------
!  Parameters for mass average edits
!-----------------------------------------------------------------------
!  intdma(i) : the number of cycles between the implementation of
!   subsection i of subroutine editma.
!
!  nedma(i)  : the number of cycles since the last implementation of
!   subsection i of subroutine editma.
!
!  idxema(i) : a parameter used in subroutine editma. Subsection i of
!   subroutine editma will print data for every idxema(i)'th radial zone,
!   starting with the outermost zone. The data corresponding to the
!   innermost radial zone will also be printed. Setting idxema(i) = 1 will
!   cause data for all radial zones to be printed.
!-----------------------------------------------------------------------

INTEGER, DIMENSION(20)                              :: intdma
INTEGER, DIMENSION(20)                              :: nedma
INTEGER, DIMENSION(20)                              :: idxema

!-----------------------------------------------------------------------
!  Parameters for hydrodynamic edits
!-----------------------------------------------------------------------
!  intedh(i) : the number of cycles between the implementation of
!   subsection i of subroutine edith.
!
!  nedh(i)   : the number of cycles since the last implementation of
!   subsection i of subroutine edith.
!
!  idxedh(i) : a parameter used in subroutine edith. Subsection i of
!   subroutine editma will print data for every idxedh(i)'th radial zone,
!   starting with the outermost zone. The data corresponding to the
!   innermost radial zone will also be printed. Setting idxedh(i) = 1 will
!   cause data for all radial zones to be printed.
!-----------------------------------------------------------------------

INTEGER, DIMENSION(20)                              :: intedh
INTEGER, DIMENSION(20)                              :: nedh
INTEGER, DIMENSION(20)                              :: idxedh

!-----------------------------------------------------------------------
!  Parameters for hydrodynamic edits
!-----------------------------------------------------------------------
!  intedh(i) : the number of cycles between the implementation of
!   subsection i of subroutine edith.
!
!  nedh(i)   : the number of cycles since the last implementation of
!   subsection i of subroutine edith.
!
!  idxedh(i) : a parameter used in subroutine edith. Subsection i of
!   subroutine editma will print data for every idxedh(i)'th radial zone,
!   starting with the outermost zone. The data corresponding to the
!   innermost radial zone will also be printed. Setting idxedh(i) = 1 will
!   cause data for all radial zones to be printed.
!-----------------------------------------------------------------------

INTEGER, DIMENSION(20)                              :: intdps
INTEGER, DIMENSION(20)                              :: nedps
INTEGER, DIMENSION(20)                              :: idxeps

!-----------------------------------------------------------------------
!  Parameters for energy edits
!-----------------------------------------------------------------------
!  intedu(i) : the number of cycles between the implementation of
!   subsection i of subroutine editu.
!
!  nedu(i)   : the number of cycles since the last implementation of
!   subsection i of subroutine editu.
!
!  idxedu(i) : a parameter used in subroutine editu. Subsection i of
!   subroutine editu will print data for every idxedu(i)'th radial zone,
!   starting with the outermost zone. The data corresponding to the
!   innermost radial zone will also be printed. Setting idxedu(i) = 1 will
!   cause data for all radial zones to be printed.
!-----------------------------------------------------------------------

INTEGER, DIMENSION(20)                              :: intedu
INTEGER, DIMENSION(20)                              :: nedu
INTEGER, DIMENSION(20)                              :: idxedu

!-----------------------------------------------------------------------
!  Parameters for composition edits
!-----------------------------------------------------------------------
!  intedy(i): the number of cycles between the implementation of
!   subsection i of subroutine edity.
!
!  nedy(i): the number of cycles since the last implementation of
!   subsection i of subroutine edity.
!
!  idxedy(i): a parameter used in subroutine edity. Subsection i of
!   subroutine edity will print data for every idxedy(i)'th radial zone,
!   starting with the outermost zone. The data corresponding to the
!   innermost radial zone will also be printed. Setting idxedy(i) = 1 will
!   cause data for all radial zones to be printed.
!-----------------------------------------------------------------------

INTEGER, DIMENSION(20)                              :: intedy
INTEGER, DIMENSION(20)                              :: nedy
INTEGER, DIMENSION(20)                              :: idxedy

!-----------------------------------------------------------------------
!  Parameters for entropy and chemical potential edits
!-----------------------------------------------------------------------
!  intdsc(i) : the number of cycles between the implementation of
!   subsection i of subroutine editsc.
!
!  nedsc(i)  : the number of cycles since the last implementation of
!   subsection i of subroutine editsc.
!
!  idxesc(i) : a parameter used in subroutine editsc. Subsection i of
!   subroutine editsc will print data for every idxesc(i)'th radial zone,
!   starting with the outermost zone. The data corresponding to the
!   innermost radial zone will also be printed. Setting idxesc(i) = 1 will
!   cause data for all radial zones to be printed.
!-----------------------------------------------------------------------

INTEGER, DIMENSION(20)                              :: intdsc
INTEGER, DIMENSION(20)                              :: nedsc
INTEGER, DIMENSION(20)                              :: idxesc

!-----------------------------------------------------------------------
!  Parameters for diferential neutrino edits
!-----------------------------------------------------------------------
!  intedn(n) : the number of cycles between the implementation of
!   subroutine editn for neutrinos of type n.
!
!  nedn(n)   : the number of cycles since the last implementation of
!   subroutine editn for neutrinos of type n.
!
!  idxedn(n) : a parameter used in subroutine editn. Subroutine editn
!   will print the data corresponding neutrinos of type n for every
!   idxedn(n)'th radial zone, starting with the outermost zone. The data
!   corresponding to the innermost radial zone will also be printed.
!   Setting idxedn(i) = 1 will cause data for all radial zones to be
!   printed.
!
!  niedn(i)  : a parameter used in subroutine editn. Subroutine editn
!   will print the i'th datum corresponding to neutrinos of type n every
!   niedn(i) implementation of editn.
!
!  neden(i)  : the number of implementations of editn since the last
!   printing of the i'th datum corresponding to neutrinos of type n.
!-----------------------------------------------------------------------

INTEGER, ALLOCATABLE, DIMENSION(:)                  :: intedn
INTEGER, ALLOCATABLE, DIMENSION(:)                  :: nedn
INTEGER, ALLOCATABLE, DIMENSION(:)                  :: idxedn
INTEGER, DIMENSION(60)                              :: niedn
INTEGER, DIMENSION(60)                              :: neden

!-----------------------------------------------------------------------
!  Parameters for integral neutrino edits
!-----------------------------------------------------------------------
!  intdng(i,n) : the number of cycles between the implementation of
!   subsection i of subroutine editng for neutrinos of type n.
!
!  nedng(i,n)  : the number of cycles since the last implementation
!   of subsection i of subroutine editng for neutrinos of type n.
!
!  idxeng(i,n) : a parameter used in subroutine editng. Subsection i of
!   subroutine editng will print the data corresponding to every
!   idxesc(i)'th radial zone, starting with the outermost zone. The
!   data corresponding to the innermost radial zone will also be printed.
!   Setting idxeng(i,n) = 1 will cause data for all radial zones to be
!   printed. 
!-----------------------------------------------------------------------

INTEGER, ALLOCATABLE, DIMENSION(:,:)                :: intdng
INTEGER, ALLOCATABLE, DIMENSION(:,:)                :: nedng
INTEGER, ALLOCATABLE, DIMENSION(:,:)                :: idxeng

!-----------------------------------------------------------------------
!  Parameters for neutrino luminosity plot edits
!-----------------------------------------------------------------------
!  ilumplt  : luminosity plot edit switch
!
!          0: bypass subroutine lumplot.
!          1: enter subroutine lumplot.
!
!  nlumplt1 : the unit number for editing lumplot data at r_lum.
!
!  nlumplt2 : the unit number for editing lumplot data at d_lum.
!
!  nlum     : number of cycles since the last implemtation of subroutine lumplot.
!
!  r_lum    : the radius at which to evaluate luminosity data.
!
!  d_lum    : the density at which to evaluate luminosity data.
!
!  intlum, ncylum:
!
!     ncycle < ncylum(1)             : subroutine lumplot is called every
!                                       intlum(1) cycles.
!     ncylum(1) < ncycle < ncylum(2) : subroutine lumplot is called every
!                                       intlum(2) cycles.
!     ncylum(2) < ncycle             : subroutine lumplot is called every
!                                       intlum(3) cycles.                                                     !-----------------------------------------------------------------------

INTEGER                                             :: ilumplt
INTEGER                                             :: nlumplt1
INTEGER                                             :: nlumplt2
INTEGER                                             :: nlum
INTEGER, DIMENSION(3)                               :: intlum
INTEGER, DIMENSION(2)                               :: ncylum

REAL(KIND=double)                                   :: r_lum
REAL(KIND=double)                                   :: d_lum

!-----------------------------------------------------------------------
!  Parameters for rms neutrino energy plot edits
!-----------------------------------------------------------------------
!  ienuplt : rms neutrino energy plot edit switch
!
!          0: bypass subroutine enuvplot.
!          1: enter subroutine enuvplot.
!
!  nenuplt : the unit number for editing enuvplot data.
!
!  nenu    : the number of cycles since the last implemtation of subroutine
!   enuvplot.
!
!  r_e_rms :  the radius at which to store neutrino energy data.
!
!  d_e_rms : the density at which to store neutrino energy data.
!
!  intenu, ncyenu:
!
!     ncycle < ncyenu(1)             : subroutine enuvplot is called every                                          !                                       intenu(1) cycles.
!     ncyenu(1) < ncycle < ncyenu(2) : subroutine enuvplot is called every
!                                       intenu(2) cycles.
!     ncyenu(2) < ncycle             : subroutine enuvplot is called every
!                                       intenu(3) cycles.
!-----------------------------------------------------------------------

INTEGER                                             :: ienuplt
INTEGER                                             :: nenuplt1
INTEGER                                             :: nenuplt2
INTEGER                                             :: nenu
INTEGER, DIMENSION(3)                               :: intenu
INTEGER, DIMENSION(2)                               :: ncyenu

REAL(KIND=double)                                   :: r_e_rms
REAL(KIND=double)                                   :: d_e_rms

!-----------------------------------------------------------------------
!  Parameters for selected radii plot edits
!-----------------------------------------------------------------------
!  irnuplt : selected radii plot edit switch
!
!          0: bypass subroutine rnuplot.
!          1: enter subroutine rnuplot.
!
!  nrnuplt : the unit number for editing rnuplot data.
!
!  nrum    : the number of cycles since the last implemtation of subroutine
!   rnuplot.
!
!  intrnu, ncyrnu:
!
!     ncycle < ncyrnu(1)             : subroutine rnuplot is called every
!                                       intrnu(1) cycles.
!     ncyrnu(1) < ncycle < ncyrnu(2) : subroutine rnuplot is called every
!                                       intrnu(2) cycles.      
!     ncyrnu(2) < ncycle             : subroutine rnuplot is called every
!                                       intrnu(3) cycles.
!-----------------------------------------------------------------------

INTEGER                                             :: irnuplt
INTEGER                                             :: nrnuplt
INTEGER                                             :: nrnu
INTEGER, DIMENSION(3)                               :: intrnu
INTEGER, DIMENSION(2)                               :: ncyrnu

!-----------------------------------------------------------------------
!  Parameters for configuration plot edits
!-----------------------------------------------------------------------
!  ivarplt   : configuration plot edit switch
!
!     0 : bypass subroutine varplot.
!     1 : enter subroutine varplot.
!
!  nvar      : the number of cycles since the last implemtation of subroutine varplot.                                         !
!
!  nvarint   : subroutine varplot is called every nvarint cycles.
!
!  nvarplt   : the unit number for editing varplot data.
!
!  nvardump  : varplot file number (varplot files are numbered sequentially.                                               !
!  dtvarplot : write to subroutine varplot every dtvarplot ms.
!-----------------------------------------------------------------------

INTEGER                                             :: ivarplt
INTEGER                                             :: nvar
INTEGER                                             :: nvarint
INTEGER                                             :: nvarplt
INTEGER                                             :: nvardump

REAL(KIND=double)                                   :: dtvarplot

!-----------------------------------------------------------------------
!  Parameters for comparison plot edits
!-----------------------------------------------------------------------
!  icomplt   : comparison plot edit switch
!
!     0 : bypass subroutine complot.
!     1 : enter subroutine complot.
!
!  ncomplt   : the unit number for editing complot data.
!
!  ncomdump  : the complot file number (complot files are numbered
!   sequentially.                                               !
!
!  dtcomplot : write to subroutine complot every dtcomplot ms.
!-----------------------------------------------------------------------

INTEGER                                             :: icomplt
INTEGER                                             :: ncomplt
INTEGER                                             :: ncomdump

REAL(KIND=double)                                   :: dtcomplot


!-----------------------------------------------------------------------
!  Parameters for bounday plot edits
!-----------------------------------------------------------------------
!  dtimeplot   : write to subroutine bnuplot every dtimeplot ms.
!
!  rinnerb     : quantities are interpolated to rinnerb for ibound data.
!
!  nplotinnerb : the unit number for ibound data.
!
!  iplotinnerb : inner boundary plot switch
!
!     0    : bypass writing ibound data.
!     ne 0 : write ibound data.
!
!  routerb     : quantities are interpolated to routerb for obound data.
!
!  nplotouterb : the unit number for obound data.
!
!  iplotouterb : outer boundary plot switch
!
!     0   : bypass writing obound data.
!     ne 0: write obound data.
!
!  r_lumerms   : luminosities and rms energies are interpolated to
!   r_lumerms for lbound data.                                            
!
!  nplotlum    : the unit number for lbound data.
!
!  iplotlum    : luminosities and rms energies plot switch
!
!     0   : bypass writing lbound data.
!     ne 0: write lbound data.
!
!  nplotshk    : the unit number for shock data.
!
!  iplotshk    : shock data plot switch
!
!     0    : bypass writing shock data.
!     ne 0 : write shock data.
!
!  nplotcnv    : the unit number for convection data.
!
!  iplotcnv    : convection plot switch
!
!     0    : bypass writing convection data.
!     ne 0 : write convection data.
!
!  nplotmss    : the unit number for mass data.
!
!  iplotmss    : mass data plot switch
!
!     0    : bypass writing mass data.
!     ne 0 : write mass data.
!-----------------------------------------------------------------------

INTEGER                                             :: nplotinnerb
INTEGER                                             :: iplotinnerb
INTEGER                                             :: nplotouterb
INTEGER                                             :: iplotouterb
INTEGER                                             :: nplotlum
INTEGER                                             :: iplotlum
INTEGER                                             :: nplotshk
INTEGER                                             :: iplotshk
INTEGER                                             :: nplotcnv
INTEGER                                             :: iplotcnv
INTEGER                                             :: nplotmss
INTEGER                                             :: iplotmss

REAL(KIND=double)                                   :: dtimeplot
REAL(KIND=double)                                   :: rinnerb
REAL(KIND=double)                                   :: routerb
REAL(KIND=double)                                   :: r_lumerms

!-----------------------------------------------------------------------
!  Parameters for nurad plot edits
!-----------------------------------------------------------------------
!  nplotnurad    : the unit number for nuradplot data.
!
!  iplotnurad    : nurad plot switch
!
!     0    : bypass writing nuradplot data.
!     ne 0 : write nuradplot data.
!
!  dtnuradplot   : write to subroutine nuradplot every dtnuradplot ms.
!
!  r_nurad       : the radius at which to evaluate neutrino radiation.
!
!  rho_nurad     : the density at which to evaluate neutrino radiation.
!
!  nu_r(k,n,i_ray)     : the number of neutrinos of energy group k
!   radiated across r_nurad in time dtnuradplot.
!
!  nu_rt(k,n,i_ray)    : the cumulative number of neutrinos of energy
!   group k radiated across r_nurad.
!
!  nu_rho(k,n,i_ray)   : the number of neutrinos of energy group k
!   radiated  across rho_nurad in time dtnuradplot.
!
!  nu_rhot(k,n,i_ray)  : the cumulative number of neutrinos of energy
!   group k radiated across rho_nurad.
!
!  unu_r(k,n,i_ray)    : the energy of neutrinos of energy group k
!   radiated across r_nurad in time dtnuradplot.
!
!  unu_rt(k,n,i_ray)   : the- cumulative energy of neutrinos of energy
!   group k radiated across r_nurad.
!
!  unu_rho(k,n,i_ray)  : the energy of neutrinos of energy group k
!   radiated across rho_nurad in time dtnuradplot.
!
!  unu_rhot(k,n,i_ray) : the cumulative energy of neutrinos of energy
!   group k radiated across rho_nurad.
!-----------------------------------------------------------------------

INTEGER                                             :: nplotnurad
INTEGER                                             :: iplotnurad

REAL(KIND=double)                                   :: dtnuradplot
REAL(KIND=double)                                   :: r_nurad
REAL(KIND=double)                                   :: rho_nurad
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)    :: nu_r
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)    :: nu_rt
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)    :: nu_rho
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)    :: nu_rhot
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)    :: unu_r
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)    :: unu_rt
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)    :: unu_rho
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)    :: unu_rhot

!-----------------------------------------------------------------------
!  Parameters for neutrino distribution plot edits
!-----------------------------------------------------------------------
!  nnudata      : the unit number for nuplot data.
!
!  inudata      : neutrino distribution plot switch
!
!     0    : bypass writing nuplot data.
!     ne 0 : write nuplot data.
!
!  dtnudata     : write to subroutine nuplot every dtnudata ms.
!
!  r_nudata     : the radius at which psi0 and psi1 are evaluated for nuplot.
!
!  t_nudata     : the time when the last data dump to nuplot was made.
!
!  psi0dat(k,n,i_ray) : the time integrated psi0, to be time averaged
!   on dumping.
!
!  psi1dat(k,n,i_ray) : the time integrated psi1, to be time averaged
!   on dumping.
!-----------------------------------------------------------------------

INTEGER                                             :: nnudata
INTEGER                                             :: inudata

REAL(KIND=double)                                   :: dtnudata
REAL(KIND=double)                                   :: r_nudata
REAL(KIND=double)                                   :: t_nudata
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)    :: psi0dat
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)    :: psi1dat

!-----------------------------------------------------------------------
!  Parameters for lagrangian plot edits
!-----------------------------------------------------------------------
!  nlagplt    : the unit number for lagplot data.
!
!  ilagplt    : lagrangian plot switch
!
!     0    : bypass writing lagplot data.
!     ne 0 : write nuplot data.
!
!  nlagdump   : lagplot file number (lagplot files are numbered sequentially.
!
!  msslag     : lagrangian mass of fluid element for lagplot data.
!
!  int_lagplt : the number of cycles between lagplt edits
!
!  ned_lagplt : the number of cycles since the last lagplt edit
!-----------------------------------------------------------------------

INTEGER                                             :: nlagplt
INTEGER                                             :: ilagplt
INTEGER                                             :: nlagdump
INTEGER                                             :: int_lagplt
INTEGER                                             :: ned_lagplt

REAL(KIND=double)                                   :: msslag

!-----------------------------------------------------------------------
!  Parameters for r-lagrangian plot edits
!-----------------------------------------------------------------------
!  nrlagplt : the unit number for rlagplot data.
!
!  irlagplt : r-lagrangian plot switch
!
!     0    : bypass writing rlagplot data.
!     ne 0 : write rlagplot data.
!
!  dmlag    : the lagrangian mass difference between adjacent lagrangian
!   points to be plotted.
!-----------------------------------------------------------------------

INTEGER                                             :: nrlagplt
INTEGER                                             :: irlagplt

REAL(KIND=double)                                   :: dmlag


!-----------------------------------------------------------------------
!  Parameters for energy-check plot edits
!-----------------------------------------------------------------------
!  n_eplt   : the unit number for e_chk data.
!
!  i_eplt   : energy-check plot switch
!
!     0    : bypass writing e_chk data.
!     ne 0 : write e_chk data.
!
!  dt_eplot : write to file e_chk.d every dt_eplot ms.
!-----------------------------------------------------------------------

INTEGER                                             :: n_eplt
INTEGER                                             :: i_eplt

REAL(KIND=double)                                   :: dt_eplot

!-----------------------------------------------------------------------
!  Parameters for pinch parameter plot edits
!-----------------------------------------------------------------------
!  r_pinch : the radius at which to evaluate the pinch parameters.
!
!  d_pinch : the density at which to evaluate the pinch parameters.
!-----------------------------------------------------------------------

REAL(KIND=double)                                   :: r_pinch
REAL(KIND=double)                                   :: d_pinch

!-----------------------------------------------------------------------
!  Parameters for 2-D u-edits
!-----------------------------------------------------------------------
!
!  ied2Du    : switch for performing 2D u-edits
!
!     0   : bypass
!     ne 0: edit when criteria satisfied
!
!  ned2Du    : edit counter
!  inted2Du  : edit every inted2Du cycles
!  n_edit2Du : model number counter
!-----------------------------------------------------------------------

INTEGER                                             :: ied2Du
INTEGER                                             :: ned2Du
INTEGER                                             :: inted2Du
INTEGER                                             :: n_edit2Du

!-----------------------------------------------------------------------
!  Parameters for 2D v (angular) edits
!-----------------------------------------------------------------------
!
!  ied2Dv    : switch for performing 2D v-edits
!
!  0    : bypass
!  ne 0 : edit when criteria satisfied
!
!  ned2Dv    : edit counter
!  inted2Dv  : edit every inted2Dv cycles
!  n_edit2Dv : model number counter
!-----------------------------------------------------------------------

INTEGER                                             :: ied2Dv
INTEGER                                             :: ned2Dv
INTEGER                                             :: inted2Dv
INTEGER                                             :: n_edit2Dv

!-----------------------------------------------------------------------
!  Parameters for 2-D w (azimuthal) edits
!-----------------------------------------------------------------------
!
!  ied2Dw    : switch for performing 2D w (azimuthal) edits
!
!  0    : bypass
!  ne 0 : edit when criteria satisfied
!
!  ned2Dw    : edit counter
!  inted2Dw  : edit every inted2Dw cycles
!  n_edit2Dw : model number counter
!-----------------------------------------------------------------------

INTEGER                                             :: ied2Dw
INTEGER                                             :: ned2Dw
INTEGER                                             :: inted2Dw
INTEGER                                             :: n_edit2Dw

!-----------------------------------------------------------------------
!  Parameters for 2D s (entropy) edits
!-----------------------------------------------------------------------
!
!  ied2Ds    : switch for performing 2D s-edits
!
!     0    : bypass
!     ne 0 : edit when criteria satisfied
!
!  ned2Ds    : edit counter
!  inted2Ds  : edit every inted2Ds cycles
!  n_edit2Ds : model number counter
!-----------------------------------------------------------------------

INTEGER                                             :: ied2Ds
INTEGER                                             :: ned2Ds
INTEGER                                             :: inted2Ds
INTEGER                                             :: n_edit2Ds

!-----------------------------------------------------------------------
!  Parameters for 2-D t (temperature) edits
!-----------------------------------------------------------------------
!
!  ied2Dt    : switch for performing 2D t (temperature) edits
!
!     0    : bypass
!     ne 0 : edit when criteria satisfied
!
!  ned2Dt    : edit counter
!  inted2Dt  : edit every inted2Dt cycles
!  n_edit2Dt : model number counter
!-----------------------------------------------------------------------

INTEGER                                             :: ied2Dt
INTEGER                                             :: ned2Dt
INTEGER                                             :: inted2Dt
INTEGER                                             :: n_edit2Dt

!-----------------------------------------------------------------------
!  Parameters for 2-D d (density) edits
!-----------------------------------------------------------------------
!
!  ied2Dd    : switch for performing 2D d (density) edits
!
!     0    : bypass
!     ne 0 : edit when criteria satisfied
!
!  ned2Dd    : edit counter
!  inted2Dd  : edit every inted2Dd cycles
!  n_edit2Dd : model number counter
!-----------------------------------------------------------------------

INTEGER                                             :: ied2Dd
INTEGER                                             :: ned2Dd
INTEGER                                             :: inted2Dd
INTEGER                                             :: n_edit2Dd

!-----------------------------------------------------------------------
!  Parameters for 2-D e (internal energy) edits
!-----------------------------------------------------------------------
!
!  ied2De    : switch for performing 2D a (internal energy) edits
!
!     0    : bypass
!     ne 0 : edit when criteria satisfied
!
!  ned2De    : edit counter
!  inted2De  : edit every inted2De cycles
!  n_edit2De : model number counter
!-----------------------------------------------------------------------

INTEGER                                             :: ied2De
INTEGER                                             :: ned2De
INTEGER                                             :: inted2De
INTEGER                                             :: n_edit2De

!-----------------------------------------------------------------------
!  Parameters for 2-D p (pressure) edits
!-----------------------------------------------------------------------
!
!  ied2Dp    : switch for performing 2D p (pressure) edits
!
!     0    : bypass
!     ne 0 : edit when criteria satisfied
!
!  ned2Dp    : edit counter
!  inted2Dp  : edit every inted2Dp cycles
!  n_edit2Dp : model number counter
!-----------------------------------------------------------------------

INTEGER                                             :: ied2Dp
INTEGER                                             :: ned2Dp
INTEGER                                             :: inted2Dp
INTEGER                                             :: n_edit2Dp

!-----------------------------------------------------------------------
!  Parameters for 2-D enu (neutrino energy density) edits
!-----------------------------------------------------------------------
!
!  ied2Denu    : switch for performing 2D enu (neutrino energy density)
!   edits
!
!     0    : bypass
!     ne 0 : edit when criteria satisfied
!
!  ned2Denu    : edit counter
!  inted2Denu  : edit every inted2De cycles
!  n_edit2Denu : model number counter
!-----------------------------------------------------------------------

INTEGER                                             :: ied2Denu
INTEGER                                             :: ned2Denu
INTEGER                                             :: inted2Denu
INTEGER                                             :: n_edit2Denu

!-----------------------------------------------------------------------
!  Parameters for 2-D fnu (neutrino flux) edits
!-----------------------------------------------------------------------
!
!  ied2Dfnu    : switch for performing 2D fnu (neutrino flux) edits
!
!     0    : bypass
!     ne 0 : edit when criteria satisfied
!
!  ned2Dfnu    : edit counter
!  inted2Dfnu  : edit every inted2Dfnu cycles
!  n_edit2Dfnu : model number counter
!-----------------------------------------------------------------------

INTEGER                                             :: ied2Dfnu
INTEGER                                             :: ned2Dfnu
INTEGER                                             :: inted2Dfnu
INTEGER                                             :: n_edit2Dfnu

!-----------------------------------------------------------------------
!  Parameters for 2-D a (mean mass number) edits
!-----------------------------------------------------------------------
!
!  ied2Da    : switch for performing 2D a (mean mass number) edits
!
!     0    : bypass
!     ne 0 : edit when criteria satisfied
!
!  ned2Da    : edit counter
!  inted2Da  : edit every inted2Da cycles
!  n_edit2Da : model number counter
!-----------------------------------------------------------------------

INTEGER                                             :: ied2Da
INTEGER                                             :: ned2Da
INTEGER                                             :: inted2Da
INTEGER                                             :: n_edit2Da

!-----------------------------------------------------------------------
!  Parameters for 2-D x (parameter) edits
!-----------------------------------------------------------------------
!
!  ied2Dx    : switch for performing 2D x (parameter) edits
!
!     0    : bypass
!     ne 0 : edit when criteria satisfied
!
!  ned2Dx    : edit counter
!  inted2Dx  : edit every inted2Dx cycles
!  n_edit2Dx : model number counter
!-----------------------------------------------------------------------

INTEGER                                             :: ied2Dx
INTEGER                                             :: ned2Dx
INTEGER                                             :: inted2Dx
INTEGER                                             :: n_edit2Dx

!-----------------------------------------------------------------------
!  Parameters for 2-D ye-edits
!-----------------------------------------------------------------------
!
!  ied2Dye   : switch for performing 2D ye-edits
!
!     0    : bypass
!     ne 0 : edit when criteria satisfied
!
!  ned2Dye    : edit counter
!  inted2Dye  : edit every inted2Dye cycles
!  n_edit2Dye : model number counter
!-----------------------------------------------------------------------

INTEGER                                             :: ied2Dye
INTEGER                                             :: ned2Dye
INTEGER                                             :: inted2Dye
INTEGER                                             :: n_edit2Dye

!-----------------------------------------------------------------------
!  Parameters for 2-D cm (composition) edits
!-----------------------------------------------------------------------
!
!  ied2Dcm   : switch for performing 2D cm (composition) edits
!
!     0    : bypass
!     ne 0 : edit when criteria satisfied
!
!  ned2Dcm    : edit counter
!  inted2Dcm  : edit every inted2Dcm cycles
!  n_edit2Dcm : model number counter
!-----------------------------------------------------------------------

INTEGER                                             :: ied2Dcm
INTEGER                                             :: ned2Dcm
INTEGER                                             :: inted2Dcm
INTEGER                                             :: n_edit2Dcm

!-----------------------------------------------------------------------
!  Parameters for 2-D nu (neutrino energy density) edits
!-----------------------------------------------------------------------
!
!  ied2Dnu   : switch for performing 2D nu (neutrino energy density) edits
!
!     0    : bypass
!     ne 0 : edit when criteria satisfied
!
!  ned2Dnu    : edit counter
!  inted2Dnu  : edit every inted2Dnu cycles
!  n_edit2Dnu : model number counter
!-----------------------------------------------------------------------

INTEGER                                             :: ied2Dnu
INTEGER                                             :: ned2Dnu
INTEGER                                             :: inted2Dnu
INTEGER                                             :: n_edit2Dnu

!-----------------------------------------------------------------------
!  Parameters for 2-D nc (nuclear energy deposition) edits
!-----------------------------------------------------------------------
!
!  ied2Dnc   : switch for performing 2D nc (neutrino energy flux) edits
!
!     0    : bypass
!     ne 0 : edit when criteria satisfied
!
!  ned2Dnc    : edit counter
!  inted2Dnc  : edit every inted2Dnc cycles
!  n_edit2Dnc : model number counter
!-----------------------------------------------------------------------

INTEGER                                             :: ied2Dnc
INTEGER                                             :: ned2Dnc
INTEGER                                             :: inted2Dnc
INTEGER                                             :: n_edit2Dnc

!-----------------------------------------------------------------------
!  Parameters for 2-D nl-edits (neutrino luminosity)
!-----------------------------------------------------------------------
!
!  ied2Dnl   : switch for performing 2D ye-edits
!
!     0    : bypass
!     ne 0 : edit when criteria satisfied
!
!  ned2Dnl    : edit counter
!  inted2Dnl  : edit every inted2Dye cycles
!  n_edit2Dnl : model number counter
!-----------------------------------------------------------------------

INTEGER                                             :: ied2Dnl
INTEGER                                             :: ned2Dnl
INTEGER                                             :: inted2Dnl
INTEGER                                             :: n_edit2Dnl

!-----------------------------------------------------------------------
!  Parameters for 2-D ne-edits (neutrino rms energy)
!-----------------------------------------------------------------------
!
!  ied2Dne   : switch for performing 2D ye-edits
!
!     0    : bypass
!     ne 0 : edit when criteria satisfied
!
!  ned2Dne    : edit counter
!  inted2Dne  : edit every inted2Dye cycles
!  n_edit2Dne : model number counter
!-----------------------------------------------------------------------

INTEGER                                             :: ied2Dne
INTEGER                                             :: ned2Dne
INTEGER                                             :: inted2Dne
INTEGER                                             :: n_edit2Dne

!-----------------------------------------------------------------------
!  Parameters for 2-D gx-edits (x-gravitational acceleration)
!-----------------------------------------------------------------------
!
!  ied2Dgx   : switch for performing 2D ye-edits
!
!     0    : bypass
!     ne 0 : edit when criteria satisfied
!
!  ned2Dgx    : edit counter
!  inted2Dgx  : edit every inted2Dye cycles
!  n_edit2Dgx : model number counter
!-----------------------------------------------------------------------

INTEGER                                             :: ied2Dgx
INTEGER                                             :: ned2Dgx
INTEGER                                             :: inted2Dgx
INTEGER                                             :: n_edit2Dgx

!-----------------------------------------------------------------------
!  Parameters for 2-D gy-edits (y-gravitational acceleration)
!-----------------------------------------------------------------------
!
!  ied2Dgy   : switch for performing 2D ye-edits
!
!     0    : bypass
!     ne 0 : edit when criteria satisfied
!
!  ned2Dgy    : edit counter
!  inted2Dgy  : edit every inted2Dye cycles
!  n_edit2Dgy : model number counter
!-----------------------------------------------------------------------

INTEGER                                             :: ied2Dgy
INTEGER                                             :: ned2Dgy
INTEGER                                             :: inted2Dgy
INTEGER                                             :: n_edit2Dgy

!-----------------------------------------------------------------------
!  Parameters for 2-D BVw-edits (Brunt-Vaisala frequency)
!-----------------------------------------------------------------------
!
!  ied2DBVw  : switch for performing 2D ye-edits
!
!     0    : bypass
!     ne 0 : edit when criteria satisfied
!
!  ned2DBVw    : edit counter
!  inted2DBVw  : edit every inted2Dye cycles
!  n_edit2DBVw : model number counter
!-----------------------------------------------------------------------

INTEGER                                             :: ied2DBVw
INTEGER                                             :: ned2DBVw
INTEGER                                             :: inted2DBVw
INTEGER                                             :: n_edit2DBVw

!-----------------------------------------------------------------------
!  Parameters for 2-D yl-edits (lepton fraction)
!-----------------------------------------------------------------------
!
!  ied2Dyl  : switch for performing 2D ye-edits
!
!     0    : bypass
!     ne 0 : edit when criteria satisfied
!
!  ned2Dyl    : edit counter
!  inted2Dyl  : edit every inted2Dyl cycles
!  n_edit2Dyl : model number counter
!-----------------------------------------------------------------------

INTEGER                                             :: ied2Dyl
INTEGER                                             :: ned2Dyl
INTEGER                                             :: inted2Dyl
INTEGER                                             :: n_edit2Dyl 

!-----------------------------------------------------------------------
!  Parameters for 2-D psi-edits (zero and first neutrino distribution
!   moments)
!-----------------------------------------------------------------------
!
!  ied2Dpsi  : switch for performing 2D psi-edits
!
!     0    : bypass
!     ne 0 : edit when criteria satisfied
!
!  ned2Dpsi    : edit counter
!  inted2Dpsi  : edit every inted2Dpsi cycles
!  n_edit2Dpsi : model number counter
!-----------------------------------------------------------------------

INTEGER                                             :: ied2Dpsi
INTEGER                                             :: ned2Dpsi
INTEGER                                             :: inted2Dpsi
INTEGER                                             :: n_edit2Dpsi 

!-----------------------------------------------------------------------
!  Parameters for 2-D xn-edits (composition mass fractions)
!-----------------------------------------------------------------------
!
!  ied2Dxn  : switch for performing 2D xn-edits
!
!     0    : bypass
!     ne 0 : edit when criteria satisfied
!
!  ned2Dxn    : edit counter
!  inted2Dxn  : edit every inted2Dxn cycles
!  n_edit2Dxn : model number counter
!-----------------------------------------------------------------------

INTEGER                                             :: ied2Dxn
INTEGER                                             :: ned2Dxn
INTEGER                                             :: inted2Dxn
INTEGER                                             :: n_edit2Dxn 

!-----------------------------------------------------------------------
!  Parameters for 2-D aux-edits (auxiliaru nucleus parameters)
!-----------------------------------------------------------------------
!
!  ied2Daux : switch for performing 2D aux-edits
!
!     0    : bypass
!     ne 0 : edit when criteria satisfied
!
!  ned2Daux    : edit counter
!  inted2Daux  : edit every inted2Daux cycles
!  n_edit2Daux : model number counter
!-----------------------------------------------------------------------

INTEGER                                             :: ied2Daux
INTEGER                                             :: ned2Daux
INTEGER                                             :: inted2Daux
INTEGER                                             :: n_edit2Daux 

!-----------------------------------------------------------------------
!  Parameters for 2-D nse-edits (Nuclear Statistical Equilibrium)
!-----------------------------------------------------------------------
!
!  ied2Dnse  : switch for performing 2D nse-edits
!
!     0    : bypass
!     ne 0 : edit when criteria satisfied
!
!  ned2Dnse    : edit counter
!  inted2Dnse  : edit every inted2Dnse cycles
!  n_edit2Dnse : model number counter
!-----------------------------------------------------------------------

INTEGER                                             :: ied2Dnse
INTEGER                                             :: ned2Dnse
INTEGER                                             :: inted2Dnse
INTEGER                                             :: n_edit2Dnse 

!-----------------------------------------------------------------------
!  Time parameters for 2-D edits
!-----------------------------------------------------------------------
!
!  i_edit2D   : switch for performing 2D edits at selected times
!  n_edit2D   : model number counter
!  dt_2Dedit1 : dump files for 2D edits every dt_2Dedit1 ms before bounce.
!  dt_2Dedit2 : dump files for 2D edits every dt_2Dedit2 ms after bounce.
!-----------------------------------------------------------------------

INTEGER                                             :: i_edit2D
INTEGER                                             :: n_edit2D

REAL(KIND=double)                                   :: dt_2Dedit1
REAL(KIND=double)                                   :: dt_2Dedit2

!-----------------------------------------------------------------------
!  Cycle parameters for global edits
!-----------------------------------------------------------------------
!
!  ied_global_n  : switch for performing global using cycle criterion
!
!     0    : bypass
!     ne 0 : edit when criteria satisfied
!
!  nned_global   : edit counter
!  inted_global  : edit every inted_global cycles
!  n_edit_global : model number counter
!-----------------------------------------------------------------------

INTEGER                                             :: ied_global_n
INTEGER                                             :: nned_global
INTEGER                                             :: inted_global
INTEGER                                             :: n_edit_global

!-----------------------------------------------------------------------
!  Time parameters for global edits
!-----------------------------------------------------------------------
!
!  ied_global_t  : switch for performing 2D edits at selected times
!
!     0    : bypass
!     ne 0 : edit when criteria satisfied
!
!  nted_global   : model number counter
!  dt_global_ed1 : dump files for global edits every dt_global_ed1 ms
!                   before bounce.
!  dt_global_ed2 : dump files for global edits every dt_global_ed2 ms
!                   after bounce.
!-----------------------------------------------------------------------

INTEGER                                             :: ied_global_t
INTEGER                                             :: nted_global

REAL(KIND=double)                                   :: dt_global_ed1
REAL(KIND=double)                                   :: dt_global_ed2

!-----------------------------------------------------------------------
!  Parameters for HDF 2-D edits
!-----------------------------------------------------------------------
!
!  nd_MDedit  : the unit number for dumping files for MD edits.
!
!  i_MDedit   : switch for dumping files for MD edits.
!
!  n_MDedit   : MD edit counter.
!
!     0    : bypass dumping files for MD edits.
!     ne 0 : dump files for MD edits.
!
!  dt_MDedit1 : dump files for MD edits every dt_MDedit1 ms before bounce.
!  dt_MDedit2 : dump files for MD edits every dt_MDedit1 ms after bounce.
!-----------------------------------------------------------------------

INTEGER                                             :: i_MDedit
INTEGER                                             :: nd_MDedit
INTEGER                                             :: n_MDedit

REAL(KIND=double)                                   :: dt_MDedit1
REAL(KIND=double)                                   :: dt_MDedit2

!-----------------------------------------------------------------------
!  Parameters for temperary restart dumps
!-----------------------------------------------------------------------
!  intrst : the number of cycles between temporary restart dumps.
!
!  nnrst  : the number of cycles since the last temporary restart dump.
!
!  nrstd1 : the unit number for file 'rstdmp1.d' to which temporary
!   restart dumps are written every intrst cycles. These restart dumps
!   are temporary in the sense that the next restart dump written to
!   file 'rstdmp1.d' writes over the preceding dump.
!
!  nrstd2 : the unit number for file 'rstdmp2.d' to which temporary
!   restart dumps can be written every other intrst cycles alternating
!   with file rstdmp1.d.
!
!  nrstfl : a temporary restart dump parameter -
!
!     2    : restart dumps are alternated between file 'rstdmp1.d' and
!      file 'rstdmp2.d'.                                       
!     ne 2 : restart dumps are written only to file 'rstdmp1.d'.
!
!  input_flag : Input flag directing the input files to be read on
!   restart.
!
!   If cycle_flag > 0 and ndim = 1 (ny = 1 and nz = 1) problem are read
!    by processor 0 from file
!
!     rst_tmp1_keys.d           : input_flag = 1
!     rst_tmp2_keys.d           : input_flag = 2
!     restart_CYCLE_FLAG.d      : input_flag = 3
!     restart_final.d           : input_flag = 4
!     restart_CYCLE_FLAG.d      : input_flag = 5
!
!   The nuclear abundance data and model configuration are read from
!    file
!
!     rst_tmp1_model.d          : input_flag = 1
!     rst_tmp2_model.d          : input_flag = 2
!     restart_modelCYCLE_FLAG.d ! input_flag = 3
!     restart_final_mod.d       : input_flag = 4
!     restart_modelCYCLE_FLAG.d : input_flag = 5
!
!   If cycle_flag > 0 and ndim > 1 (ny > 1) problem keys are read by
!    processor 0 and broadcast to all processors from file
!
!     rst_tmp1_keys.d           : input_flag = 1
!     rst_tmp2_keys.d           : input_flag = 2
!     restart_CYCLE_FLAG.d      : input_flag = 3
!     restart_final.d           : input_flag = 4
!     restart_CYCLE_FLAG.d      : input_flag = 5
!
!   The nuclear abundance data and model configuration are read by all
!    processors from file
!
!     rst_tmp1_MYID.d                       : input_flag = 1
!     rst_tmp2_MYID.d                       : input_flag = 2
!     restart_modelMYID_CYCLE_FLAG.d        : input_flag = 3
!     restart_final_modMYID.d               : input_flag = 4
!-----------------------------------------------------------------------

INTEGER                                             :: intrst
INTEGER                                             :: nnrst
INTEGER                                             :: nrstd1
INTEGER                                             :: nrstd2
INTEGER                                             :: nrstfl
INTEGER                                             :: input_flag

!-----------------------------------------------------------------------
!  Parameters for permanent restart dumps
!-----------------------------------------------------------------------
!  noutpmt   : the unit number of the permanent restart dump file
!   (typicallyfile 'restartxxxxxxx.d', where xxxxxxx is the cycle number
!   right-justified).
!
!  intprm    : the number of cycles between permanent restart dumps.
!
!  nprm      : the number of cycles since the last permanent restart dump.
!
!  ncyrst(i) : the cycle number at which the ith prescribed permanent
!   restart dump is implemented.
!
!  ncychg    : at termination, write to restart file with
!     ncymax = ncymax + ncychg.
!
!  irstbgn   : restart dump parameter -
!
!     1   : temporary and permanent restart dumps are implemented at
!            the beginning of the run (ncycle = 0)
!     ne 1: temporary and permanent restart dumps are bypassed at the
!            beginning of the run (ncycle = 0)
!
!  noutfl    :
!-----------------------------------------------------------------------

INTEGER                                             :: noutpmt
INTEGER                                             :: intprm
INTEGER                                             :: nprm
INTEGER, DIMENSION(100)                             :: ncyrst
INTEGER                                             :: ncychg
INTEGER                                             :: irstbgn
INTEGER                                             :: noutfl

!-----------------------------------------------------------------------
!  Work
!-----------------------------------------------------------------------
!  pdv(j)   : work done per unit mass by zone j-1/2 in compression or
!   expansion (i.e., at the expense of its internal energy) since the
!   beginning of the calculation.
!
!  twrk(j): total work done per unit mass by zone j-1/2 minus pdv(j)
!   (i.e., at the expense of its kinetic energy) since the beginning
!   of the calculation.
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)        :: pdv
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)        :: twrk

!-----------------------------------------------------------------------
!  Stresses
!-----------------------------------------------------------------------
!  nustrss(j) : force per unit mass at the boundary between radial zone
!   j and j+1 due to neutrinos of all types.
!
!  pstrss(j)  : force per unit mass at the boundary between radial zone
!   j and j+1 due to material pressure.
!
!  gstrss(j)  : force per unit mass at the boundary between radial  zone
!   j and j+1 due to gravity.
!
!  gstrss_cx(j) : x-component of force per unit mass at the center of
!   radial zone j due to gravity.
!
!  gstrss_cy(j) : x-component of force per unit mass at the center of
!   radial zone j due to gravity.
!
!  rstrss(j)  : force per unit mass at the boundary between radial zone
!   j and j+1 due to  relativistic terms.
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)        :: nustrss
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)        :: pstrss
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)        :: gstrss
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)        :: gstrss_cx
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)        :: gstrss_cy
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)        :: rstrss

!-----------------------------------------------------------------------
!  Gravotational potential
!-----------------------------------------------------------------------
!  g_pot(j) : gravitational potential of radal zone j
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)        :: g_pot

!-----------------------------------------------------------------------
!  Neutron star parameters
!-----------------------------------------------------------------------
!  mass_ns : mass of the neutron star (solar masses)
!
!  vel_ns  : velocityh of the neutron star (km s^{-1})
!-----------------------------------------------------------------------

REAL(KIND=double)                                   :: mass_ns
REAL(KIND=double)                                   :: vel_ns

!-----------------------------------------------------------------------
!  Solid angles subtended by rays
!-----------------------------------------------------------------------
!  d_omega : solid angles subtended by radial rays
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)        :: d_omega

!-----------------------------------------------------------------------
!  Angularly averaged quantities
!-----------------------------------------------------------------------
!  rhobar : zone shifted angular averaged density [g cm^{-3}]
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)        :: rhobar

!-----------------------------------------------------------------------
!  Neutrino-matter energy transfer rates
!-----------------------------------------------------------------------
!  dudt_ABEM(j,n) : energy transferred to matter due to the emission
!   and absorption of n-neutrinos [MeV/nucleon]
!
!  dudt_NES(j,n)  : energy transferred to matter due to the neutrino-
!   electron scattering of n-neutrinos [MeV/nucleon]
!
!  dudt_NPS(j,n)  : energy transferred to matter due to the neutrino-
!   positron scattering of n-neutrinos [MeV/nucleon]
!
!  dudt_NNS(j,n)  : energy transferred to matter due to the neutrino-
!   nucleon elastic scattering of n-neutrinos [MeV/nucleon]
!
!  dudt_NNNS(j,n) : energy transferred to matter due to the neutrino-
!   nucleon inelastic scattering of n-neutrinos [MeV/nucleon]
!
!  dudt_NAS(j,n)  : energy transferred to matter due to the neutrino-
!   nucleus inelastic scattering of n-neutrinos [MeV/nucleon]
!
!  dudt_PR(j,n)   : energy transferred to matter due to the electron-
!   positron pair annihilation into n-neutrinos-antineutrino pairs
!   [MeV/nucleon]
!
!  dudt_Brem(j,n) : energy transferred to matter due to the nucleon-
!   nucleon pair production into n-neutrinos-antineutrino pairs
!   [MeV/nucleon]
!
!  dudt_NET(j) : net energy transferred to matter due to the neutrino-
!   matter interactions [MeV/nucleon]
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)      :: dudt_ABEM
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)      :: dudt_NES
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)      :: dudt_NPS
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)      :: dudt_NNS
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)      :: dudt_NAS
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)      :: dudt_PR
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)      :: dudt_Brem
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)        :: dudt_NET

!-----------------------------------------------------------------------
! File per processor output control
!-----------------------------------------------------------------------
!  io_per_processor_mode : determines which processes write individual
!                          I/O files (models_n, intl_model, etc.)
!                   Values:  0: write on root processor only (mype=0,default)
!                            1: write on all processes
!   write_local          : logical set to true if above condition is met
!                             locally for each processor
!
!-----------------------------------------------------------------------

INTEGER :: io_per_processor_mode = 0
LOGICAL :: write_local           = .true.

!-----------------------------------------------------------------------
! Per cycle output control
!-----------------------------------------------------------------------
! write_cycle_d      : (logical) write cycle number to file "cycle.d"
! no_cycle_d         : (key) set to 1 will stop write of cycle.d
! superdump_interval : cycles between flushes of "superdump.d"
!-----------------------------------------------------------------------

INTEGER :: superdump_interval = 100
INTEGER :: no_cycle_d         = 1
LOGICAL :: write_cycle_d      = .false.

END module edit_module
