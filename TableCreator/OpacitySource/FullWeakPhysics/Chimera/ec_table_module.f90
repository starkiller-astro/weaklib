MODULE ec_table_module

  USE wlKindModule, ONLY: dp

  USE, INTRINSIC :: iso_fortran_env, only : stdin=>input_unit, &
                                            stdout=>output_unit, &
                                            stderr=>error_unit

  implicit none

  INTEGER, PARAMETER  :: npts=200, ntemp=37, nrho=11, nye=8 ! Revised Table
  !INTEGER, PARAMETER  :: npts=200, ntemp=35, nrho=11, nye=7 ! Original Table
  INTEGER, PARAMETER  :: tmx=40, tmn=(tmx+1-ntemp), rhomx=26, rhomn=(rhomx+1-nrho)

  REAL(dp), PARAMETER :: deltaE=0.5d0,deltalrho_ec=0.5d0,deltaT_ec=0.1d0,deltaYe_ec=0.02d0

  REAL(dp)            :: energy(0:npts)
  REAL(dp)            :: spectab(0:npts,rhomn:rhomx,tmn:tmx,0:(nye-1))
  REAL(dp)            :: ratetab(rhomn:rhomx,tmn:tmx,0:(nye-1))

  INTEGER             :: yefloor(rhomn:rhomx)

  CONTAINS

    SUBROUTINE read_ec_table( ec_table )
    !-----------------------------------------------------------------------
    !
    !    Author:       W. R. Hix, ORNL
    !                  E. J. Lentz, UTK
    !
    !    Date:         7/19/07
    !
    !    Purpose:
    !      To read the tabular electron capture grid and data.
    !
    !    Subprograms called:
    !      none
    !
    !    Input arguments:
    !      none
    !
    !    Output arguments:
    !      none
    !
    !    Input arguments (module):
    !
    !  iaenct              : 0, e-neutrino emission from nuclei
    !                          (tabular rates) omitted
    !                      : 1, e-neutrino emission from nuclei
    !                          (tabular rates) included
    !  ntemp, nrho,nye     : number of grid points in respective directions.
    !  tmn,tmx             : minimum and maximum temperature grid points
    !  rhomn,rhomx         : minimum and maximum density grid points
    !
    !    Output arguments (module):
    !  spectab              : nuclear electron capture emissivity table
    !  energy              : table energy grid
    !
    !-----------------------------------------------------------------------

    USE numerical_module, ONLY: pi
    USE physcnst_module, ONLY: cvel, hbar

    !USE ec_table_module
    !USE edit_module, ONLY: nlog, nprint
    USE prb_cntl_module, ONLY: iaenct

    !USE parallel_module
    !USE mpi

    IMPLICIT NONE

    character(len=*),intent(in)   :: ec_table

    REAL(dp)                       :: spec(0:npts)
    INTEGER, parameter             :: nrdata = 31   ! unit number for accessing the data files
    INTEGER                        :: istat         ! open-close file flag
    INTEGER                        :: i_extent      ! length of broadcast data
    INTEGER                        :: j,k
    integer :: itemp,irho,iye,it,iy,ir
    character(len(ec_table)) :: tablefile
    character (LEN=120) :: desc(2)
    integer :: izone,nline,npart
    real(dp) :: rho,ye,temp,rate,dum
    logical :: ispart
    integer, parameter :: nperline=7, outperline=6
    real(dp), parameter :: gigK=11.60451 !MeV
    real(dp), parameter :: eps=0.01 ! 1% tolerance on grid readin
    real(dp), parameter :: const= 2.0*(pi*cvel)**2*hbar**3

    !-----------------------------------------------------------------------
    ! Return if tabular ec rates are not used.
    !-----------------------------------------------------------------------
    IF (iaenct /=1 ) THEN
      RETURN
    ENDIF

    !-----------------------------------------------------------------------
    !  Determine spectrum input file geometry
    !-----------------------------------------------------------------------

    nline=(npts+1)/nperline
    npart=(npts+1)-nline*nperline
    ispart= npart.gt.0

    !-----------------------------------------------------------------------
    !  Initialize table energy grid data
    !-----------------------------------------------------------------------
    DO k=0,npts
      energy(k)=dble(k)*deltaE
    ENDDO

    !IF ( myid == 0 ) THEN

    spectab=0.d0

    !-----------------------------------------------------------------------
    !  Open file
    !-----------------------------------------------------------------------

    !tablefile = "ec_table.d"
    tablefile = ec_table

    OPEN (nrdata,file=trim(tablefile),status='old',iostat=istat)

    IF ( istat /= 0 ) THEN
      WRITE (*,*) 'EC_table data file ',trim(tablefile),' not found.'
      !CALL MPI_ABORT( MPI_COMM_WORLD, 333, ierr )
    ENDIF

    !-----------------------------------------------------------------------
    !  Read table data
    !-----------------------------------------------------------------------
    READ(nrdata,"(a)") desc(1)
    READ(nrdata,"(a)") desc(2)

    DO it=tmn,tmx
     DO ir=rhomn,rhomx
      DO iy=(nye-1),0,-1
        READ(nrdata,"(6x,i2,6x,f12.6,5x,f9.6,5x,f12.6,5x,f9.6,5x,f12.6)") &
    &     izone,rho,ye,temp,dum,rate
        irho=int(log10(rho)/deltalrho_ec+eps)
        itemp=int(temp/gigK/deltaT_ec + eps)
        iye=int(ye/deltaYe_ec + eps)
        IF (irho /= ir) THEN
          WRITE(*,*) 'EC table misread: density ',irho, ' /= ',ir
          stop
        Elseif (itemp /= it) THEN
          WRITE(*,*) 'EC table misread: temperature ',itemp, ' /= ',it
          stop
        ENDIF
        READ(nrdata,*)
        DO j=1,nline
         READ(nrdata,*) (spec(k),k=(j-1)*nperline,j*nperline-1)
    !    WRITE(nprint,'(7es13.6)') (spec(k),k=(j-1)*nperline,j*nperline-1)
        ENDDO
        IF (ispart) READ(nrdata,*) (spec(k),k=nline*nperline,npts)
    !   IF (ispart) WRITE(6,'(7es13.6)') (spec(k),k=nline,nperline,npts)
        !spectab(:,ir,it,iy)=rate*spec
        spectab(:,ir,it,iy) = spec
        ratetab( ir,it,iy) = rate
      ENDDO
      IF ( it == tmx ) THEN
       yefloor(ir)=iye
    !   WRITE(nprint,'(a,es10.3,a,f5.3)') "yefloor for rho=",rho,': ',.02*yefloor(ir)
      ENDIF
     ENDDO
    ENDDO

    CLOSE(nrdata)

    !-----------------------------------------------------------------------
    !  Multiply spectrum by rate and take log of table
    !-----------------------------------------------------------------------

    WHERE (spectab <= 0.d0)
      spectab = -200.d0
    ELSEWHERE
      !spectab = log10(const*spectab)
      spectab = log10(spectab)
    ENDWHERE

    ratetab = log10(ratetab)

    !-----------------------------------------------------------------------
    !  Signal completion
    !-----------------------------------------------------------------------

    WRITE(*,*) ' Finished reading ',trim(desc(1))
    WRITE(*,*) ' Table grid:',trim(desc(2))
    !ENDIF ! myid == 0

    RETURN

    END SUBROUTINE read_ec_table

    !SUBROUTINE interp_ec( spec, rate, unui, unubi, dunui, temp, rho, Ye, nez )
    SUBROUTINE interp_ec( spec, rate, unui, temp, rho, Ye, nez )
    !-----------------------------------------------------------------------
    !
    !    Author:       W.R. Hix, Physics Division
    !                  Oak Ridge National Laboratory, Oak Ridge TN 37831
    !                  E.J. Lentz, Department of Physics and Astronomy
    !                  University of Tennessee, Knoxville TN 37919
    !
    !    Date:         7/17/07
    !
    !    Purpose:
    !      To compute the inverse mean free paths for the absorption and emission
    !       of neutrinos using the NSE folded table of electron capture
    !       rates.
    !
    !    Variables that must be passed through common:
    !        spec, energy, spectab
    !
    !    Subprograms called:
    !        none
    !
    !    Input arguments:
    !  rho         : matter density [g cm^{-3}]
    !  temp        : matter temperature [K]
    !  ye          : electron fraction
    !
    !    Output arguments:
    !  spec        :  emission spectrum inverse mean free path (/cm)
    !  rate        :  electron capture rate 
    !
    !    Input arguments (modules):
    !
    !  iaenct     : 0; inverse mean free paths set to zero.
    !               1; inverse mean free paths computed.
    !  roaenct    : density above which rates are set to zero.
    !
    !-----------------------------------------------------------------------

    USE physcnst_module, ONLY: kmev, cvel, hbar
    USE numerical_module, ONLY: pi

    IMPLICIT NONE

    !-----------------------------------------------------------------------
    !  Output variables
    !-----------------------------------------------------------------------

    INTEGER, intent(in) :: nez
    REAL(dp), dimension(nez), intent(out) :: spec
    REAL(dp), intent(out)                 :: rate

    !-----------------------------------------------------------------------
    !  Local variables
    !-----------------------------------------------------------------------

    REAL(dp), dimension(nez), intent(in)   :: unui !, dunui
    !REAL(dp), dimension(nez+1), intent(in) :: unubi

    REAL(dp) :: temp, rho, ye
    REAL(dp), DIMENSION(0:npts) :: jecfine
    REAL(dp) :: rate_interp
    REAL(dp) :: deltaupper(nez),deltalower(nez)
    REAL(dp) :: c1,c2,c3,c3base,loctot,tmev

    INTEGER :: kfmin(nez),kfmax(nez)
    INTEGER :: k, ii
    INTEGER :: rhoup,rhodown,tup,tdown
    INTEGER :: irho,iye,itemp
    INTEGER :: kmax,ktop,yeoffset(2)

    LOGICAL :: askew

    REAL(dp) :: loc_integral(nez)

    real(dp), parameter :: const= 2.0*(pi*cvel)**2*hbar**3

    !REAL(dp), DIMENSION(nez) :: unubi, dunui

    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------

    kmax = nez
    spec  = 0.d0
    rate  = 0.d0

    !-----------------------------------------------------------------------
    !  Find temperature in table.  Temperature tabulated in MeV.
    !  For temperature outside of table, set rates to small value
    !-----------------------------------------------------------------------

    tmev  = temp*kmev
    itemp = int(tmev/deltaT_ec)
    tdown = itemp
    tup   = itemp+1

    IF (tup > tmx .or. tdown < tmn ) then
      RETURN
    endif

    !-----------------------------------------------------------------------
    !  Find density in table.
    !  For density outside of table, set rates to small value
    !-----------------------------------------------------------------------

    irho=int(log10(rho)/deltalrho_ec)
    rhodown=irho
    rhoup = irho+1

    !IF (rhoup > rhomx .or. rhodown < rhomn) then
    !  RETURN
    !endif

    !Allow for extrapolation at lower density end of the table
    IF (rhoup > rhomx) RETURN

    IF (rhodown < rhomn) then
      rhodown = rhomn
      rhoup   = rhomn+1
    endif

    !-----------------------------------------------------------------------
    !  Find ye in table.  Ye grid is density dependent
    !  For ye outside of table, allow extrapolation
    !-----------------------------------------------------------------------

    iye         = INT(Ye/deltaYe_ec)
    yeoffset(1) = iye-yefloor(rhodown)
    yeoffset(2) = iye-yefloor(rhoup)

    !Comment below to allow for extrapolation in Ye
    !IF (sum(yeoffset) < -2) RETURN

    !Comment below to allow for extrapolation in Ye
    !IF (sum(yeoffset) > (2*nye-2)) RETURN

    !-----------------------------------------------------------------------
    !  Because of density dependence, Ye may require askew interpolation.
    !-----------------------------------------------------------------------

    askew=.false.
    DO ii=1,2
      IF (yeoffset(ii).lt.0) THEN
        askew        = .true.
        yeoffset(ii) = 0
      ENDIF
    ENDDO

    DO ii=1,2
      IF (yeoffset(ii).gt.(nye-2)) THEN
        askew        = .true.
        yeoffset(ii) = nye-2
      END IF
    ENDDO

    !-----------------------------------------------------------------------
    !  Compute interpolation forefactors.
    !-----------------------------------------------------------------------
    c1=(log10(rho)-dble(rhodown)*deltalrho_ec)/deltalrho_ec

    c2=(log10(tmev/(dble(itemp)*deltaT_ec)))/(log10(dble(tup)/dble(tdown)))

    c3base=(ye-dble(iye)*deltaYe_ec)/deltaYe_ec
    c3=c3base + &
    & c1*dble(yefloor(rhoup)+yeoffset(2)-yefloor(rhodown)-yeoffset(1)) !correction for skew If present

    !-----------------------------------------------------------------------
    ! Interpolate rate on table energy grid
    !-----------------------------------------------------------------------

    jecfine(:) = c3 & !note on indices: yeoffset=1 for rhodown, 2 for rhoup
    &           *( (1.0-c1)*(1.0-c2) * spectab(:,rhodown,tdown,yeoffset(1)+1) &  !   add one in upper block for old yeup
    &           +  c1*(1.0-c2)       * spectab(:,rhoup,tdown,yeoffset(2)+1) &
    &           +  c2*(1.0-c1)       * spectab(:,rhodown,tup,yeoffset(1)+1) &
    &           +  c1*c2             * spectab(:,rhoup,tup,yeoffset(2)+1) ) &
    &           + (1.0-c3) &
    &           *( (1.0-c1)*(1.0-c2) * spectab(:,rhodown,tdown,yeoffset(1)) &
    &           +  c1*(1.0-c2)       * spectab(:,rhoup,tdown,yeoffset(2)) &
    &           +  c2*(1.0-c1)       * spectab(:,rhodown,tup,yeoffset(1)) &
    &           +  c1*c2             * spectab(:,rhoup,tup,yeoffset(2)) )

    jecfine=10.0d0**(jecfine)

    rate_interp = c3 & !note on indices: yeoffset=1 for rhodown, 2 for rhoup
    &           *( (1.0-c1)*(1.0-c2) * ratetab(rhodown,tdown,yeoffset(1)+1) &  !   add one in upper block for old yeup
    &           +  c1*(1.0-c2)       * ratetab(rhoup,tdown,yeoffset(2)+1) &
    &           +  c2*(1.0-c1)       * ratetab(rhodown,tup,yeoffset(1)+1) &
    &           +  c1*c2             * ratetab(rhoup,tup,yeoffset(2)+1) ) &
    &           + (1.0-c3) &
    &           *( (1.0-c1)*(1.0-c2) * ratetab(rhodown,tdown,yeoffset(1)) &
    &           +  c1*(1.0-c2)       * ratetab(rhoup,tdown,yeoffset(2)) &
    &           +  c2*(1.0-c1)       * ratetab(rhodown,tup,yeoffset(1)) &
    &           +  c1*c2             * ratetab(rhoup,tup,yeoffset(2)) )

    rate_interp=10.0d0**(rate_interp)
!write(*,*) 'rate_interp', rate_interp

    !renormalize spectrum to 1 after interpolation

    loctot = 0.0d0
    do k=0,npts
      loctot = loctot + jecfine(k) * deltaE
    enddo

    !if(abs(loctot-1.0d0) > 0.05d0) then
    !  write(*,*) 'rho,T,Ye ', rho, temp, ye
    !  write(*,*) 'rate ', rate_interp
    !  write(*,*) 'integral ', loctot
    !endif

    !jecfine(:) = jecfine(:) * const * rate_interp / loctot
    !jecfine(:) = jecfine(:) * const / loctot
    jecfine(:) = jecfine(:) / loctot

    do k=0,npts
    spec(k+1) = jecfine(k)
    enddo
    rate = rate_interp

    END SUBROUTINE interp_ec

    SUBROUTINE abem_nuclei_EC_table_weaklib( n, unui, &
               rho, t, ye, xh, ah, cmpn, cmpp, cmpe,  &
               EC_table_spec, EC_table_rate, nse, nez )
    !-----------------------------------------------------------------------
    !
    !    Author:       W.R. Hix, Physics Division
    !                  Oak Ridge National Laboratory, Oak Ridge TN 37831
    !
    !    Date:         7/17/07
    !
    !    Purpose:
    !      To compute the inverse mean free paths for the absorption and emission
    !       of n-type neutrinos using the NSE folded table of electron capture
    !       rates.
    !
    !    Subprograms called:
    !  interp_ec_table
    !
    !    Input arguments:
    !  n           : neutrino type (1, e-neutrino; 2, e-antineutrino; 3, t-neutrino)
    !  rho         : matter density [g cm^{-3}]
    !  t           : matter temperature [K]
    !  ye          : electron fraction
    !  xh          : heavy nucleus mass fraction
    !  ah          : heavy nucleus mass number
    !  cmpn        : free neutron chemical potential (excluding rest mass) [MeV]
    !  cmpp        : free proton chemical potential (excluding rest mass) [MeV]
    !  cmpe        : electron chemical potential (including rest mass) [MeV]
    !
    !    Output arguments:
    !  EC_table_spec : electron neutrino emission spectrum normalised to 1
    !  EC_table_rate : electron capture rate at given rho,T,Ye  
    !
    !    Input arguments (module):
    !
    !  iaenct     : 0; inverse mean free paths set to zero.
    !               1; inverse mean free paths computed.
    !  roaenct    : density above which rates are set to zero.
    !
    !-----------------------------------------------------------------------

    !USE kind_module
    !USE array_module, ONLY: nez, nezext
    USE numerical_module, ONLY: zero, one
    USE physcnst_module, ONLY: rmu, kmev, dmnp

    !USE math_functions_module, ONLY: fexp
    !USE nu_energy_grid_module, ONLY: unui
    USE prb_cntl_module, ONLY: iaenct, roaenct

    IMPLICIT none

    !-----------------------------------------------------------------------
    !        Input variables.
    !-----------------------------------------------------------------------

    INTEGER, INTENT(in)     :: n               ! neutrino flavor index
    INTEGER, INTENT(in)     :: nse             ! NSE flag
    INTEGER, INTENT(in)     :: nez             !number of neutrino energy bins
    REAL(dp), dimension(nez), intent(in)   :: unui !neutrino energies
    !REAL(dp), dimension(nez+1), intent(in) :: unubi !neutrino energies, bin edges
    !REAL(dp), dimension(nez), intent(in)   :: dunui !neutrino energies, dE for edges

    REAL(dp), INTENT(in)    :: rho             ! density (g/cm^3)
    REAL(dp), INTENT(in)    :: t               ! temperature [K]
    REAL(dp), INTENT(in)    :: ye              ! electron fraction
    REAL(dp), INTENT(in)    :: xh              ! heavy nuclei mass fraction
    REAL(dp), INTENT(in)    :: ah              ! heavy nuclei mass number
    REAL(dp), INTENT(in)    :: cmpn            ! neutron chemical porential
    REAL(dp), INTENT(in)    :: cmpp            ! proton chemical porential
    REAL(dp), INTENT(in)    :: cmpe            ! electron chemical porential

    !-----------------------------------------------------------------------
    !        Output variables.
    !-----------------------------------------------------------------------

    REAL(dp), INTENT(out) :: EC_table_spec(nez)
    REAL(dp), INTENT(out) :: EC_table_rate

    !-----------------------------------------------------------------------
    !        Local variables
    !-----------------------------------------------------------------------

    INTEGER                 :: je              ! do index

    REAL(dp)                :: tmev            ! temperature [MeV]
    REAL(dp)                :: xnuc            ! number density of nuclei (cm^{-3})

    REAL(dp), EXTERNAL      :: fexp

    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    !  Set emitnc and absrnc to zero and return if
    !               iaenct = 0
    !  or
    !               rho > roaenct
    !  or
    !               n /= 1
    !-----------------------------------------------------------------------
    !IF ( iaenct == 0  .or.  rho  >  roaenct  .or.  n /= 1  .or.  ah < 40 .or. nse == 0 ) THEN
    !IF ( rho  >  roaenct  .or.  ah < 40 .or. nse == 0 ) THEN
    !IF ( rho  >  roaenct  .or.  ah < 40.0d0 ) THEN
    !  EC_table_spec = zero
    !  EC_table_rate = zero
    !  RETURN
    !END IF

    !-----------------------------------------------------------------------
    !  Initialize
    !-----------------------------------------------------------------------

    !tmev              = kmev * t
    !if(ah == 0.0d0 ) then
    !  xnuc = 0.0d0
    !else
    !  xnuc              = ( xh/( rmu * ah ) ) * rho
    !endif

    !-----------------------------------------------------------------------
    !
    !       \\\\\ E-NEUTRINO-NUCLEUS  EMISSION  AND ABSORPTION/////
    !
    !-----------------------------------------------------------------------
    !  Calculate values for rate/per heavy nucleus by interpolation
    !-----------------------------------------------------------------------

    CALL interp_ec( EC_table_spec, EC_table_rate, unui, t, rho, ye, nez )

    !-----------------------------------------------------------------------
    !  Multiply by heavy nucleus abundance
    !-----------------------------------------------------------------------

    !emitnc = emitnc * xnuc
    !EC_table_spec = EC_table_spec * xnuc
    !-----------------------------------------------------------------------
    !  Compute the absorption rate using detailed balance
    !-----------------------------------------------------------------------
    !DO je = 1,nez
    !  absrnc(je) = emitnc(je) * fexp((unui(je) + dmnp + cmpn - cmpp - cmpe)/tmev)
    !END DO

    END SUBROUTINE abem_nuclei_EC_table_weaklib

END MODULE ec_table_module
