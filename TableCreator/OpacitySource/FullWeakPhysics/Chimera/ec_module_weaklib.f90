MODULE ec_table_module

  USE wlKindModule, ONLY: dp

  USE, INTRINSIC :: iso_fortran_env, only : stdin=>input_unit, &
                                            stdout=>output_unit, &
                                            stderr=>error_unit

  implicit none

  INTEGER, PARAMETER  :: npts=200, ntemp=37, nrho=11, nye=8 ! Revised Table
  !INTEGER, PARAMETER  :: npts=200, ntemp=35, nrho=11, nye=7 ! Original Table
  INTEGER, PARAMETER  :: tmx=40, tmn=(tmx+1-ntemp), rhomx=26, rhomn=(rhomx+1-nrho)

  REAL(dp), PARAMETER :: deltaE=0.5,deltalrho_ec=.5,deltaT_ec=.1,deltaYe_ec=.02

  REAL(dp)            :: energy(0:npts)
  REAL(dp)            :: jectab(0:npts,rhomn:rhomx,tmn:tmx,0:(nye-1))

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
    !  jectab              : nuclear electron capture emissivity table
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
    character(80) :: tablefile
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

    jectab=0.d0

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
        jectab(:,ir,it,iy)=rate*spec
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

    WHERE (jectab <= 0.d0)
      jectab = -200.d0
    ELSEWHERE
      jectab = log10(const*jectab)
    ENDWHERE

    !-----------------------------------------------------------------------
    !  Signal completion
    !-----------------------------------------------------------------------

    WRITE(*,*) ' Finished reading ',trim(desc(1))
    WRITE(*,*) ' Table grid:',trim(desc(2))

    !ENDIF ! myid == 0

    RETURN

    END SUBROUTINE read_ec_table

    SUBROUTINE interp_ec( jec, e_in, unubi, dunui, temp, rho, Ye, nez )
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
    !        spec, energy, jectab
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
    !  jec         :  emission inverse mean free path (/cm)
    !
    !    Input arguments (modules):
    !
    !  iaenct     : 0; inverse mean free paths set to zero.
    !               1; inverse mean free paths computed.
    !  roaenct    : density above which rates are set to zero.
    !
    !-----------------------------------------------------------------------

    !USE kind_module
    !USE array_module, ONLY: nez, nezext
    USE physcnst_module, ONLY: kmev
    !USE nu_energy_grid_module, ONLY: unui, unubi, dunui
    !USE ec_table_module

    IMPLICIT NONE

    !-----------------------------------------------------------------------
    !  Output variables
    !-----------------------------------------------------------------------

    INTEGER, intent(in) :: nez
    REAL(dp), dimension(nez), intent(out) :: jec

    !-----------------------------------------------------------------------
    !  Local variables
    !-----------------------------------------------------------------------

    REAL(dp), dimension(nez), intent(in)   :: e_in, dunui
    REAL(dp), dimension(nez+1), intent(in) :: unubi

    REAL(dp) :: temp, rho, ye
    REAL(dp), DIMENSION(0:npts) :: jecfine
    REAL(dp) :: deltaupper(nez),deltalower(nez)
    REAL(dp) :: c1,c2,c3,c3base,loctot,tmev

    INTEGER :: kfmin(nez),kfmax(nez)
    INTEGER :: k, ii
    INTEGER :: rhoup,rhodown,tup,tdown
    INTEGER :: irho,iye,itemp
    INTEGER :: kmax,ktop,yeoffset(2)

    LOGICAL :: askew

    !REAL(dp), DIMENSION(nez) :: unubi, dunui

    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------

    kmax = nez
    jec  = 0.d0

    !unubi(1) = e_in(1)
    !dunui(1) = e_in(1)
    !do k=2,nez
    !  unubi(k) = e_in(k)
    !  dunui(k) = e_in(k) - e_in(k-1)
    !enddo

    !-----------------------------------------------------------------------
    !  Find temperature in table.  Temperature tabulated in MeV.
    !  For temperature outside of table, set rates to small value
    !-----------------------------------------------------------------------

    tmev  = temp*kmev
    itemp = int(tmev/deltaT_ec)
    tdown = itemp
    tup   = itemp+1

    IF (tup > tmx .or. tdown < tmn ) RETURN

    !-----------------------------------------------------------------------
    !  Find density in table.
    !  For density outside of table, set rates to small value
    !-----------------------------------------------------------------------

    irho=int(log10(rho)/deltalrho_ec)
    rhodown=irho
    rhoup = irho+1

    IF (rhoup > rhomx .or. rhodown < rhomn) RETURN

    !-----------------------------------------------------------------------
    !  Find ye in table.  Ye grid is density dependent
    !  For ye outside of table, allow extrapolation
    !-----------------------------------------------------------------------

    iye         = INT(Ye/deltaYe_ec)
    yeoffset(1) = iye-yefloor(rhodown)
    yeoffset(2) = iye-yefloor(rhoup)

    IF (sum(yeoffset) < -2) RETURN

    IF (sum(yeoffset) > (2*nye-2)) RETURN

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
    c1=(log10(rho)-dble(rhodown)*deltalrho_ec)
    !  /deltalrho_ec

    c2=(log10(tmev/(dble(itemp)*deltaT_ec)))
    ! c2=(log10(tmev/(dble(itemp)*deltaT_ec)))
    ! /(log10(dble(tup)/dble(tdown)))  !a deltaT_ec factor is ommitted from numerator and denominator

    c3base=(ye-dble(iye)*deltaYe_ec)/deltaYe_ec
    c3=c3base + &
    & c1*dble(yefloor(rhoup)+yeoffset(2)-yefloor(rhodown)-yeoffset(1)) !correction for skew If present

    !-----------------------------------------------------------------------
    ! Interpolate rate on table energy grid
    !-----------------------------------------------------------------------

    jecfine(:) = c3 & !note on indices: yeoffset=1 for rhodown, 2 for rhoup
    &           *( (1.0-c1)*(1.0-c2) * jectab(:,rhodown,tdown,yeoffset(1)+1) &  !   add one in upper block for old yeup
    &           +  c1*(1.0-c2)       * jectab(:,rhoup,tdown,yeoffset(2)+1) &
    &           +  c2*(1.0-c1)       * jectab(:,rhodown,tup,yeoffset(1)+1) &
    &           +  c1*c2             * jectab(:,rhoup,tup,yeoffset(2)+1) ) &
    &           + (1.0-c3) &
    &           *( (1.0-c1)*(1.0-c2) * jectab(:,rhodown,tdown,yeoffset(1)) &
    &           +  c1*(1.0-c2)       * jectab(:,rhoup,tdown,yeoffset(2)) &
    &           +  c2*(1.0-c1)       * jectab(:,rhodown,tup,yeoffset(1)) &
    &           +  c1*c2             * jectab(:,rhoup,tup,yeoffset(2)) )

    jecfine=10.0**(jecfine)

    !-----------------------------------------------------------------------
    ! Match Table grid to MGFLD energy grid
    !-----------------------------------------------------------------------

    ktop=kmax ! top energy bin that has a EC contribution
    DO k=1,kmax
      kfmin(k)=int(e_in(k)/deltaE)
      kfmax(k)=kfmin(k) +1
      IF (kfmax(k).gt.npts) kfmax(k)=npts
      IF (kfmin(k).gt.npts) kfmin(k)=npts
      IF (unubi(k+1) .gt. energy(npts)) ktop=min(k,ktop)
    ENDDO

    !-----------------------------------------------------------------------
    ! Integrate onto MGFLD energy grid
    !-----------------------------------------------------------------------

    DO k=1,ktop

      deltalower(k) = (e_in(k)-deltaE*dble(kfmin(k)))/deltaE
      deltaupper(k) = (deltaE*dble(kfmax(k))-e_in(k))/deltaE
      loctot        = dunui(k)*(deltaupper(k)*jecfine(kfmin(k))+deltalower(k)*jecfine(kfmax(k)))
      jec(k)        = loctot/(dunui(k)*e_in(k)**2)

    ENDDO

    RETURN

    END SUBROUTINE interp_ec

    SUBROUTINE abem_nuclei_EC_table_weaklib( n, e_in, unubi, dunui, rho, t, ye, xh, ah, cmpn, cmpp, cmpe, absrnc, emitnc, nse, nez )
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
    !  absrnc      :  absorption inverse mean free path (/cm)
    !  emitnc      :  emission inverse mean free path (/cm)
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
    REAL(dp), dimension(nez), intent(in)   :: e_in !neutrino energies
    REAL(dp), dimension(nez+1), intent(in) :: unubi !neutrino energies, bin edges
    REAL(dp), dimension(nez), intent(in)   :: dunui !neutrino energies, dE for edges

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

    REAL(dp), INTENT(out) :: absrnc(nez) ! inverse mean free path for absorption on free nucleons
    REAL(dp), INTENT(out) :: emitnc(nez) ! inverse mean free path for emission from free nucleons

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
    IF ( rho  >  roaenct  .or.  ah < 40 .or. nse == 0 ) THEN
      emitnc           = zero
      absrnc           = zero
      RETURN
    END IF

    !-----------------------------------------------------------------------
    !  Initialize
    !-----------------------------------------------------------------------

    tmev              = kmev * t
    xnuc              = ( xh/( rmu * ah ) ) * rho

    !-----------------------------------------------------------------------
    !
    !       \\\\\ E-NEUTRINO-NUCLEUS  EMISSION  AND ABSORPTION/////
    !
    !-----------------------------------------------------------------------
    !  Calculate values for rate/per heavy nucleus by interpolation
    !-----------------------------------------------------------------------

    CALL interp_ec( emitnc, e_in, unubi, dunui, t, rho, ye, nez )

    !-----------------------------------------------------------------------
    !  Multiply by heavy nucleus abundance
    !-----------------------------------------------------------------------

    emitnc = emitnc * xnuc

    !-----------------------------------------------------------------------
    !  Compute the absorption rate using detailed balance
    !-----------------------------------------------------------------------

    DO je = 1,nez
      absrnc(je) = emitnc(je) * fexp((e_in(je) + dmnp + cmpn - cmpp - cmpe)/tmev)
    END DO

    RETURN
    END SUBROUTINE abem_nuclei_EC_table_weaklib

END MODULE ec_table_module
