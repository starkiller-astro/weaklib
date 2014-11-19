PROGRAM table_get
!       
!    File:         table_get
!    Module:       table_get
!    Type:         Program
!
!    Date:         7/30/13
!
!    Purpose:
!        To create (rho, T, Ye) table of outputs from eos_get.f90
!        From these inputs and EOS variables, opacity kernels are called
!        and opacity tables are populated.    
!        
!        Coming updates:
!        -Extension to spanning all EOS space (and sewing the EOS space together)    
!        -Changing all module references to CHIMERA/GenASiS linked modules
!
!
!--------------------------------------------------
USE HDF5
USE kind_module, ONLY: double
USE edit_module, ONLY: nout, nlog
USE numerical_module, ONLY: zero, one, epsilon
USE physcnst_module, ONLY: kmev
USE ABEM_MODULE
USE EC_TABLE_MODULE 

IMPLICIT NONE

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

Integer                          :: nout2                   ! write index
Integer                          :: nout3                   ! write index
Integer                          :: nout4                   ! write index
Integer                          :: nout5                   ! write index
Integer                          :: nout6                   ! write index
Integer                          :: i                       ! rho index
Integer                          :: j                       ! T index
Integer                          :: k                       ! Ye index
Integer                          :: m                       ! nu-energy group index
Integer                          :: d                       ! quadrature index
Integer                          :: mbuff                   ! nu-energy group buffer index
Integer                          :: n                       ! # of data values per rho, T, Ye point
Integer, Parameter               :: nq = 3                  ! Degree of Gaussian Quadrature        
Integer, Parameter               :: Nrho = 2                ! Resolution in rho (# of grid points) 
Integer, Parameter               :: Nt = 2                  ! Resolution in T (# of grid points) 
Integer, Parameter               :: Nye = 2                 ! Resolution in ye (# of grid points) 
Integer, Parameter               :: nez = 20                ! # of neutrino energy groups
Integer, Parameter               :: nezbuff = (nez*nq)    ! # of neutrino energy buffer groups
Integer, Parameter               :: nse = 1                 ! NSE = 1, non-NSE = 0
!REAL(KIND=double), Parameter     :: rhomin = 1.000e+08     ! Minimum density, [g/cm**3] 
REAL(KIND=double), Parameter     :: rhomin = 1.000e+09      ! Minimum density, [g/cm**3] 
REAL(KIND=double), Parameter     :: rhomax = 1.000e+09      ! Maximum Density, [g/cm**3]
!REAL(KIND=double), Parameter     :: Tmin = 1.0d-01         ! Minimum Temperature, [MeV]
REAL(KIND=double), Parameter     :: Tmin = 4.5d00           ! Minimum Temperature, [MeV]
REAL(KIND=double), Parameter     :: Tmax = 4.5d00           ! Maximum Temperature, [MeV]
!EAL(KIND=double), Parameter     :: Yemin = 1.0d-02         ! Min Electron fraction   
REAL(KIND=double), Parameter     :: Yemin = 3.5d-01         ! Min Electron fraction   
REAL(KIND=double), Parameter     :: Yemax = 3.5d-01         ! Max Electron fraction
REAL(KIND=double), Parameter     :: Emin = 0.0d00           ! Minimum Neutrino Energy
REAL(KIND=double), Parameter     :: Emax = 3.0d02           ! Maximum Neutrino Energy
REAL(KIND=double), Parameter     :: alpha = 1.23902         ! Delta-E ratio
REAL(KIND=double), Parameter     :: delta_E_min = 1.0d00    ! Minimum Neutrino Energy bin width


REAL(double), DIMENSION(nezbuff)     :: E_quad              ! neutrino energy buffer from E_edge to E_quad; ebuff is gonna become E_quad [MeV]
!REAL(double), DIMENSION(nezbuff)     :: e_in                ! quadrature grid neutrino energy [MeV]
REAL(double), DIMENSION(nez)         :: delta_E             ! neutrino energy zone width [MeV]
REAL(double), DIMENSION(nez+1)       :: E_edge              ! neutrino energy zone edge [MeV]
REAL(double), DIMENSION(nez)         :: E_cent              ! neutrino energy zone center [MeV]
CHARACTER(LEN=15), PARAMETER  :: filename   = "opacitytable.h5"
CHARACTER(LEN=12), PARAMETER  :: dataset0   = "rhogrid_dset"
CHARACTER(LEN=10), PARAMETER  :: dataset00  = "Tgrid_dset"
CHARACTER(LEN=11), PARAMETER  :: dataset000 = "Yegrid_dset"
CHARACTER(LEN=32), PARAMETER  :: dataset1   = "Pressure                        "
CHARACTER(LEN=32), PARAMETER  :: dataset2   = "InternalEnergyDensity           "
CHARACTER(LEN=32), PARAMETER  :: dataset3   = "EntropyPerBaryon                "
CHARACTER(LEN=32), PARAMETER  :: dataset4   = "NeutronChemicalPotential        "
CHARACTER(LEN=32), PARAMETER  :: dataset5   = "ProtonChemicalPotential         "
CHARACTER(LEN=32), PARAMETER  :: dataset6   = "ElectronChemicalPotential       "
CHARACTER(LEN=32), PARAMETER  :: dataset7   = "NeutronMassFraction             "
CHARACTER(LEN=32), PARAMETER  :: dataset8   = "ProtonMassFraction              "
CHARACTER(LEN=32), PARAMETER  :: dataset9   = "HeliumMassFraction              "
CHARACTER(LEN=32), PARAMETER  :: dataset10  = "HeavyMassFraction               "
CHARACTER(LEN=32), PARAMETER  :: dataset11  = "HeavyMassNumber                 "
CHARACTER(LEN=32), PARAMETER  :: dataset12  = "HeavyChargeNumber               "
CHARACTER(LEN=32), PARAMETER  :: dataset13  = "Absorption                      "
CHARACTER(LEN=32), PARAMETER  :: dataset14  = "Emission                        "
CHARACTER(LEN=8), PARAMETER   :: dataset99  = "metadata"
CHARACTER(LEN=1), PARAMETER   :: L          = "L"              ! EoS flag

INTEGER, PARAMETER    :: dim00    = 12         ! EOS Fields dimension
INTEGER, PARAMETER    :: dim0     = 2          ! placeholder for Nrho
INTEGER, PARAMETER    :: dim1     = 2          ! placeholder for Nt
INTEGER, PARAMETER    :: dim2     = 2          ! placeholder for Nye 
INTEGER, PARAMETER    :: dim3     = nezbuff    ! placeholder for E-Nu bins , nezbuff when doing quad grid
INTEGER, PARAMETER    :: dim4     = 1          ! 
INTEGER, PARAMETER    :: dim5     = 15         ! # of metadata values 

INTEGER        :: hdferr                                 ! Error flag
INTEGER(HID_T) :: file, space0, space00, space000, space1, space2, space3 
INTEGER(HID_T) :: space4, space5, space6, space7, space8, space9, space10 
INTEGER(HID_T) :: space11, space12, space13, space14, space15, space99
INTEGER(HID_T) :: dset0, dset00, dset000 
INTEGER(HID_T) :: dset1, dset2, dset3, dset4, dset5, dset6, dset7, dset8, dset9
INTEGER(HID_T) :: dset10, dset11, dset12,  dset13, dset14, dset15, dset99
INTEGER(HSIZE_T), DIMENSION(1:3)   :: dims = (/dim0, dim1, dim2/) ! size individual write buffers
INTEGER(HSIZE_T), DIMENSION(1:4)   :: bdims = (/dim3, dim0, dim1, dim2/) ! size main read/write buffer
INTEGER(HSIZE_T), DIMENSION(1:2)   :: mddims = (/dim4, dim5/) ! size metadata write buffer

REAL, DIMENSION(1:dim0,1:dim1,1:dim2) :: rhogrid         ! Rho grid write buffer
REAL, DIMENSION(1:dim0,1:dim1,1:dim2) :: Tgrid           ! T grid write buffer
REAL, DIMENSION(1:dim0,1:dim1,1:dim2) :: Yegrid          ! Ye grid write buffer
REAL, DIMENSION(1:dim0,1:dim1,1:dim2) :: pressure        ! Pressure write buffer
REAL, DIMENSION(1:dim0,1:dim1,1:dim2) :: internalenergy  ! Internal Energy write buffer
REAL, DIMENSION(1:dim0,1:dim1,1:dim2) :: entropy         ! Entropy Write buffer
REAL, DIMENSION(1:dim0,1:dim1,1:dim2) :: neutchempot     ! Neutron chemical potential write buffer
REAL, DIMENSION(1:dim0,1:dim1,1:dim2) :: protchempot     ! Proton chemical [otential write buffer
REAL, DIMENSION(1:dim0,1:dim1,1:dim2) :: elecchempot     ! Electron chemical potential write buffer
REAL, DIMENSION(1:dim0,1:dim1,1:dim2) :: neutmassfrac    ! Neutron mass fraction write buffer
REAL, DIMENSION(1:dim0,1:dim1,1:dim2) :: protmassfrac    ! Proton mass fraction write buffer
REAL, DIMENSION(1:dim0,1:dim1,1:dim2) :: helmassfrac     ! Helium mass fraction write buffer
REAL, DIMENSION(1:dim0,1:dim1,1:dim2) :: heavymassfrac   ! Heavy mass fraction write buffer
REAL, DIMENSION(1:dim0,1:dim1,1:dim2) :: heavymassnum    ! Heavy mass number write buffer
REAL, DIMENSION(1:dim0,1:dim1,1:dim2) :: heavyqnum       ! Heavy charge # write buffer
REAL, DIMENSION(1:dim3,1:dim0,1:dim1,1:dim2) :: absorption ! Absorption inverse mean free path write buffer, all energy groups
REAL, DIMENSION(1:dim3,1:dim0,1:dim1,1:dim2) :: emission   ! Emission inverse mean free path write buffer, all energy groups
REAL, DIMENSION(1:dim4,1:dim5) :: metadata                ! Metadata write buffer

REAL(KIND=double), DIMENSION(Nrho*Nt*Nye)     :: rho      ! density [g/cm**3]
REAL(KIND=double), DIMENSION(Nrho*Nt*Nye)     :: ye       ! electron fraction
REAL(KIND=double), DIMENSION(Nrho*Nt*Nye)     :: t        ! temperature [K]
REAL(KIND=double)                             :: t_mev    ! temperature [MeV]
REAL(KIND=double), DIMENSION(nezbuff)         :: absornp  ! absorption inverse mean free path (/cm)
REAL(KIND=double), DIMENSION(nezbuff)         :: emitnp   ! emission inverse mean free path
REAL(KIND=double), DIMENSION(12, Nrho*Nt*Nye) :: Output   ! (rho, T, Ye, E, P, S,..)        
REAL(KIND=double), DIMENSION(nq)              :: wt       ! Gaussian Quadrature weighting
REAL(KIND=double), DIMENSION(nq)              :: x        ! Gaussian Quadrature points 
nout  = 3216
nout2 = 3217
nout3 = 3218
nout4 = 3219
nout5 = 3220
nout6 = 3221
nlog  = 4124
!------------------------------
! Initialize HDF5 FORTRAN interface.
!--------------------------------

  CALL h5open_f(hdferr)

!-----------------------------------------------------------------------
!  Do Loops to Create table
!-----------------------------------------------------------------------

OPEN(nout2, FILE="PressureOutput")
OPEN(nout, FILE="OutputFile")
DO k=0,Nye-1
  DO j=0,Nt-1
    DO i=0,Nrho-1
      rho(i+1 + j*Nrho + k*Nrho*Nt) = rhomin * 10.**((DBLE(i)/DBLE(Nrho))*LOG10(rhomax/rhomin)) 
      t_mev = (Tmin)*(10**((DBLE(j)/DBLE(Nt))*LOG10(Tmax/Tmin)))  
      t(i+1 + j*Nrho + k*Nrho*Nt) = t_mev / kmev
      ye(i+1 + j*Nrho + k*Nrho*Nt) = Yemin + ((j/(Nye-1))*(Yemax - Yemin))
      rhogrid(i+1, j+1, k+1) = rho(i+1 + j*Nrho + k*Nrho*Nt) 
      Tgrid(i+1, j+1, k+1) = t(i+1 + j*Nrho + k*Nrho*Nt) 
      Yegrid(i+1, j+1, k+1) = ye(i+1 + j*Nrho + k*Nrho*Nt) 

      END DO
   END DO
END DO
WRITE(*,*) "Sample rho t ye"

!-----------------------------------------------------------------------
!  Do Loops to Create table
!
!-----------------------------------------------------------------------

!DO k=0,50
!  DO j=0,50
!    DO i=0,80
        Call eos_get(Nrho, Nt, Nye, rho, t, ye, Output) ! Output to memory
        
!      END DO
!   END DO
!END DO

!-----------------------------------------------------------------------
!  Do Loops to Create table
!
!-----------------------------------------------------------------------

DO k=0,Nye-1  
  DO j=0,Nt-1  
    DO i=0,Nrho-1  
     pressure(i+1, j+1, k+1)=(Output(1, i+1+j*Nrho+k*Nrho*Nt))
     internalenergy(i+1, j+1, k+1)=(Output(2, i+1+j*Nrho+k*Nrho*Nt))
     entropy(i+1, j+1, k+1)=(Output(3, i+1+j*Nrho+k*Nrho*Nt))
     neutchempot(i+1, j+1, k+1)=(Output(4, i+1+j*Nrho+k*Nrho*Nt))
     elecchempot(i+1, j+1, k+1)=(Output(5, i+1+j*Nrho+k*Nrho*Nt))
     protchempot(i+1, j+1, k+1)=(Output(6, i+1+j*Nrho+k*Nrho*Nt))
     neutmassfrac(i+1, j+1, k+1)=(Output(7, i+1+j*Nrho+k*Nrho*Nt))
     protmassfrac(i+1, j+1, k+1)=(Output(8, i+1+j*Nrho+k*Nrho*Nt))
     helmassfrac(i+1, j+1, k+1)=(Output(9, i+1+j*Nrho+k*Nrho*Nt))
     heavymassfrac(i+1, j+1, k+1)=(Output(10, i+1+j*Nrho+k*Nrho*Nt)) 
     heavymassnum(i+1, j+1, k+1)=(Output(11, i+1+j*Nrho+k*Nrho*Nt))
     heavyqnum(i+1, j+1, k+1)=(Output(12, i+1+j*Nrho+k*Nrho*Nt))
       WRITE(nout2,"(es13.3)") (pressure(i+1, j+1, k+1))
       WRITE(nout,"(12es13.3)") (Output(n, i+1+j*Nrho+k*Nrho*Nt),n=1,12) ! to output file
    END DO
  END DO
END DO

!----------------------------------------------------------------------
! Fill holes due to lack of NR convergence
!
!----------------------------------------------------------------------
!DO k=0,Nye-1 
!  DO j=0,Nt-1 
!    DO i=0,Nrho-1
!      IF (pressure(i+1, j+1, k+1) < 0), THEN 
!       (pressure(i+1, j+1, k+1) = 1)
!       (internalenergy(i+1, j+1, k+1) = 1)
!       (entropy(i+1, j+1, k+1) = 1)
!       (neutchempot(i+1, j+1, k+1) = 1)
!       (elecchempot(i+1, j+1, k+1) = 1)
!       (protchempot(i+1, j+1, k+1) = 1)
!       (neutmassfrac(i+1, j+1, k+1) = 1)
!       (protmassfrac(i+1, j+1, k+1) = 1)
!       (helmassfrac(i+1, j+1, k+1) = 1)
!       (heavymassfrac(i+1, j+1, k+1) = 1)
!       (heavymassnum(i+1, j+1, k+1) = 1)
!       (heavyqnum(i+1, j+1, k+1)= 1 )
!      ELSE IF (pressure(i+1, j+1, k+1) = 1), THEN 
!      
!      END IF 
!       WRITE(nout2,"(es13.3)") (pressure(i+1, j+1, k+1))
!       WRITE(nout,"(12es13.3)") (Output(n, i+1+j*Nrho+k*Nrho*Nt),n=1,12) ! to output file
!    END DO
!  END DO
!END DO

!----------------------------------------------------------------------
!  Write fine energy grid for Gaussian quadrature
!----------------------------------------------------------------------

  CALL gquad(nq,x,wt,nq)

!----------------------------------------------------------------------
! Set up neutrino energy bins and add quadrature points
!----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Delta-E loop
!-----------------------------------------------------------------------
   delta_E(1) = delta_E_min

   DO m=2,nez
     delta_E(m) = (delta_E(m-1))*DBLE(alpha) ! Add Newton-Raphson solve of (alpha - 1)/(alpha^nez -1) = delta_E_min/Emax
   END DO

!-----------------------------------------------------------------------
! Energy Group Edge Loop
!-----------------------------------------------------------------------
   E_edge(1)=Emin

   DO m=2,nez+1
     E_edge(m) = E_edge(m-1) + delta_E(m-1)
   END DO

!-----------------------------------------------------------------------
! Energy Group Center Loop
!-----------------------------------------------------------------------

DO m=1,nez
   E_cent(m)=(E_edge(m+1) - E_edge(m))/2 +E_edge(m)
END DO

!-----------------------------------------------------------------------
! Quadrature Grid Loop
!-----------------------------------------------------------------------
DO m=1,nez
  DO d=1,nq
    mbuff = nq*(m-1) + d
    E_quad(mbuff) = E_edge(m) + ((E_edge(m+1)-E_edge(m))/2.d0)*(x(d)+1.d0) 
  END DO
END DO

!-----------------------------------------------------------------------
!  Do Loops to Populate Opacities
!  -Absorption/Emission
!-----------------------------------------------------------------------

DO k=0,Nye-1 
  DO j=0,Nt-1 
    DO i=0,Nrho-1  
   CALL abem_cal( 1, E_quad, rhogrid(i+1, j+1, k+1), Tgrid(i+1, j+1, k+1), &
 & neutmassfrac(i+1, j+1, k+1), protmassfrac(i+1, j+1, k+1), heavymassfrac(i+1, j+1, k+1), &
 & heavymassnum(i+1, j+1, k+1), heavyqnum(i+1, j+1, k+1), neutchempot(i+1, j+1, k+1), &
 & protchempot(i+1, j+1, k+1), elecchempot(i+1, j+1, k+1), absornp, emitnp, Tgrid(i+1, j+1, k+1), &
 & nezbuff, nse, L )

      DO m=0,nezbuff-1  
      absorption(m+1,i+1, j+1, k+1)=absornp(m+1)
      emission(m+1,i+1, j+1, k+1)=emitnp(m+1)
      END DO
    END DO
  END DO
END DO


!--------------------------------------------------------------------
! Do loop to write Metadata
!--------------------------------------------------------------------

! Nrho, Nt, Nye
! # of Eos Fields
! # of Opacity fields
! Max/Min of each rho, t, ye, E)
! Species, EOS type, etc.? 
! Metadata guide
  metadata(1,1) = Nrho
  metadata(1,2) = Nt  
  metadata(1,3) = Nye  
  metadata(1,4) = rhomin    ! Min rho
  metadata(1,5) = rhomax    !Max rho
  metadata(1,6) = Tmin      !Min T  
  metadata(1,7) = Tmax      !Max T  
  metadata(1,8) = Yemin     !Min Ye  
  metadata(1,9) = Yemax     !Max Ye
  metadata(1,10) = 12.d0    ! Number of EOS fields
  metadata(1,11) = 2.d0     !Number of Opacity fields
  metadata(1,12) = 1.d0     !Number of Neutrino Species
  metadata(1,13) = nezbuff  !Number of Energy groups
  metadata(1,14) = Emin                 
  metadata(1,15) = Emax


 
!--------------------------------------------------------------------
! Write Emissivity/Absorption vs. Energy
!--------------------------------------------------------------------

OPEN(nout3, FILE="AbsPlot.txt")
OPEN(nout4, FILE="EmPlot.txt")
OPEN(nout5, FILE="Grid.txt")
OPEN(nout6, FILE="QuadGrid.txt")

DO m=1,nezbuff
  WRITE(nout3,*) m, E_quad(m), absornp(m) 
  WRITE(nout4,*) m, E_quad(m), emitnp(m)
  WRITE(nout6,*) m, E_quad(m)
END DO

DO m=0,nez
  WRITE(nout5,*) m+1, E_edge(m+1), E_cent(m+1) 

END DO

CLOSE(nout3)
CLOSE(nout4)
CLOSE(nout5)
CLOSE(nout6)

CLOSE(nout)
CLOSE(nout2)

!-------------------------------------------------
! Create a new HDF5 file.
!-------------------------------------------------

   CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file, hdferr)

!---------------------------------------------------------
! Create dataspace.  Setting size to be the current size.
!---------------------------------------------------------

   CALL h5screate_simple_f(3, dims, space0, hdferr)
   CALL h5screate_simple_f(3, dims, space00, hdferr)
   CALL h5screate_simple_f(3, dims, space000, hdferr)
   CALL h5screate_simple_f(3, dims, space1, hdferr)
   CALL h5screate_simple_f(3, dims, space2, hdferr)
   CALL h5screate_simple_f(3, dims, space3, hdferr)
   CALL h5screate_simple_f(3, dims, space4, hdferr)
   CALL h5screate_simple_f(3, dims, space5, hdferr)
   CALL h5screate_simple_f(3, dims, space6, hdferr)
   CALL h5screate_simple_f(3, dims, space7, hdferr)
   CALL h5screate_simple_f(3, dims, space8, hdferr)
   CALL h5screate_simple_f(3, dims, space9, hdferr)
   CALL h5screate_simple_f(3, dims, space10, hdferr)
   CALL h5screate_simple_f(3, dims, space11, hdferr)
   CALL h5screate_simple_f(3, dims, space12, hdferr)
   CALL h5screate_simple_f(3, dims, space99, hdferr)
   CALL h5screate_simple_f(4, bdims, space13, hdferr)
   CALL h5screate_simple_f(4, bdims, space14, hdferr)
   CALL h5screate_simple_f(2, mddims, space99, hdferr)

!--------------------------------------------------------------------
! Create the datasets
!--------------------------------------------------------------------

   CALL h5dcreate_f(file, dataset0, H5T_NATIVE_DOUBLE, space0, dset0, hdferr)
   CALL h5dcreate_f(file, dataset00, H5T_NATIVE_DOUBLE, space00, dset00, hdferr)
   CALL h5dcreate_f(file, dataset000, H5T_NATIVE_DOUBLE, space000, dset000, hdferr)
   CALL h5dcreate_f(file, dataset1, H5T_NATIVE_DOUBLE, space1, dset1, hdferr)
   CALL h5dcreate_f(file, dataset2, H5T_NATIVE_DOUBLE, space2, dset2, hdferr)
   CALL h5dcreate_f(file, dataset3, H5T_NATIVE_DOUBLE, space3, dset3, hdferr)
   CALL h5dcreate_f(file, dataset4, H5T_NATIVE_DOUBLE, space4, dset4, hdferr)
   CALL h5dcreate_f(file, dataset5, H5T_NATIVE_DOUBLE, space5, dset5, hdferr)
   CALL h5dcreate_f(file, dataset6, H5T_NATIVE_DOUBLE, space6, dset6, hdferr)
   CALL h5dcreate_f(file, dataset7, H5T_NATIVE_DOUBLE, space7, dset7, hdferr)
   CALL h5dcreate_f(file, dataset8, H5T_NATIVE_DOUBLE, space8, dset8, hdferr)
   CALL h5dcreate_f(file, dataset9, H5T_NATIVE_DOUBLE, space9, dset9, hdferr)
   CALL h5dcreate_f(file, dataset10, H5T_NATIVE_DOUBLE, space10, dset10, hdferr)
   CALL h5dcreate_f(file, dataset11, H5T_NATIVE_DOUBLE, space11, dset11, hdferr)
   CALL h5dcreate_f(file, dataset12, H5T_NATIVE_DOUBLE, space12, dset12, hdferr)
   CALL h5dcreate_f(file, dataset13, H5T_NATIVE_DOUBLE, space13, dset13, hdferr)
   CALL h5dcreate_f(file, dataset14, H5T_NATIVE_DOUBLE, space14, dset14, hdferr)
   CALL h5dcreate_f(file, dataset99, H5T_NATIVE_DOUBLE, space99, dset99, hdferr)
!---------------------------------------------------------------
! Write the data to the dataset.
!---------------------------------------------------------------

   CALL h5dwrite_f(dset0, H5T_NATIVE_DOUBLE, rhogrid, dims, hdferr)
   CALL h5dwrite_f(dset00, H5T_NATIVE_DOUBLE, Tgrid, dims, hdferr)
   CALL h5dwrite_f(dset000, H5T_NATIVE_DOUBLE, Yegrid, dims, hdferr)
   CALL h5dwrite_f(dset1, H5T_NATIVE_DOUBLE, pressure, dims, hdferr)
   CALL h5dwrite_f(dset2, H5T_NATIVE_DOUBLE, internalenergy, dims, hdferr)
   CALL h5dwrite_f(dset3, H5T_NATIVE_DOUBLE, entropy, dims, hdferr)
   CALL h5dwrite_f(dset4, H5T_NATIVE_DOUBLE, neutchempot, dims, hdferr)
   CALL h5dwrite_f(dset5, H5T_NATIVE_DOUBLE, protchempot, dims, hdferr)
   CALL h5dwrite_f(dset6, H5T_NATIVE_DOUBLE, elecchempot, dims, hdferr)
   CALL h5dwrite_f(dset7, H5T_NATIVE_DOUBLE, neutmassfrac, dims, hdferr)
   CALL h5dwrite_f(dset8, H5T_NATIVE_DOUBLE, protmassfrac, dims, hdferr)
   CALL h5dwrite_f(dset9, H5T_NATIVE_DOUBLE, helmassfrac, dims, hdferr)
   CALL h5dwrite_f(dset10, H5T_NATIVE_DOUBLE, heavymassfrac, dims, hdferr)
   CALL h5dwrite_f(dset11, H5T_NATIVE_DOUBLE, heavymassnum, dims, hdferr)
   CALL h5dwrite_f(dset12, H5T_NATIVE_DOUBLE, heavyqnum, dims, hdferr)
   CALL h5dwrite_f(dset13, H5T_NATIVE_DOUBLE, absorption, bdims, hdferr)
   CALL h5dwrite_f(dset14, H5T_NATIVE_DOUBLE, emission, bdims, hdferr)
   CALL h5dwrite_f(dset99, H5T_NATIVE_DOUBLE, metadata, mddims, hdferr)
!---------------------------------------------------------------
! Close and release resources.
!----------------------------------------------------------------

   CALL h5dclose_f(dset0, hdferr)
   CALL h5dclose_f(dset00, hdferr)
   CALL h5dclose_f(dset000, hdferr)
   CALL h5dclose_f(dset1, hdferr)
   CALL h5dclose_f(dset2, hdferr)
   CALL h5dclose_f(dset3, hdferr)
   CALL h5dclose_f(dset4, hdferr)
   CALL h5dclose_f(dset5, hdferr)
   CALL h5dclose_f(dset6, hdferr)
   CALL h5dclose_f(dset7, hdferr)
   CALL h5dclose_f(dset8, hdferr)
   CALL h5dclose_f(dset9, hdferr)
   CALL h5dclose_f(dset10, hdferr)
   CALL h5dclose_f(dset11, hdferr)
   CALL h5dclose_f(dset12, hdferr)
   CALL h5dclose_f(dset13, hdferr)
   CALL h5dclose_f(dset14, hdferr)
   CALL h5dclose_f(dset99, hdferr)
   CALL h5sclose_f(space0, hdferr)
   CALL h5sclose_f(space00, hdferr)
   CALL h5sclose_f(space000, hdferr)
   CALL h5sclose_f(space1, hdferr)
   CALL h5sclose_f(space2, hdferr)
   CALL h5sclose_f(space3, hdferr)
   CALL h5sclose_f(space4, hdferr)
   CALL h5sclose_f(space5, hdferr)
   CALL h5sclose_f(space6, hdferr)
   CALL h5sclose_f(space7, hdferr)
   CALL h5sclose_f(space8, hdferr)
   CALL h5sclose_f(space9, hdferr)
   CALL h5sclose_f(space10, hdferr)
   CALL h5sclose_f(space11, hdferr)
   CALL h5sclose_f(space12, hdferr)
   CALL h5sclose_f(space13, hdferr)
   CALL h5sclose_f(space14, hdferr)
   CALL h5sclose_f(space99, hdferr)
   CALL h5fclose_f(file , hdferr)




!-----------------------------------------------------------------------
!  Edit
!-----------------------------------------------------------------------

!  WRITE (nprint)

STOP
END PROGRAM table_get
