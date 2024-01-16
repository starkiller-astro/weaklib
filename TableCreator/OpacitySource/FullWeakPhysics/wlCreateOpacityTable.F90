PROGRAM wlCreateOpacityTable
!---------------------------------------------------------------------
!
!    Author:       V. Mewes, ORNL
!
!    Created:      08/02/21
!
!    WeakLib ver:
!		   08/02/21	V. Mewes
!
!    Purpose:
!      Create table for neutrino opacities using state-of-art weak
!      physics.
!      Existing EoS table is readin.
!      Function were created to fill OpacityTable and Chimera
!      routines were called.
!
!---------------------------------------------------------------------
!                         Four Opacity Types
!
! OpacityType A for emission and absorption EmAb( rho, T, Ye, E)
!
!       e- + p/A <--> v_e + n/A*
!       v + n <--> e- + p
!       vb + p <--> e+ + n
!       vb + p + e- <--> n
!
!
! OpacityType B for isoenergetic scattering Iso( e, l, rho, T, Ye )
!
!       v_i/anti(v_i) + A --> v_i/anti(v_i) + A
!       v_i/anti(v_i) + e+/e-/n/p  <-->  v_i/anti(v_i) + e+/e-/n/p
!
! OpacityType C for non-iso scattering NES/Pair( e_in, e_out, l, T, eta)
!
!       e+ + e-  <--> v_i + anti(v_i);   i=e, muon, tau
!       N + N   <--> N + N + v_i + anti(v_i)
!
! OpacityType D for nucleon-nucleon Bremsstrahlung( e_p, e, rho, T)
!
! The full spin-density autocorrelation function S_sigma(eps+eps') 
! in units [1/MeV] from Hannestad and Raffelt 1998
! for a generic rho. The composition can later be taken into account
! by interpolating and adding the contribution from
! xp*rho, xn*rho and sqrt(xp+xn)*rho
!
!       nu nubar NN <--> NN   i=e, mu, tau
!---------------------------------------------------------------------


  USE wlKindModule, ONLY: dp
  USE wlInterpolationUtilitiesModule, ONLY: &
      Index1D_Lin, &
      GetIndexAndDelta_Lin, &
      GetIndexAndDelta_Log
  USE wlGridModule, ONLY: &
      MakeLogGrid,        &
      MakeLinearGrid,     &
      MakeGeometricGrid
  USE wlIOModuleHDF, ONLY: &
      InitializeHDF,       &
      FinalizeHDF
  USE wlOpacityTableModule, ONLY: &
      OpacityTableType,     &
      AllocateOpacityTable, &
      DescribeOpacityTable, &
      DeAllocateOpacityTable
  USE wlOpacityTableIOModuleHDF, ONLY: &
      WriteOpacityTableHDF
  USE wlExtPhysicalConstantsModule, ONLY: kMeV
  USE wlExtNumericalModule, ONLY: epsilon
  USE HR98_Bremsstrahlung
  USE prb_cntl_module, ONLY: &
      i_aeps, iaefnp, rhoaefnp, iaence, iaenct, roaenct, &
      edmpa, edmpe, iaenca, in, ip, ietann

  USE abem_module_weaklib, ONLY: &
      nleg_a, x_a, wt_a, &
      nleg_e, x_e, wt_e

  USE CC_module_weaklib

  USE ec_table_module

  USE, INTRINSIC :: iso_fortran_env, only : stdin=>input_unit, &
                                            stdout=>output_unit, &
                                            stderr=>error_unit
  USE omp_lib

  USE wlExtEOSWrapperModule, ONLY: &
      wlGetFullEOS, wlGetElectronEOS

IMPLICIT NONE

   INTEGER*4 today(3), now(3)

!---------------------------------------------------------------------
! Set Eos table
!---------------------------------------------------------------------
   CHARACTER(256) :: EOSTableName = "wl-EOS-SFHo-15-25-50.h5"
   !CHARACTER(256) :: EOSTableName = "wl-EOS-SFHo-25-50-100.h5"
   !CHARACTER(256) :: EOSTableName = "wl-EOS-LS220-25-50-100.h5"
   !CHARACTER(256) :: EOSTableName = "wl-EOS-LS220-15-25-50-Lower-T.h5"

!---------------------------------------------------------------------
! Set neutrino interaction type
!---------------------------------------------------------------------
   CHARACTER(256)          :: WriteTableName
   CHARACTER(256)          :: WriteTableNameBrem

   TYPE(OpacityTableType)  :: OpacityTable
   INTEGER, PARAMETER      :: nOpac_EmAb = 0  ! 2 for electron type

   INTEGER, PARAMETER      :: EmAb_np_FK &
                              = 0 
                              !Use Fischer et al 2020 full kinematics rates for
                              !EmAb on free nucleons

   INTEGER, PARAMETER      :: EmAb_np_FK_inv_n_decay &
                              = 0 
                              !Use Fischer et al 2020 full kinematics rates for
                              !inverse neutron decay 

   INTEGER, PARAMETER      :: EmAb_np_isoenergetic &
                              = 0
                              !EmAb on free nucleons using isoenergetic approximation
                              !Bruenn 1985
                              !Mezzacappa & Bruenn (1993)

   INTEGER, PARAMETER      :: EmAb_np_non_isoenergetic &
                              = 1
                              !EmAb on free nucleons taking into account recoil,
                              !nucleon final-state blocking, and special relativity
                              !Reddy et al 1998
                              !Only used for rho > 1e9, for low densities 
                              !the isoenergetic approximation of Bruenn 1985 
                              !and Mezzacappa & Bruenn (1993) is used

   INTEGER, PARAMETER      :: EmAb_np_weak_magnetism &
                              = 1
                              !Weak magnetism corrections for EmAb on free nucleons
                              !Horowitz 2002

   INTEGER, PARAMETER      :: EmAb_nuclei_EC_FFN &
                              = 0
                              !EmAb on nuclei using the Fuller, Fowler, 
                              !Neuman 1982, Ap. J. 252, 715 approximation for the 
                              !heavy nucleus matrix element as given in Bruenn 1985

   INTEGER, PARAMETER      :: EmAb_nuclei_EC_table &
                              = 1 
                              !EmAb on nuclei using a NSE-folded LMSH table
                              !Langanke et al. (2003), Hix et al. (2003)
   !tabulated EC table range
   REAL(dp), PARAMETER     :: EC_rho_min =   1.0d8
   REAL(dp), PARAMETER     :: EC_rho_max =   1.0d13
   REAL(dp), PARAMETER     :: EC_T_min   =   0.4d0 / kMeV
   REAL(dp), PARAMETER     :: EC_T_max   =   4.0d0 / kMeV
   REAL(dp), PARAMETER     :: EC_Ye_min  =   0.24d0
   REAL(dp), PARAMETER     :: EC_Ye_max  =   0.58d0
   REAL(dp), PARAMETER     :: EC_E_min   =   0.0d0
   REAL(dp), PARAMETER     :: EC_E_max   = 100.0d0

   INTEGER, PARAMETER      :: nOpac_Iso  = 2  ! 2 for electron type
                                              !   ( flavor identical )
   INTEGER, PARAMETER      :: nMom_Iso   = 2  ! 2 for 0th & 1st order
                                              !   legendre coff.
   INTEGER, PARAMETER      :: Iso_weak_magnetism &
                              = 1
                              !Weak magnetism corrections for isoenergetic scattering
                              !Horowitz 2002

   INTEGER, PARAMETER      :: Iso_ion_ion_corrections &
                              = 1
                              !Ion-ion correlation corrections to isoenergetic scattering
                              !Horowitz 1997, Bruenn and Mezzacappa 1997

   INTEGER, PARAMETER      :: Iso_many_body_corrections &
                              = 1
                              !Modification to neutral current scattering due to many-body effects
                              !Horowitz et al 2017

   REAL(DP), PARAMETER     :: Iso_ga_strange &
                              = -0.1d0
                              !Include strange-quark contributions to the axial vector coupling constant ga
                              !Value from Hobbs et al 2016

   INTEGER, PARAMETER      :: nOpac_NES  = 0  ! 1 ( either 0 or 1 )
   INTEGER, PARAMETER      :: nMom_NES   = 4  ! 4 for H1l, H2l
                                              !   ( either 0 or 4 )

   INTEGER, PARAMETER      :: NPS        = 1  !Include neutrino-positron scattering as well

   INTEGER, PARAMETER      :: nOpac_Pair = 0  ! 1 ( either 0 or 1 )
   INTEGER, PARAMETER      :: nMom_Pair  = 4  ! 4 for J1l, J2l
                                              !   ( either 0 or 4 )

   INTEGER, PARAMETER      :: nOpac_Brem = 0  !Only S_sigma(eps+eps') is needed for all
   INTEGER, PARAMETER      :: nMom_Brem  = 1  !species and moments

!---------------------------------------------------------------------
! Set E grid limits
!---------------------------------------------------------------------
   INTEGER,  PARAMETER     :: nPointsE = 40
   REAL(dp), PARAMETER     :: Emin     = 1.0d-01   !lower face of first energy cell
!   REAL(dp), PARAMETER     :: Emin = 0.0d+00 !lower face of first energy cell
   REAL(dp), PARAMETER     :: Emax = 3.0d+02 !upper face of last energy cell

!---------------------------------------------------------------------
! Set Eta grid limits
!---------------------------------------------------------------------
   INTEGER                 :: nPointsEta = 120
   REAL(dp), PARAMETER     :: Etamin     = 1.0d-3
   REAL(dp), PARAMETER     :: Etamax     = 2.5d03

   ! --- other inner variables
   INTEGER                 :: i_r, i_rb, i_e, j_rho, k_t, l_ye, &
                              t_m, i_eta, i_ep, stringlength, stringlengthBrem
   REAL(dp)                :: rho, T, TMeV, ye, Z, A, &
                              chem_e, chem_n, chem_p, xheavy, xn, &
                              xp, xhe, bb, eta, minvar

   REAL(dp), DIMENSION(nPointsE,2) :: cok
   REAL(dp), DIMENSION(nPointsE, nPointsE) :: H0i, H0ii, H1i, H1ii
   REAL(dp)                                :: j0i, j0ii, j1i, j1ii

   REAL(dp), DIMENSION(nPointsE, nPointsE) :: S_sigma !the Bremsstrahlung annihilation kernel

   REAL(dp) :: total_fraction, press, ent, eps, GAMMA1, el_press_e, el_entrop_e, el_energ_e, el_chem_e

   CALL idate(today)
   CALL itime(now)
   WRITE ( *, 10000 )  today(2), today(1), today(3), now

   10000 format ( ' Date ', i2.2, '/', i2.2, '/', i4.4, '; Time ',&
     &         i2.2, ':', i2.2, ':', i2.2 )

#if defined(GIT_HASH)
        write(*,*) 'Git version: ', GIT_HASH
#else
        write(*,*) 'Git version: unknown'
#endif

#if defined(GIT_BRANCH)
        write(*,*) 'Git branch: ', GIT_BRANCH
#else
        write(*,*) 'Git branch: unknown'
#endif

#if defined(GIT_DATE)
        write(*,*) 'Git date: ', GIT_DATE
#else
        write(*,*) 'Git date: unknown'
#endif


! -- Set Write-Out FileName
   stringlength = LEN(TRIM(EOSTableName))
   WRITE(WriteTableName,'(A,A,A,A,I2)') &
         EOSTableName(1:3),'Op-',EOSTableName(8:stringlength-3),'-E',nPointsE
   stringlength = len(TRIM(WriteTableName))

   stringlengthBrem = LEN(TRIM(EOSTableName))
   WRITE(WriteTableNameBrem,'(A,A,A,A,I2)') &
         EOSTableName(1:3),'Op-',EOSTableName(8:stringlengthBrem-3),'-E',nPointsE

   stringlengthBrem = len(TRIM(WriteTableNameBrem))

   CALL InitializeHDF( )

   !This needs to be here so that we allocate the EC table 
   !spectrum and rate arrays in EmAb
   OpacityTable % EmAb % nuclei_EC_table = EmAb_nuclei_EC_table

   CALL AllocateOpacityTable &
            ( OpacityTable, nOpac_EmAb, nOpac_Iso, nMom_Iso, &
              nOpac_NES, nMom_NES, nOpac_Pair, nMom_Pair, &
              nOpac_Brem, nMom_Brem, &
              nPointsE, nPointsEta, &
              EquationOfStateTableName_Option = EOSTableName )
   CALL FinalizeHDF( )

! -- Set OpacityTableTypeEmAb

   IF( nOpac_EmAb .gt. 0 ) THEN

   if(EmAb_np_FK + EmAb_np_isoenergetic .gt. 1) then
     write (*,*) "Must choose either full kinematics or Bruenn85 for EmAb on free nucleons"
     return
   endif
   if(EmAb_np_FK + EmAb_np_non_isoenergetic .gt. 1) then
     write (*,*) "Must choose either full kinematics or Reddy et al 98 for EmAb on free nucleons"
     return
   endif
   if(EmAb_np_isoenergetic + EmAb_np_non_isoenergetic .gt. 1) then
     write (*,*) "Must choose either Bruenn85 or Reddy et al 98 for EmAb on free nucleons"
     write (*,*) "When choosing Reddy et al 98, Bruenn85 is still used for rho < 1e9!" 
     return
   endif
   if(EmAb_nuclei_EC_FFN + EmAb_nuclei_EC_table .gt. 1) then
     write (*,*) "Must choose either Bruenn85-FFN formalism or EC table for EmAb on nuclei"
     return
   endif

   OpacityTable % EmAb % nOpacities   = nOpac_EmAb

   OpacityTable % EmAb % nPoints(1)   = nPointsE

   OpacityTable % EmAb % nPoints(2:4) = OpacityTable % nPointsTS

   OpacityTable % EmAb % Names = &
                                          (/'Electron Neutrino           ',  &
                                            'Electron Antineutrino       '/)
   OpacityTable % EmAb % EC_table_Names = (/'Electron Neutrino           '/)

   OpacityTable % EmAb % Units = &
                                          (/'Per Centimeter              ',  &
                                            'Per Centimeter              '/)
   OpacityTable % EmAb % EC_table_Units = &
                                        (/'Per Centimeter              '/)

   OpacityTable % EmAb % np_FK               = EmAb_np_FK
   OpacityTable % EmAb % np_FK_inv_n_decay   = EmAb_np_FK_inv_n_decay
   OpacityTable % EmAb % np_isoenergetic     = EmAb_np_isoenergetic
   OpacityTable % EmAb % np_non_isoenergetic = EmAb_np_non_isoenergetic
   OpacityTable % EmAb % np_weak_magnetism   = EmAb_np_weak_magnetism
   OpacityTable % EmAb % nuclei_EC_FFN       = EmAb_nuclei_EC_FFN

   OpacityTable % EmAb % EC_table_rho_min    = EC_rho_min
   OpacityTable % EmAb % EC_table_rho_max    = EC_rho_max
   OpacityTable % EmAb % EC_table_T_min      = EC_T_min
   OpacityTable % EmAb % EC_table_T_max      = EC_T_max
   OpacityTable % EmAb % EC_table_Ye_min     = EC_Ye_min
   OpacityTable % EmAb % EC_table_Ye_max     = EC_Ye_max
   OpacityTable % EmAb % EC_table_E_min      = EC_T_min
   OpacityTable % EmAb % EC_table_E_max      = EC_T_max

   END IF
! -- Set OpacityTableTypeScat Iso
   IF( nOpac_Iso .gt. 0 ) THEN

   OpacityTable % Scat_Iso % nOpacities   = nOpac_Iso

   OpacityTable % Scat_Iso % nMoments     = nMom_Iso

   OpacityTable % Scat_Iso % nPoints(1)   = nPointsE
   OpacityTable % Scat_Iso % nPoints(2)   = nMom_Iso
   OpacityTable % Scat_Iso % nPoints(3:5) = OpacityTable % nPointsTS

   OpacityTable % Scat_Iso % Names = &
                                (/'Electron Neutrino           ',  &
                                  'Electron Antineutrino       '/)

   OpacityTable % Scat_Iso % Units = &
                                (/'Per Centimeter              ',  &
                                  'Per Centimeter              '/)

   OpacityTable % Scat_Iso % weak_magnetism_corrections = &
                  Iso_weak_magnetism

   OpacityTable % Scat_Iso % ion_ion_corrections = &
                  Iso_ion_ion_corrections

   OpacityTable % Scat_Iso % many_body_corrections = &
                  Iso_many_body_corrections

   OpacityTable % Scat_Iso % ga_strange = &
                  Iso_ga_strange

   END IF

! -- Set OpacityTableTypeScat NES
   IF( nOpac_NES .gt. 0 ) THEN

   OpacityTable % Scat_NES % nOpacities = nOpac_NES

   OpacityTable % Scat_NES % nMoments   = nMom_NES

   OpacityTable % Scat_NES % nPoints(1) = nPointsE
   OpacityTable % Scat_NES % nPoints(2) = nPointsE
   OpacityTable % Scat_NES % nPoints(3) = nMom_NES
   OpacityTable % Scat_NES % nPoints(4) = OpacityTable % nPointsTS(2)
   OpacityTable % Scat_NES % nPoints(5) = nPointsEta

   OpacityTable % Scat_NES % Names = &
                                (/'Kernels'/)

   OpacityTable % Scat_NES % Units = &
                                (/'Per Centimeter Per MeV^3'/)

   OpacityTable % Scat_NES % NPS = NPS

   END IF
! -- Set OpacityTableTypeScat Pair

   IF( nOpac_Pair .gt. 0 ) THEN

   OpacityTable % Scat_Pair % nOpacities = nOpac_Pair

   OpacityTable % Scat_Pair % nMoments   = nMom_Pair

   OpacityTable % Scat_Pair % nPoints(1) = nPointsE
   OpacityTable % Scat_Pair % nPoints(2) = nPointsE
   OpacityTable % Scat_Pair % nPoints(3) = nMom_Pair
   OpacityTable % Scat_Pair % nPoints(4) = OpacityTable % nPointsTS(2)
   OpacityTable % Scat_Pair % nPoints(5) = nPointsEta

   OpacityTable % Scat_Pair % Names = (/'Kernels'/)

   OpacityTable % Scat_Pair % Units = (/'Per Centimeter Per MeV^3'/)

   END IF

! -- Set OpacityTableTypeBrem Brem
   IF( nOpac_Brem .gt. 0 ) THEN

   OpacityTable % Scat_Brem % nOpacities = nOpac_Brem

   OpacityTable % Scat_Brem % nMoments   = nMom_Brem

   OpacityTable % Scat_Brem % nPoints(1) = nPointsE
   OpacityTable % Scat_Brem % nPoints(2) = nPointsE
   OpacityTable % Scat_Brem % nPoints(3) = nMom_Brem
   OpacityTable % Scat_Brem % nPoints(4) = OpacityTable % nPointsTS(OpacityTable % TS % Indices % iRho)
   OpacityTable % Scat_Brem % nPoints(5) = OpacityTable % nPointsTS(OpacityTable % TS % Indices % iT)

   OpacityTable % Scat_Brem % Names      = (/'S_sigma'/)

   OpacityTable % Scat_Brem % Units      = (/'Per MeV'/)

   END IF

!-----------------------------
! Generate E grid from limits
!-----------------------------
PRINT*, "Making Energy Grid ... "

   ASSOCIATE( EnergyGrid => OpacityTable % EnergyGrid )

   EnergyGrid % Name &
     = 'Comoving Frame Neutrino Energy'

   EnergyGrid % Unit &
     = 'MeV                           '

   EnergyGrid % MinValue = Emin
   EnergyGrid % MaxValue = Emax
   EnergyGrid % LogInterp = 1
!   EnergyGrid % Zoom = 1.26603816071016d0

   CALL MakeLogGrid &
          ( EnergyGrid % MinValue, EnergyGrid % MaxValue, &
            EnergyGrid % nPoints,  EnergyGrid % Values )
!   CALL MakeGeometricGrid &
!          ( EnergyGrid % MinValue, EnergyGrid % MaxValue, EnergyGrid % Zoom, &
!            EnergyGrid % nPoints, EnergyGrid % Values, EnergyGrid % Width, &
!            EnergyGrid % Edge )

   END ASSOCIATE ! EnergyGrid

!-----------------------------
! Generate Eta grid from limits
!-----------------------------
PRINT*, "Making Eta Grid ... "

   ASSOCIATE( EtaGrid => OpacityTable % EtaGrid )

   EtaGrid % Name &
     = 'Elect. Chem. Pot. / Temperature'

   EtaGrid % Unit &
     = 'DIMENSIONLESS                   '

   EtaGrid % MinValue = Etamin
   EtaGrid % MaxValue = Etamax
   EtaGrid % LogInterp = 1

   CALL MakeLogGrid &
          ( EtaGrid % MinValue, EtaGrid % MaxValue, &
            EtaGrid % nPoints, EtaGrid % Values )

   END ASSOCIATE ! EtaGrid

!---------------------------------------------------------------------
!              Fill OpacityTable
!---------------------------------------------------------------------

PRINT*, 'Filling OpacityTable ...'
   ASSOCIATE(  &
       iRho    => OpacityTable % TS % Indices % iRho                           , &
       nRho    => OpacityTable % nPointsTS(OpacityTable % TS % Indices % iRho) , &
       iT      => OpacityTable % TS % Indices % iT                             , &
       nT      => OpacityTable % nPointsTS(OpacityTable % TS % Indices % iT)   , &
       iYe     => OpacityTable % TS % Indices % iYe                            , &
       nYe     => OpacityTable % nPointsTS(OpacityTable % TS % Indices % iYe)  , &
       E_cells => OpacityTable % EnergyGrid % Values                           , &
       Indices => OpacityTable % EOSTable % DV % Indices                       , &
       DVOffs  => OpacityTable % EOSTable % DV % Offsets                       , &
       DVar    => OpacityTable % EOSTable % DV % Variables  )

!-----------------  ECAPEM --------------------
  IF( ( nOpac_EmAb + nOpac_Iso ) .gt. 0 ) THEN
    WRITE(*,*) 'Calculating EmAb and Elastic Scattering Kernel ...'

    EmAb_Scat_Iso: BLOCK
   
    REAL(dp), PARAMETER :: mass_e  =   0.510998950d0 !electron restmass
    REAL(dp), PARAMETER :: mass_mu = 105.6583755d0   !muon restmass

    REAL(dp), PARAMETER :: mass_p =  938.272046d0  !proton restmass
    REAL(dp), PARAMETER :: mass_n =  939.565379d0  !neutron restmass

    REAL(dp), DIMENSION(nPointsE) :: absor, emit

    REAL(dp), DIMENSION(nPointsE) :: em_nucleons,        ab_nucleons
    REAL(dp), DIMENSION(nPointsE) :: em_inv_n_decay,     ab_inv_n_decay
    REAL(dp), DIMENSION(nPointsE) :: em_nuclei_EC_FFN,      ab_nuclei_EC_FFN

    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: EOS_Un, EOS_Up !neutron and proton mean field potentials from nuclear EOS
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: EOS_massn, EOS_massp !effective neutron and proton masses from nuclear EOS

    REAL(dp) :: Un_loc, Up_loc
    REAL(dp) :: massn_loc, massp_loc

    !arrays for weak magnetism corrections
    REAL(dp), DIMENSION(nPointsE) :: xi_n_wm, xib_p_wm

    IF(EmAb_np_FK .gt. 0 .or. &
       EmAb_np_FK_inv_n_decay   .gt. 0) THEN

      BLOCK

        integer :: SFHo_EOS, LS220_EOS
        integer, parameter :: ls220 = 1
        integer, parameter :: sfho  = 2

        WRITE(*,*) 'Need to get Un and Up from EOS for full kinematics EmAb'

        SFHo_EOS  = index(EOSTableName,"SFHo")
        LS220_EOS = index(EOSTableName,"LS220")

        ALLOCATE(EOS_Un   (nRho, nT, nYe))
        ALLOCATE(EOS_Up   (nRho, nT, nYe))
        ALLOCATE(EOS_massn(nRho, nT, nYe))
        ALLOCATE(EOS_massp(nRho, nT, nYe))

        IF(SFHo_EOS .gt. 0) THEN

          WRITE(*,*) 'Interpolate Un and Up from sfho_frdm_v1.03'

          CALL Interp_UnUp("UnUp_LS220_SFHO.h5", &
                           OpacityTable % TS % States (iRho) % Values, &
                           OpacityTable % TS % States (iT) % Values,   &
                           OpacityTable % TS % States (iYe) % Values,  &
                           EOS_Un, EOS_Up, EOS_massn, EOS_massp, sfho)

        ELSE IF(LS220_EOS .gt. 0) THEN
 
          WRITE(*,*) 'Interpolate Un and Up from CompOSE_LS220wl'

          CALL Interp_UnUp("UnUp_LS220_SFHO.h5", &
                           OpacityTable % TS % States (iRho) % Values, &
                           OpacityTable % TS % States (iT) % Values,   &
                           OpacityTable % TS % States (iYe) % Values,  &
                           EOS_Un, EOS_Up, EOS_massn, EOS_massp, ls220)

        ENDIF

      END BLOCK

    ENDIF

   !if we are tabulating EmAb on nuclei with the Hix2003
   !formalism, read in EC capture table and construct the
   !face-centered neutrino energy grid
   IF ( EmAb_nuclei_EC_table .gt. 0) THEN
             
     iaenct = EmAb_nuclei_EC_table

     build_EC_table: BLOCK

       CHARACTER(len=1000)     :: ec_table
       INTEGER                 :: k     
       INTEGER                 :: iD_min, iT_min, iYe_min
       INTEGER                 :: iD_max, iT_max, iYe_max
       INTEGER                 :: jD, jT, jYe
       INTEGER                 :: iD_EOS, iT_EOS, iYe_EOS
       INTEGER                 :: nPointsD, nPointsT, nPointsYe
       INTEGER, PARAMETER      :: nE = npts + 1
       REAL(dp), DIMENSION(nE) :: EC_table_spec
       REAL(dp)                :: EC_table_rate
       
       REAL(dp)                :: dD, dT, dYe

       CALL GetIndexAndDelta_Lin( LOG10(OpacityTable % EmAb % EC_table_rho_min), &
                                  LOG10(OpacityTable % TS % States (iRho) % Values), iD_min, dD )
       CALL GetIndexAndDelta_Lin( LOG10(OpacityTable % EmAb % EC_table_T_min), &
                                  LOG10(OpacityTable % TS % States (iT) % Values),   iT_min, dT )
       CALL GetIndexAndDelta_Lin( OpacityTable % EmAb % EC_table_Ye_min, &
                                  OpacityTable % TS % States (iYe) % Values,  iYe_min, dYe )

       CALL GetIndexAndDelta_Lin( LOG10(OpacityTable % EmAb % EC_table_rho_max), &
                                  LOG10(OpacityTable % TS % States (iRho) % Values), iD_max, dD )
       CALL GetIndexAndDelta_Lin( LOG10(OpacityTable % EmAb % EC_table_T_max), &
                                  LOG10(OpacityTable % TS % States (iT) % Values),   iT_max, dT )
       CALL GetIndexAndDelta_Lin( OpacityTable % EmAb % EC_table_Ye_max, &
                                  OpacityTable % TS % States (iYe) % Values,  iYe_max, dYe )

       iD_min  = iD_min  - 1
       iT_min  = iT_min  - 1
       iYe_min = iYe_min - 1

       iD_max  = iD_max  + 1
       iT_max  = iT_max  + 1
       iYe_max = iYe_max + 1

       nPointsD  = iD_max  - iD_min  + 1
       nPointsT  = iT_max  - iT_min  + 1
       nPointsYe = iYe_max - iYe_min + 1

       OpacityTable % EmAb % EC_table_nRho = nPointsD
       OpacityTable % EmAb % EC_table_nT   = nPointsT
       OpacityTable % EmAb % EC_table_nYe  = nPointsYe
       OpacityTable % EmAb % EC_table_nE   = nE

       ALLOCATE(OpacityTable % EmAb % EC_table_rho(nPointsD) )
       ALLOCATE(OpacityTable % EmAb % EC_table_T  (nPointsT) )
       ALLOCATE(OpacityTable % EmAb % EC_table_Ye (nPointsYe))
       ALLOCATE(OpacityTable % EmAb % EC_table_E  (nE)       )

      DO i_r = 1, OpacityTable % EmAb % EC_table_nOpacities
        ALLOCATE( OpacityTable % EmAb % EC_table_spec(i_r) % Values &
                ( nPointsD, nPointsT, nPointsYe, nE) )
        ALLOCATE( OpacityTable % EmAb % EC_table_rate(i_r) % Values &
                ( nPointsD, nPointsT, nPointsYe) )
      ENDDO

       CALL get_environment_variable('WEAKLIB_DIR', ec_table)

       IF( LEN_TRIM(ec_table) .le. 0) THEN

         WRITE(*,*) 'Environment variable $WEAKLIB_DIR not set, please setup your weaklib environment!'
 
         STOP

       ENDIF

       ec_table = TRIM(ec_table) // TRIM("/TableCreator/OpacitySource/FullWeakPhysics/Executables/")

       ec_table = TRIM(ec_table) // TRIM("ec_table.d")

       CALL read_ec_table(ec_table)

       WRITE(*,*) 'Read in EC table from file ', ec_table 

       WRITE(*,*) 'Building EC table in weaklib'

       !Setup EC table matching energy grid
       OpacityTable % EmAb % EC_table_E(1) = 0.0d0
       DO i_e = 2, nE
         OpacityTable % EmAb % EC_table_E(i_e) = OpacityTable % EmAb % EC_table_E(i_e-1) + deltaE
       ENDDO

       DO l_ye = 1, nPointsYe

         iYe_EOS = l_ye + iYe_min - 1

         OpacityTable % EmAb % EC_table_Ye(l_ye) = OpacityTable % TS % States (iYe) % Values (iYe_EOS)
         ye = OpacityTable % EmAb % EC_table_Ye(l_ye)

         DO k_t = 1, nPointsT

           iT_EOS = k_t + iT_min - 1

           OpacityTable % EmAb % EC_table_T(k_t) = OpacityTable % TS % States (iT) % Values (iT_EOS)
           T = OpacityTable % EmAb % EC_table_T(k_t)
           TMeV = T * kMeV

           DO j_rho = 1, nPointsD

             iD_EOS = j_rho + iD_min - 1

             OpacityTable % EmAb % EC_table_rho(j_rho) = OpacityTable % TS % States (iRho) % Values (iD_EOS)
             rho = OpacityTable % EmAb % EC_table_rho(j_rho)

             chem_e = 10**DVar(Indices % iElectronChemicalPotential) % &
                      Values(iD_EOS, iT_EOS, iYe_EOS) &
                    - DVOffs(Indices % iElectronChemicalPotential)   &
                    - epsilon

             chem_p = 10**DVar(Indices % iProtonChemicalPotential) % &
                      Values(iD_EOS, iT_EOS, iYe_EOS) &
                    - DVOffs(Indices % iProtonChemicalPotential)   &
                    - epsilon

             chem_n = 10**DVar(Indices % iNeutronChemicalPotential) % &
                      Values(iD_EOS, iT_EOS, iYe_EOS) &
                    - DVOffs(Indices % iNeutronChemicalPotential)   &
                    - epsilon

                xp  = 10**DVar(Indices % iProtonMassFraction) % &
                      Values(iD_EOS, iT_EOS, iYe_EOS) &
                    - DVOffs(Indices % iProtonMassFraction)   &
                    - epsilon

                xn  = 10**DVar(Indices % iNeutronMassFraction) % &
                      Values(iD_EOS, iT_EOS, iYe_EOS) &
                    - DVOffs(Indices % iNeutronMassFraction)   &
                    - epsilon

                xhe  = 10**DVar(Indices % iAlphaMassFraction) % &
                       Values(iD_EOS, iT_EOS, iYe_EOS) &
                     - DVOffs(Indices % iAlphaMassFraction)   &
                     - epsilon

             xheavy  = 10**DVar(Indices % iHeavyMassFraction) % &
                       Values(iD_EOS, iT_EOS, iYe_EOS) &
                     - DVOffs(Indices % iHeavyMassFraction)   &
                     - epsilon

                 Z   = 10**DVar(Indices % iHeavyChargeNumber) % &
                       Values(iD_EOS, iT_EOS, iYe_EOS) &
                     - DVOffs(Indices % iHeavyChargeNumber)   &
                     - epsilon

                 A   = 10**DVar(Indices % iHeavyMassNumber) % &
                       Values(iD_EOS, iT_EOS, iYe_EOS) &
                     - DVOffs(Indices % iHeavyMassNumber)   &
                     - epsilon

                 bb  = (chem_e + chem_p - chem_n)/(T*kMev)


             EC_table_spec = 0.0d0
             EC_table_rate = 0.0d0

             CALL abem_nuclei_EC_table_weaklib ( 1, E_cells, rho, T, Ye, xheavy, A, chem_n, &
                                                 chem_p, chem_e, EC_table_spec, EC_table_rate, &
                                                 1, nE)
             DO i_e = 1, nE
               OpacityTable % EmAb % EC_table_spec(1) % &
                            Values (j_rho, k_t, l_ye, i_e) = EC_table_spec(i_e)
             END DO  !i_e

             OpacityTable % EmAb % EC_table_rate(1) % &
                          Values (j_rho, k_t, l_ye) = EC_table_rate

           END DO  !j_rho
         END DO  !k_t
       END DO  !l_ye

     END BLOCK build_EC_table

   ENDIF

   WRITE(*,*) 'Now building EmAb on nucleons'

   IF(EmAb_np_non_isoenergetic == 1) THEN

     CALL gquad(nleg_a,x_a,wt_a,nleg_a)
     CALL gquad(nleg_e,x_e,wt_e,nleg_e)

     CALL load_polylog_weaklib

   END IF

   !these are multiplicative corrections
   xi_n_wm  = 1.0d0
   xib_p_wm = 1.0d0

   IF(EmAb_np_weak_magnetism == 1) THEN
                 
     CALL cc_weak_mag_weaklib( E_cells, xi_n_wm, xib_p_wm, nPointsE )

   ENDIF

   DO l_ye = 1, nYe

     ye = OpacityTable % TS % States (iYe) % Values (l_ye)

     DO k_t = 1, nT

       T = OpacityTable % TS % States (iT) % Values (k_t)
       TMeV = T * kMeV

       DO j_rho = 1, nRho

         rho = OpacityTable % TS % States (iRho) % &
               Values (j_rho)

         chem_e = 10**DVar(Indices % iElectronChemicalPotential) % &
                  Values(j_rho, k_t, l_ye) &
                  - DVOffs(Indices % iElectronChemicalPotential)   &
                  - epsilon

         chem_p = 10**DVar(Indices % iProtonChemicalPotential) % &
                  Values(j_rho, k_t, l_ye) &
                  - DVOffs(Indices % iProtonChemicalPotential)   &
                  - epsilon

         chem_n = 10**DVar(Indices % iNeutronChemicalPotential) % &
                  Values(j_rho, k_t, l_ye) &
                  - DVOffs(Indices % iNeutronChemicalPotential)   &
                  - epsilon

            xp  = 10**DVar(Indices % iProtonMassFraction) % &
                  Values(j_rho, k_t, l_ye) &
                  - DVOffs(Indices % iProtonMassFraction)   &
                  - epsilon

            xn  = 10**DVar(Indices % iNeutronMassFraction) % &
                  Values(j_rho, k_t, l_ye) &
                  - DVOffs(Indices % iNeutronMassFraction)   &
                  - epsilon

           xhe  = 10**DVar(Indices % iAlphaMassFraction) % &
                  Values(j_rho, k_t, l_ye) &
                  - DVOffs(Indices % iAlphaMassFraction)   &
                  - epsilon

        xheavy  = 10**DVar(Indices % iHeavyMassFraction) % &
                  Values(j_rho, k_t, l_ye) &
                  - DVOffs(Indices % iHeavyMassFraction)   &
                  - epsilon

            Z   = 10**DVar(Indices % iHeavyChargeNumber) % &
                  Values(j_rho, k_t, l_ye) &
                  - DVOffs(Indices % iHeavyChargeNumber)   &
                  - epsilon

            A   = 10**DVar(Indices % iHeavyMassNumber) % &
                  Values(j_rho, k_t, l_ye) &
                  - DVOffs(Indices % iHeavyMassNumber)   &
                  - epsilon

        press   = 10**DVar(Indices % iPressure) % &
                  Values(j_rho, k_t, l_ye) &
                  - DVOffs(Indices % iPressure)   &
                  - epsilon

          ent   = 10**DVar(Indices % iEntropyPerBaryon) % &
                  Values(j_rho, k_t, l_ye) &
                  - DVOffs(Indices % iEntropyPerBaryon)   &
                  - epsilon

          eps   = 10**DVar(Indices % iInternalEnergyDensity) % &
                  Values(j_rho, k_t, l_ye) &
                  - DVOffs(Indices % iInternalEnergyDensity)   &
                  - epsilon

       GAMMA1   = 10**DVar(Indices % iGamma1) % &
                  Values(j_rho, k_t, l_ye) &
                  - DVOffs(Indices % iGamma1)   &
                  - epsilon

            bb  = (chem_e + chem_p - chem_n)/(T*kMev)

            IF(press < 0.0d0) WRITE(*,'(A22,es12.3,A2,es12.3,A2,f5.2,A2,es12.3)') &
                              'pressure is negative! ', rho, ' ',TMeV, ' ',ye, ' ',press

            IF(ent < 0.0d0) WRITE(*,'(A21,es12.3,A2,es12.3,A2,f5.2,A2,es12.3)') &
                              'entropy is negative! ', rho, ' ',TMeV, ' ',ye, ' ',ent

            IF(GAMMA1 < 0.0d0) WRITE(*,'(A29,es12.3,A2,es12.3,A2,f5.2,A2,es12.3)') &
                               'adiabatic index is negative! ', rho, ' ', TMeV, ' ', ye, ' ', GAMMA1

            IF(xn     < 0.0d0 .or. xn     > 1.0d0) write(*,'(A21,es12.3,A2,es12.3,A2,f5.2,A2,es15.3)') &
                                                   'xn either < 0 or > 1 ', rho, ' ',TMeV, ' ',Ye, ' ',xn
            IF(xp     < 0.0d0 .or. xp     > 1.0d0) write(*,'(A21,es12.3,A2,es12.3,A2,f5.2,A2,es15.3)') &
                                                   'xp either < 0 or > 1 ', rho, ' ',TMeV, ' ',Ye, ' ',xp
            IF(xhe    < 0.0d0 .or. xhe    > 1.0d0) write(*,'(A22,es12.3,A2,es12.3,A2,f5.2,A2,es15.3)') &
                                                   'xhe either < 0 or > 1 ', rho, ' ',TMeV, ' ',Ye, ' ',xhe
            IF(xheavy < 0.0d0 .or. xheavy > 1.0d0) write(*,'(A25,es12.3,A2,es12.3,A2,f5.2,A2,es15.3)') &
                                                   'xheavy either < 0 or > 1 ', rho, ' ',TMeV, ' ',Ye, ' ',xheavy

            total_fraction = xn+xp+xhe+xheavy
            IF(ABS(1.0d0-total_fraction) >= 1.0d-10) WRITE(*,'(A34,es12.3,A2,es12.3,A2,f5.2,A2,es15.3)') &
                                                     'mass fractions do not add up to 1 ', &
                                                     rho, ' ',TMeV, ' ',Ye, ' ',total_fraction

            if (     EmAb_np_FK == 1 &
                .or. EmAb_np_FK_inv_n_decay == 1 ) then
 
              Un_loc = EOS_Un(j_rho, k_t, l_ye)
              Up_loc = EOS_Up(j_rho, k_t, l_ye)

              massn_loc = EOS_massn(j_rho, k_t, l_ye)
              massp_loc = EOS_massp(j_rho, k_t, l_ye)

            endif

         DO i_r = 1, nOpac_EmAb

             iaefnp = EmAb_np_isoenergetic
             i_aeps = EmAb_np_non_isoenergetic
             rhoaefnp = HUGE(1.d0) ! (?)
             iaence = EmAb_nuclei_EC_FFN !e-nu
             edmpe = 3.d0
             iaenca = EmAb_nuclei_EC_FFN !e-nubar
             edmpa = 3.d0
             iaenct = EmAb_nuclei_EC_table

! Initialise the opacity arrays to zero
             ab_nucleons = 0.0d0
             em_nucleons = 0.0d0

             ab_inv_n_decay = 0.0d0
             em_inv_n_decay = 0.0d0

             ab_nuclei_EC_FFN = 0.0d0
             em_nuclei_EC_FFN = 0.0d0

             IF(EmAb_np_FK == 1) THEN

               WRITE(*,*) "EmAb on nucleons using full kinematics approach"

               IF(i_r <= 2) THEN
                 CALL CC_EmAb( OpacityTable % EnergyGrid % Values,     &
                               TMeV, mass_e, chem_e, chem_n+mass_n, chem_p+mass_p,      &
                               massn_loc, massp_loc, Un_loc, Up_loc,   &
                               i_r, EmAb_np_FK_inv_n_decay,       &
                               ab_nucleons, em_nucleons,               &
                               ab_inv_n_decay, em_inv_n_decay, nPointsE)

               ENDIF
             
             ELSE IF(EmAb_np_non_isoenergetic == 1) THEN

               IF(rho >= 1.0d9) THEN
                 CALL abem_np_non_isoenergetic_weaklib & 
                      ( i_r, E_cells, &
                        rho, T, ye, xn, xp, xheavy, A, Z, chem_n, chem_p, chem_e, &
                        ab_nucleons, em_nucleons, nPointsE)

               ELSE
                 CALL abem_nucleons_isoenergetic_weaklib & 
                      ( i_r, E_cells, &
                        rho, T, ye, xn, xp, xheavy, A, Z, chem_n, chem_p, chem_e, &
                        ab_nucleons, em_nucleons, nPointsE)
               END IF
             
               IF (i_r == 1) THEN !electron neutrino weak mag correction

                 ab_nucleons = ab_nucleons * xi_n_wm
                 em_nucleons = em_nucleons * xi_n_wm

               ELSE IF (i_r == 2) THEN !electron anti-neutrino weak mag correction

                 ab_nucleons = ab_nucleons * xib_p_wm
                 em_nucleons = em_nucleons * xib_p_wm

               ENDIF

             ELSE IF(EmAb_np_isoenergetic == 1) THEN

               CALL abem_nucleons_isoenergetic_weaklib & 
                    ( i_r, E_cells, &
                      rho, T, ye, xn, xp, xheavy, A, Z, chem_n, chem_p, chem_e, &
                      ab_nucleons, em_nucleons, nPointsE)

               IF (i_r == 1) THEN !electron neutrino weak mag correction

                 ab_nucleons = ab_nucleons * xi_n_wm
                 em_nucleons = em_nucleons * xi_n_wm

               ELSE IF (i_r == 2) THEN !electron anti-neutrino weak mag correction

                 ab_nucleons = ab_nucleons * xib_p_wm
                 em_nucleons = em_nucleons * xib_p_wm

               ENDIF

             ENDIF

             IF(EmAb_nuclei_EC_FFN == 1 .and. i_r == 1) THEN

               CALL abem_nuclei_EC_FFN_weaklib( E_cells, rho, T, xheavy, A, Z, chem_n, chem_p, &
                                             chem_e, ab_nuclei_EC_FFN, em_nuclei_EC_FFN, &
                                             nPointsE )

             ENDIF

             DO i_e = 1, OpacityTable % nPointsE
                OpacityTable % EmAb % Opacity(i_r) % &
                        Values (i_e, j_rho, k_t, l_ye) &
                = ab_nucleons(i_e)    + em_nucleons(i_e) &
                + ab_inv_n_decay(i_e) + em_inv_n_decay(i_e) &
                + ab_nuclei_EC_FFN(i_e)  + em_nuclei_EC_FFN(i_e) 
             END DO  !i_e
         END DO !i_r

!----------------  Scat_Iso -----------------------
         DO i_rb = 1, nOpac_Iso

           in = 1
           ip = 1
           ietann = 1

           CALL scatical_weaklib &
           ( i_rb, OpacityTable % EnergyGrid % Values,  &
             nPointsE, rho, T, ye, xn, xp, xhe, xheavy, A, Z, &
             Iso_weak_magnetism, Iso_ion_ion_corrections,     &
             Iso_many_body_corrections, Iso_ga_strange, cok )

           DO t_m = 1, nMom_Iso

             OpacityTable % Scat_Iso % Kernel(i_rb) % Values &
             ( :, t_m, j_rho, k_t, l_ye )  &
             = cok(:,t_m) * ( (t_m - 1) * 2.d0 + 1.d0 ) / 2.d0

           END DO !t_m

         END DO !i_rb

       END DO  !j_rho
     END DO  !k_t
   END DO  !l_ye

!------- EmAb % Offsets
   DO i_r = 1, nOpac_EmAb
     minvar = MINVAL( OpacityTable % EmAb % Opacity(i_r) % Values )
     OpacityTable % EmAb % Offsets(i_r) = -2.d0 * MIN( 0.d0, minvar )

     IF (i_r == 1 .and. EmAb_nuclei_EC_table == 1) THEN
       minvar = MINVAL( OpacityTable % EmAb % EC_table_spec(i_r) % Values )
       OpacityTable % EmAb % EC_table_spec_Offsets(i_r) = -2.d0 * MIN( 0.d0, minvar ) 
       minvar = MINVAL( OpacityTable % EmAb % EC_table_rate(i_r) % Values )
       OpacityTable % EmAb % EC_table_rate_Offsets(i_r) = -2.d0 * MIN( 0.d0, minvar ) 
     ENDIF
   END DO

!------- Scat_Iso % Offsets
   DO i_r = 1, nOpac_Iso
     DO t_m = 1, nMom_Iso
       minvar = MINVAL( OpacityTable % Scat_Iso % Kernel(i_r) &
                       % Values(:,t_m,:,:,:) )
       OpacityTable % Scat_Iso % Offsets(i_r, t_m) =           &
                      -2.d0 * MIN( 0.d0, minvar )
     END DO
   END DO

   END BLOCK EmAb_Scat_Iso

   END IF
!----------------  Scat_NES -----------------------
  IF( nOpac_NES .gt. 0 ) THEN
  PRINT*, 'Calculating Scat_NES Kernel ... '

      DO i_eta = 1, nPointsEta

        eta = OpacityTable % EtaGrid % Values(i_eta)

        DO k_t = 1, OpacityTable % nPointsTS(iT)

          T = OpacityTable % TS % States (iT) % Values (k_t)
          TMeV = T * kMeV

          CALL scatergn_weaklib &
               ( nPointsE, OpacityTable % EnergyGrid % Values, &
                 TMeV, eta, H0i, H0ii, H1i, H1ii, NPS )

          OpacityTable % Scat_NES % Kernel(1) % Values &
               ( :, :, 1, k_t, i_eta )       &
          = 0.5_DP * TRANSPOSE(H0i(:,:))  ! H0i was saved as H0i(e,ep)

          OpacityTable % Scat_NES % Kernel(1) % Values &
               ( :, :, 2, k_t, i_eta)       &
          = 0.5_DP * TRANSPOSE(H0ii(:,:)) ! H0ii was saved as H0ii(e,ep)

          OpacityTable % Scat_NES % Kernel(1) % Values &
               ( :, :, 3, k_t, i_eta )       &
          = 1.5_DP * TRANSPOSE(H1i(:,:))  ! H1i was saved as H1i(e,ep)
          OpacityTable % Scat_NES % Kernel(1) % Values &
               ( :, :, 4, k_t, i_eta )       &
          = 1.5_DP * TRANSPOSE(H1ii(:,:)) ! H1ii was saved as H1ii(e,ep)

        END DO  !k_t

      END DO !i_eta

!------- Scat_NES % Offsets
    DO i_rb = 1, nOpac_NES
      DO t_m = 1, nMom_NES

       minvar = MINVAL( OpacityTable % Scat_NES % Kernel(i_rb) &
                       % Values(:,:,t_m,:,: ) )
       OpacityTable % Scat_NES % Offsets(i_rb, t_m) =          &
                      -2.d0 * MIN( 0.d0, minvar )

    END DO ! t_m
  END DO ! i_rb
  END IF

!----------------  Scat_Pair -----------------------

  IF( nOpac_Pair .gt. 0 ) THEN
  PRINT*, 'Calculating Scat_Pair Kernel ... '

      DO i_eta = 1, nPointsEta

        eta = OpacityTable % EtaGrid % Values(i_eta)

        DO k_t = 1, OpacityTable % nPointsTS(iT)

          T = OpacityTable % TS % States (iT) % Values (k_t)
          TMeV = T * kMeV
          chem_e = TMeV * eta

          DO i_e = 1, nPointsE

            DO i_ep = 1, nPointsE

             CALL paircal_weaklib &
                  ( OpacityTable % EnergyGrid % Values(i_e),  &
                    OpacityTable % EnergyGrid % Values(i_ep), &
                    chem_e, T, j0i, j0ii, j1i, j1ii )

             OpacityTable % Scat_Pair % Kernel(1) % Values  &
                          ( i_ep, i_e, 1, k_t, i_eta )       &
              = 0.5_DP * j0i

             OpacityTable % Scat_Pair % Kernel(1) % Values  &
                          ( i_ep, i_e, 2, k_t, i_eta )       &
              = 0.5_DP * j0ii

             OpacityTable % Scat_Pair % Kernel(1) % Values  &
                          ( i_ep, i_e, 3, k_t, i_eta )       &
              = 1.5_DP * j1i

             OpacityTable % Scat_Pair % Kernel(1) % Values  &
                          ( i_ep, i_e, 4, k_t, i_eta )       &
              = 1.5_DP * j1ii

            END DO ! i_ep
          END DO ! i_e
        END DO  !k_t
      END DO !i_eta

!------- Scat_Pair % Offsets
  DO i_rb = 1, nOpac_Pair
    DO t_m = 1, nMom_Pair

       minvar = MINVAL( OpacityTable % Scat_Pair % Kernel(i_rb) &
                       % Values(:,:,t_m,:,:) )
       OpacityTable % Scat_Pair % Offsets(i_rb, t_m) =          &
                      -2.d0 * MIN( 0.d0, minvar )
    END DO ! t_m
  END DO ! i_rb
  END IF

!-----------------  Bremsstrahlung --------------------
  IF( nOpac_Brem .gt. 0 ) THEN
  PRINT*, 'Calculating Nucleon-Nucleon Bremsstrahlung Scattering Kernel ...'

   DO k_t = 1, OpacityTable % nPointsTS(iT)

    DO j_rho = 1, OpacityTable % nPointsTS(iRho)

      T = OpacityTable % EOSTable % TS % States (iT) % Values (k_t)
      rho = OpacityTable % EOSTable % TS % States (iRho) % Values (j_rho)

       CALL bremcal_weaklib &
            (nPointsE, OpacityTable % EnergyGrid % Values, &
             rho, T, S_sigma)

            OpacityTable % Scat_Brem % Kernel(1) % Values (:, :, 1, j_rho, k_t) &
            = TRANSPOSE(S_sigma(:,:))

    END DO  !j_rho
  END DO  !k_t

!------- Brem % Offsets

       minvar = MINVAL( OpacityTable % Scat_Brem % Kernel(1) % Values )
       OpacityTable % Scat_Brem % Offsets(1, 1) = -2.d0 * MIN( 0.d0, minvar )

  END IF

   END ASSOCIATE ! rho-T-Ye
!---------------------------------------------------------------------
!      Describe the Table ( give the real physical value )
!---------------------------------------------------------------------

  CALL DescribeOpacityTable( OpacityTable )

! --------------------------------------------------------------------
!          LOG the WHOLE table for storage
! --------------------------------------------------------------------

  WRITE(*,*) 'LOG the whole table with relevant offset for storage'

  DO i_r = 1, nOpac_EmAb
     OpacityTable % EmAb % Opacity(i_r) % Values&
     = LOG10( OpacityTable % EmAb % Opacity(i_r) % &
              Values + OpacityTable % EmAb % Offsets(i_r) + epsilon )
     IF (i_r == 1 .and. EmAb_nuclei_EC_table == 1) THEN
       OpacityTable % EmAb % EC_table_spec(i_r) % Values&
       = LOG10( OpacityTable % EmAb % EC_table_spec(i_r) % &
                Values + OpacityTable % EmAb % EC_table_spec_Offsets(i_r) + 1.0d-299 )
       OpacityTable % EmAb % EC_table_rate(i_r) % Values&
       = LOG10( OpacityTable % EmAb % EC_table_rate(i_r) % &
                Values + OpacityTable % EmAb % EC_table_rate_Offsets(i_r) + 1.0d-299 )
     ENDIF
  END DO  !i_r

  DO i_rb = 1, nOpac_Iso
    DO t_m = 1, nMom_Iso
      OpacityTable % Scat_Iso % Kernel(i_rb) % Values(:,t_m,:,:,:) &
      = LOG10 ( OpacityTable % Scat_Iso % Kernel(i_rb) % Values(:,t_m,:,:,:) &
                + OpacityTable % Scat_Iso % Offsets(i_rb,t_m) + epsilon )
    END DO ! t_m
  END DO ! i_rb

  DO i_rb = 1, nOpac_NES
    DO t_m = 1, nMom_NES
      OpacityTable % Scat_NES % Kernel(i_rb) % Values(:,:,t_m,:,:) &
      = LOG10 ( OpacityTable % Scat_NES % Kernel(i_rb) % Values(:,:,t_m,:,:) &
                + OpacityTable % Scat_NES % Offsets(i_rb, t_m) + epsilon )
    END DO
  END DO !i_rb

  DO i_rb = 1, nOpac_Pair
    DO t_m = 1, nMom_Pair
      OpacityTable % Scat_Pair % Kernel(i_rb) % Values(:,:,t_m,:,:) &
      = LOG10 ( OpacityTable % Scat_Pair % Kernel(i_rb) % Values(:,:,t_m,:,:) &
                + OpacityTable % Scat_Pair % Offsets(i_rb, t_m) + epsilon )
    END DO
  END DO !i_rb

  IF( nOpac_Brem > 0 ) THEN
        OpacityTable % Scat_Brem % Kernel(1) % Values &
      = LOG10 ( OpacityTable % Scat_Brem % Kernel(1) % Values &
              + OpacityTable % Scat_Brem % Offsets(1, 1) + epsilon )
  ENDIF

! -- write into hdf5 file

  IF( nOpac_EmAb > 0 ) THEN
    WriteTableName(stringlength+1:stringlength+8) = '-EmAb.h5'
    CALL InitializeHDF( )
    WRITE(*,*) 'Write EmAb data into file ', TRIM(WriteTableName)
    CALL WriteOpacityTableHDF &
         ( OpacityTable, TRIM(WriteTableName), WriteOpacity_EmAb_Option = .true. )
    CALL FinalizeHDF( )
  END IF

  IF( nOpac_Iso > 0 ) THEN
    WriteTableName(stringlength+1:stringlength+8) = '-Iso.h5 '
    CALL InitializeHDF( )
    WRITE(*,*) 'Write Iso data into file ', TRIM(WriteTableName)
    CALL WriteOpacityTableHDF &
         ( OpacityTable, TRIM(WriteTableName), WriteOpacity_Iso_Option = .true. )
    CALL FinalizeHDF( )
  END IF

  IF( nOpac_NES > 0 ) THEN
    WriteTableName(stringlength+1:stringlength+8) = '-NES.h5 '
    CALL InitializeHDF( )
    WRITE(*,*) 'Write NES data into file ', TRIM(WriteTableName)
    CALL WriteOpacityTableHDF &
         ( OpacityTable, TRIM(WriteTableName), WriteOpacity_NES_Option = .true. )
    CALL FinalizeHDF( )
  END IF

  IF( nOpac_Pair > 0 ) THEN
    WriteTableName(stringlength+1:stringlength+8) = '-Pair.h5'
    CALL InitializeHDF( )
    WRITE(*,*) 'Write Pair data into file ', TRIM(WriteTableName)
    CALL WriteOpacityTableHDF &
         ( OpacityTable, TRIM(WriteTableName), WriteOpacity_Pair_Option = .true. )
    CALL FinalizeHDF( )
  END IF

  IF( nOpac_Brem > 0 ) THEN
    WriteTableNameBrem(stringlengthBrem+1:stringlengthBrem+8) = '-Brem.h5'
    CALL InitializeHDF( )
    WRITE(*,*) 'Write Brem data into file ', TRIM(WriteTableNameBrem)
    CALL WriteOpacityTableHDF &
         ( OpacityTable, TRIM(WriteTableNameBrem), WriteOpacity_Brem_Option = .true. )
    CALL FinalizeHDF( )
  END IF

  WRITE (*,*) "HDF write successful"

!  IF ( EmAb_nuclei_EC_table .gt. 0) THEN
!    DO i_r = 1, OpacityTable % EmAb % EC_table_nOpacities
!      DEALLOCATE( OpacityTable % EmAb % EC_table_spec(i_r) % Values )
!      DEALLOCATE( OpacityTable % EmAb % EC_table_rate(i_r) % Values )
!    ENDDO
!    DEALLOCATE(OpacityTable % EmAb % EC_table_spec)
!    DEALLOCATE(OpacityTable % EmAb % EC_table_rate)
!  ENDIF

  CALL DeAllocateOpacityTable( OpacityTable )
    
  !=============================================================

  CALL itime(now)
  WRITE ( *, 10000 )  today(2), today(1), today(3), now

!warn the user if they have uncommited changes to the code when building tables.
#if defined(GIT_HASH)
        IF(index(GIT_HASH,"dirty") .gt. 0) THEN
          write(*,*) 'WARNING: There are uncommited changes in the repo; &
                      consider commiting your latest changes if you are &
                      building production opacity tables in order to store &
                      a clean git hash allowing you to later see the exact &
                      version of the code the tables were built with..'
        ENDIF
#endif

END PROGRAM wlCreateOpacityTable
