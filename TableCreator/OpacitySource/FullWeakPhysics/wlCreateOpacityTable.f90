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
! Stored is the zeroth order annihilation kernel from Hannestad and Raffelt 1998
! for a generic rho. The composition can later be taken into account
! by interpolating and adding the contribution from
! xp*rho, xn*rho and sqrt(xp+xn)*rho
!
!       nu nubar NN <--> NN   i=e, muon, tau
!---------------------------------------------------------------------


  USE wlKindModule, ONLY: dp
  USE wlGridModule, ONLY: MakeLogGrid
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
      edmpa, edmpe, iaenca

  USE abem_module_weaklib

  USE CC_module_weaklib

  USE ec_table_module

  USE, INTRINSIC :: iso_fortran_env, only : stdin=>input_unit, &
                                            stdout=>output_unit, &
                                            stderr=>error_unit
  USE omp_lib

IMPLICIT NONE

   INTEGER*4 today(3), now(3)

!---------------------------------------------------------------------
! Set Eos table
!---------------------------------------------------------------------
   CHARACTER(256) :: EOSTableName = "wl-EOS-SFHo-15-25-50.h5"
   !CHARACTER(256) :: EOSTableName = "wl-EOS-SFHo-25-50-100.h5"
   !CHARACTER(256) :: EOSTableName = "wl-EOS-LS220-25-50-100.h5"

!---------------------------------------------------------------------
! Set neutrino interaction type
!---------------------------------------------------------------------
   CHARACTER(256)          :: WriteTableName
   CHARACTER(256)          :: WriteTableNameBrem

   TYPE(OpacityTableType)  :: OpacityTable
   INTEGER                 :: nOpac_EmAb = 2  ! 2 for electron type

   INTEGER                 :: EmAb_nucleons_full_kinematics = 0 !Use Fischer et al 2020 full kinematics rates for
                                                                !EmAb on free nucleons

   INTEGER                 :: inv_n_decay_full_kinematics   = 0 !Use Fischer et al 2020 full kinematics rates for
                                                                !inverse neutron decay 

   INTEGER                 :: EmAb_nucleons_isoenergetic    = 1 !EmAb on free nucleons using isoenergetic approximation
                                                                !Bruenn 1985
                                                                !Mezzacappa & Bruenn (1993)

   INTEGER                 :: EmAb_nucleons_recoil          = 1 !EmAb on free nucleons taking into account recoil,
                                                                !nucleon final-state blocking, and special relativity
                                                                !Reddy et al 1998
                                                                !Only used for rho > 1e9

   INTEGER                 :: EmAb_nucleons_weak_magnetism  = 1 !Weak magnetism corrections for EmAb on free nucleons
                                                                !Horowitz 1997

   INTEGER                 :: EmAb_nuclei_FFN               = 0 !EmAb on nuclei using FFN formalism
                                                                !Fuller, Fowler, Neuman 1982, Ap. J. 252, 715
                                                                !Bruenn 1985

   INTEGER                 :: EmAb_nuclei_Hix               = 1 !EmAb on nuclei using NSE-folded tabular data
                                                                !Langanke et al. (2003), Hix et al. (2003)

   INTEGER                 :: nOpac_Iso  = 2  ! 2 for electron type
                                              !   ( flavor identical )
   INTEGER                 :: nMom_Iso   = 2  ! 2 for 0th & 1st order
                                              !   legendre coff.

   INTEGER                 :: nOpac_NES  = 1  ! 1 ( either 0 or 1 )
   INTEGER                 :: nMom_NES   = 4  ! 4 for H1l, H2l
                                              !   ( either 0 or 4 )

   INTEGER                 :: nOpac_Pair = 1  ! 1 ( either 0 or 1 )
   INTEGER                 :: nMom_Pair  = 4  ! 4 for J1l, J2l
                                              !   ( either 0 or 4 )

   INTEGER                 :: nOpac_Brem = 1  !1 just the zeroth order annihilation kernel
   INTEGER                 :: nMom_Brem  = 1  !1 only need to calculate zeroth order kernel

!---------------------------------------------------------------------
! Set E grid limits
!---------------------------------------------------------------------
   INTEGER,  PARAMETER     :: nPointsE = 40
   REAL(dp), PARAMETER     :: Emin = 1.0d-1
   REAL(dp), PARAMETER     :: Emax = 3.0d02

!---------------------------------------------------------------------
! Set Eta grid limits
!---------------------------------------------------------------------
   INTEGER                 :: nPointsEta = 120
   REAL(dp), PARAMETER     :: Etamin = 1.0d-3
   REAL(dp), PARAMETER     :: Etamax = 2.5d03

   ! --- other inner variables
   INTEGER                 :: i_r, i_rb, i_e, j_rho, k_t, l_ye, &
                              t_m, i_eta, i_ep, stringlength, stringlengthBrem
   REAL(dp)                :: rho, T, TMeV, ye, Z, A, &
                              chem_e, chem_n, chem_p, xheavy, xn, &
                              xp, xhe, bb, eta, minvar

   REAL(dp), DIMENSION(nPointsE,2) :: cok
   REAL(dp), DIMENSION(nPointsE, nPointsE) :: H0i, H0ii, H1i, H1ii
   REAL(dp)                                :: j0i, j0ii, j1i, j1ii

   REAL(dp), DIMENSION(nPointsE, nPointsE) :: s_a !the Bremsstrahlung annihilation kernel
   REAL(dp), PARAMETER                     :: brem_rho_min = 1.0d+07 !switch Bremsstrahlung off below rho_min
   REAL(dp), PARAMETER                     :: brem_rho_max = 1.0d+15 !switch Bremsstrahlung off above rho_max

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
   WRITE(WriteTableName,'(A,A,A,A,I2,A)') &
         EOSTableName(1:3),'Op-',EOSTableName(8:stringlength-3),'-E',nPointsE, &
         '-B85'
   stringlength = len(TRIM(WriteTableName))

   stringlengthBrem = LEN(TRIM(EOSTableName))
   WRITE(WriteTableNameBrem,'(A,A,A,A,I2,A)') &
         EOSTableName(1:3),'Op-',EOSTableName(8:stringlengthBrem-3),'-E',nPointsE, &
         '-HR98'
   stringlengthBrem = len(TRIM(WriteTableNameBrem))

   CALL InitializeHDF( )
   CALL AllocateOpacityTable &
            ( OpacityTable, nOpac_EmAb, nOpac_Iso, nMom_Iso, &
              nOpac_NES, nMom_NES, nOpac_Pair, nMom_Pair, &
              nOpac_Brem, nMom_Brem, &
              nPointsE, nPointsEta, &
              EquationOfStateTableName_Option = EOSTableName )
   CALL FinalizeHDF( )

! -- Set OpacityTableTypeEmAb

   IF( nOpac_EmAb .gt. 0 ) THEN

   if(EmAb_nucleons_full_kinematics + EmAb_nucleons_isoenergetic .gt. 1) then
     write (*,*) "Must choose either full kinematics or isoenergetic approach for EmAb on free nucleons"
     return
   endif
   if(EmAb_nuclei_FFN + EmAb_nuclei_Hix .gt. 1) then
     write (*,*) "Must choose either FFN formalism or table for EmAb on nuclei"
     return
   endif

   OpacityTable % EmAb % nOpacities   = nOpac_EmAb

   OpacityTable % EmAb % nPoints(1)   = nPointsE

   OpacityTable % EmAb % nPoints(2:4) = OpacityTable % nPointsTS

   OpacityTable % EmAb % Names = &
                                (/'Electron Neutrino           ',  &
                                  'Electron Antineutrino       '/)

   OpacityTable % EmAb % Units = &
                                (/'Per Centimeter              ',  &
                                  'Per Centimeter              '/)

   OpacityTable % EmAb % nucleons_full_kinematics    = EmAb_nucleons_full_kinematics
   OpacityTable % EmAb % inv_n_decay_full_kinematics = inv_n_decay_full_kinematics
   OpacityTable % EmAb % nucleons_isoenergetic       = EmAb_nucleons_isoenergetic
   OpacityTable % EmAb % nucleons_recoil             = EmAb_nucleons_recoil
   OpacityTable % EmAb % nucleons_weak_magnetism     = EmAb_nucleons_weak_magnetism
   OpacityTable % EmAb % nuclei_FFN                  = EmAb_nuclei_FFN
   OpacityTable % EmAb % nuclei_Hix                  = EmAb_nuclei_Hix

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

   OpacityTable % Scat_Brem % Names      = (/'Kernels'/)

   OpacityTable % Scat_Brem % Units      = (/'Per Centimeter Per MeV^3'/)

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

   CALL MakeLogGrid &
          ( EnergyGrid % MinValue, EnergyGrid % MaxValue, &
            EnergyGrid % nPoints,  EnergyGrid % Values )

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
       nu_E    => OpacityTable % EnergyGrid % Values                           , &
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

    REAL(dp), DIMENSION(nPointsE) :: em_nucleons,    ab_nucleons
    REAL(dp), DIMENSION(nPointsE) :: em_inv_n_decay, ab_inv_n_decay
    REAL(dp), DIMENSION(nPointsE) :: em_nuclei_FFN,  ab_nuclei_FFN
    REAL(dp), DIMENSION(nPointsE) :: em_nuclei_Hix,  ab_nuclei_Hix

    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: EOS_Un, EOS_Up !neutron and proton mean field potentials from nuclear EOS
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: EOS_massn, EOS_massp !effective neutron and proton masses from nuclear EOS

    REAL(dp) :: Un_loc, Up_loc
    REAL(dp) :: massn_loc, massp_loc

    !arrays for weak magnetism corrections
    REAL(dp), DIMENSION(nPointsE) :: xi_n_wm, xib_p_wm

    REAL(dp), DIMENSION(nPointsE+1) :: nu_E_edge
    REAL(dp), DIMENSION(nPointsE)   :: dnu_E_edge

    IF(EmAb_nucleons_full_kinematics .gt. 0 .or. &
       inv_n_decay_full_kinematics   .gt. 0) then

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
   IF ( EmAb_nuclei_Hix .gt. 0) THEN

     prepare_EC_table: BLOCK

       CHARACTER(len=255) :: ec_table
       integer            :: k     

       CALL get_environment_variable('WEAKLIB_DIR', ec_table)

       IF( LEN_TRIM(ec_table) .le. 0) THEN

         WRITE(*,*) 'Environment variable $WEAKLIB_DIR not set, please setup your weaklib environment!'
 
         STOP

       ENDIF

       ec_table = TRIM(ec_table) // TRIM("/TableCreator/OpacitySource/FullWeakPhysics/Executables/")

       ec_table = TRIM(ec_table) // TRIM("ec_table.d")

       CALL read_ec_table(ec_table)

       WRITE(*,*) 'Read in EC table from file ', ec_table 

       ASSOCIATE(Egrid => OpacityTable % EnergyGrid)

       IF( EGrid % LogInterp == 1) THEN

         create_unubi_dunui: BLOCK
           REAL(dp) :: runul, runu, runuh

           runul = DLOG(EGrid % MaxValue / EGrid % MinValue) &
                   / DBLE(EGrid % nPoints - 1) 
           runu  = DEXP(runul)
           runuh = DSQRT(runu) 

           nu_E_edge(1) = EGrid % MinValue / runuh

           DO k=2,nPointsE+1
   
             nu_E_edge(k) = runu * nu_E_edge(k-1)

           ENDDO

           dnu_E_edge(1:nPointsE) = nu_E_edge(2:nPointsE+1) - nu_E_edge(1:nPointsE)

         END BLOCK create_unubi_dunui

       ELSE

         WRITE(*,*) 'Needs to be implemented for non-log grid still!'
         STOP 

       ENDIF

       END ASSOCIATE

     END BLOCK prepare_EC_table

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

            bb  = (chem_e + chem_p - chem_n)/(T*kMev)

            if (     EmAb_nucleons_full_kinematics .gt. 0 &
                .or. inv_n_decay_full_kinematics   .gt. 0 ) then
 
              Un_loc = EOS_Un(j_rho, k_t, l_ye)
              Up_loc = EOS_Up(j_rho, k_t, l_ye)

              massn_loc = EOS_massn(j_rho, k_t, l_ye)
              massp_loc = EOS_massp(j_rho, k_t, l_ye)

            endif

         DO i_r = 1, nOpac_EmAb

             !iaefnp = 1
             iaefnp = EmAb_nucleons_isoenergetic
             !i_aeps = 0
             i_aeps = EmAb_nucleons_recoil
             rhoaefnp = HUGE(1.d0) ! (?)
             !iaence = 1
             iaence = EmAb_nuclei_FFN !e-nu
             edmpe = 3.d0
             !iaenca = 1
             iaenca = EmAb_nuclei_FFN !e-nubar
             edmpa = 3.d0
             !iaenct = 0
             iaenct = EmAb_nuclei_Hix
             roaenct = 1.0d+13 !TINY(1.d0)

! Initialise the opacity arrays to zero
             ab_nucleons = 0.0d0
             em_nucleons = 0.0d0

             ab_inv_n_decay = 0.0d0
             em_inv_n_decay = 0.0d0

             ab_nuclei_FFN = 0.0d0
             em_nuclei_FFN = 0.0d0

             ab_nuclei_Hix = 0.0d0
             em_nuclei_Hix = 0.0d0

             IF(EmAb_nucleons_isoenergetic .gt. 0) THEN

               CALL abem_nucleons_isoenergetic_weaklib & 
                    ( i_r, nu_E, &
                      rho, T, ye, xn, xp, xheavy, A, Z, chem_n, chem_p, chem_e, &
                      ab_nucleons, em_nucleons, nPointsE)

               !recoil corrections to EmAb on nucleons
               IF(EmAb_nucleons_recoil .gt. 0 .and. rho .ge. 1.d+09) THEN

                 CALL abem_nucleons_recoil_weaklib & 
                      ( i_r, nu_E, &
                        rho, T, ye, xn, xp, xheavy, A, Z, chem_n, chem_p, chem_e, &
                        ab_nucleons, em_nucleons, nPointsE)

               ENDIF !EmAb_nucleons_recoil .gt. 0 .and. rho .gt. 1.d+09

               !weak magnetism corrections to EmAb on nucleons
               IF(EmAb_nucleons_weak_magnetism .gt. 0) THEN

                 CALL cc_weak_mag_weaklib( nu_E, xi_n_wm, xib_p_wm, nPointsE )

                 IF (i_r .eq. 1) THEN !electron neutrino weak mag correction

                   ab_nucleons = ab_nucleons * xi_n_wm
                   em_nucleons = em_nucleons * xi_n_wm

                 ELSE IF (i_r .eq. 2) THEN !electron anti-neutrino weak mag correction

                   ab_nucleons = ab_nucleons * xib_p_wm
                   em_nucleons = em_nucleons * xib_p_wm

                 ENDIF !i_r .eq. 1 or 2

               ENDIF !EmAb_nucleons_weak_magnetism .gt. 0

             ENDIF !EmAb_nucleons_isoenergetic .gt. 0

             IF(EmAb_nucleons_full_kinematics .gt. 0) THEN

!write(*,*) 'i=', i_r, 'j=', j_rho, 'k=', k_t, 'l=', l_ye
!write(*,*) 'T =',TMeV,'Density =',rho,'Ye =',Ye
!write(*,*) 'n =',xn,'p =',xp,'heavy =',xheavy,'A =',A,'Z =',Z
!write(*,*) 'chemn =',chem_n,'chemp =',chem_p,'Un =',Un_loc,'Up =',Up_loc
!write(*,*) 'mn =',massn_loc,'mp =',massp_loc
!write(*,*) 'cheme =',chem_e

!stop
               IF(i_r .le. 2) THEN
                 CALL CC_EmAb( OpacityTable % EnergyGrid % Values,     &
                               TMeV, mass_e, chem_e, chem_n+mass_n, chem_p+mass_p,      &
                               massn_loc, massp_loc, Un_loc, Up_loc,   &
                               i_r, inv_n_decay_full_kinematics,       &
                               ab_nucleons, em_nucleons,               &
                               ab_inv_n_decay, em_inv_n_decay, nPointsE)
                 !CALL CC_EmAb( nu_E,     &
                 !              0.1d0, mass_e, 0.3148d-03, -0.1656d+01+mass_n, -0.2020d+02+mass_p,      &
                 !              0.9396d+03, 0.9383d+03, 0.1935d-08, -0.1317d-03,   &
                 !              i_r, inv_n_decay_full_kinematics,       &
                 !              ab_nucleons, em_nucleons,               &
                 !              ab_inv_n_decay, em_inv_n_decay, nPointsE)

               ENDIF

             ENDIF

             IF(EmAb_nuclei_FFN .gt. 0 .and. i_r .eq. 1) THEN

               CALL abem_nuclei_FFN_weaklib( nu_E, rho, T, xheavy, A, Z, chem_n, chem_p, &
                                             chem_e, ab_nuclei_FFN, em_nuclei_FFN, &
                                             nPointsE )

             ENDIF

             IF(EmAb_nuclei_Hix .gt. 0 .and. i_r .eq. 1) THEN
  
               CALL abem_nuclei_EC_table_weaklib ( i_r, nu_E, nu_E_edge, dnu_E_edge, &
                                                   rho, T, Ye, xheavy, A, chem_n, &
                                                   chem_p, chem_e, ab_nuclei_Hix, em_nuclei_Hix, &
                                                   1, nPointsE)

             ENDIF

             DO i_e = 1, OpacityTable % nPointsE
                OpacityTable % EmAb % Opacity(i_r) % &
                        Values (i_e, j_rho, k_t, l_ye) &
                = ab_nucleons(i_e)    + em_nucleons(i_e) &
                + ab_inv_n_decay(i_e) + em_inv_n_decay(i_e) &
                + ab_nuclei_FFN(i_e)  + em_nuclei_FFN(i_e) &
                + ab_nuclei_Hix(i_e)  + em_nuclei_Hix(i_e)
             END DO  !i_e

         END DO !i_r

!----------------  Scat_Iso -----------------------
         DO i_rb = 1, nOpac_Iso

           CALL scatical_weaklib &
           ( i_rb, OpacityTable % EnergyGrid % Values, &
             nPointsE, rho, T, xn, xp, xhe, xheavy, A, Z, cok )

           DO t_m = 1, nMom_Iso

             OpacityTable % Scat_Iso % Kernel(i_rb) % Values &
             ( :, t_m, j_rho, k_t, l_ye )  = cok(:,t_m)

           END DO !t_m

         END DO !i_rb

       END DO  !j_rho
     END DO  !k_t
   END DO  !l_ye

!stop

!------- EmAb % Offsets
   DO i_r = 1, nOpac_EmAb
     minvar = MINVAL( OpacityTable % EmAb % Opacity(i_r) % Values )
     OpacityTable % EmAb % Offsets(i_r) = -2.d0 * MIN( 0.d0, minvar )
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
                 TMeV, eta, H0i, H0ii, H1i, H1ii )

          OpacityTable % Scat_NES % Kernel(1) % Values &
               ( :, :, 1, k_t, i_eta )       &
          = TRANSPOSE(H0i(:,:)) ! H0i was saved as H0i(e,ep)

          OpacityTable % Scat_NES % Kernel(1) % Values &
               ( :, :, 2, k_t, i_eta)       &
          = TRANSPOSE(H0ii(:,:)) ! H0ii was saved as H0ii(e,ep)

          OpacityTable % Scat_NES % Kernel(1) % Values &
               ( :, :, 3, k_t, i_eta )       &
          = TRANSPOSE(H1i(:,:)) ! H1i was saved as H1i(e,ep)
          OpacityTable % Scat_NES % Kernel(1) % Values &
               ( :, :, 4, k_t, i_eta )       &
          = TRANSPOSE(H1ii(:,:)) ! H1ii was saved as H1ii(e,ep)

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
              = j0i

             OpacityTable % Scat_Pair % Kernel(1) % Values  &
                          ( i_ep, i_e, 2, k_t, i_eta )       &
              = j0ii

             OpacityTable % Scat_Pair % Kernel(1) % Values  &
                          ( i_ep, i_e, 3, k_t, i_eta )       &
              = j1i

             OpacityTable % Scat_Pair % Kernel(1) % Values  &
                          ( i_ep, i_e, 4, k_t, i_eta )       &
              = j1ii

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

      IF (rho < brem_rho_min .or. rho > brem_rho_max) THEN

          OpacityTable % Scat_Brem % Kernel(1) % Values (:, :, 1, j_rho, k_t) &
          = 0.d0
          cycle

      ELSE

        CALL bremcal_weaklib &
             (nPointsE, OpacityTable % EnergyGrid % Values, &
              rho, T, s_a)

             OpacityTable % Scat_Brem % Kernel(1) % Values (:, :, 1, j_rho, k_t) &
             = TRANSPOSE(s_a(:,:))

      END IF

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

  CALL DeAllocateOpacityTable( OpacityTable )
  !=============================================================

  CALL itime(now)
  WRITE ( *, 10000 )  today(2), today(1), today(3), now

END PROGRAM wlCreateOpacityTable
