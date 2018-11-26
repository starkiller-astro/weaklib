PROGRAM wlCreateOpacityTable
!---------------------------------------------------------------------
!
!    Author:       R. Chu, Dept. Phys. & Astronomy
!                  U. Tennesee, Knoxville
!
!    Created:      10/23/18
!    WeakLib ver:  
!
!    Purpose:
!      Create table for opacity containing EoS table with BoltzTran.
!      Existing EoS table is readin. 
!      Function were created to fill OpacityTable and Chimera
!      routines were called.
!!   
!    CONTAINS:
!    
!
!    Modules used:
!
!---------------------------------------------------------------------
!  NOTE: Only Type A interaction applied. Type B and Type C 
!        interaction needs to be added for future use.
!---------------------------------------------------------------------
!                         Three Opacity Type
!
! OpacityType A for  ABEM( rho, T, Ye, E)
!
!       e- + p/A <--> v_e + n/A*
!
! OpacityType B for  ISO( e, rho, T, Ye, l)
!                        
!       v_i/anti(v_i) + A --> v_i/anti(v_i) + A 
!       v_i/anti(v_i) + e+/e-/n/p  <-->  v_i/anti(v_i) + e+/e-/n/p
!
! OpacityType C for  NISO( e_in, e_out, rho, T, Ye, l)
!
!       e+ + e-  <--> v_i + anti(v_i);   i=e, muon, tau
!       N + N   <--> N + N + v_i + anti(v_i)
!                 
!---------------------------------------------------------------------


  USE wlKindModule, ONLY: dp
  USE wlGridModule, ONLY: MakeLogGrid 
  USE wlIOModuleHDF, ONLY: &
      InitializeHDF,       &
      FinalizeHDF
  USE wlOpacityTableModule, ONLY: &
      OpacityTableType,     &
      AllocateOpacityTable, &
      DescribeOpacityTable
  USE wlOpacityTableIOModuleHDF, ONLY: &
      WriteOpacityTableHDF
  USE wlExtPhysicalConstantsModule, ONLY: kMeV
  USE wlExtNumericalModule, ONLY: epsilon
  USE B85_scattIso
  USE B85_scattNES
  USE B85_pair
  USE prb_cntl_module, ONLY: &
      i_aeps, iaefnp, rhoaefnp, iaence, iaenct, roaenct, &
      edmpa, edmpe, iaenca

IMPLICIT NONE

   INTEGER*4 today(3), now(3)    
   TYPE(OpacityTableType)  :: OpacityTable
   INTEGER                 :: nOpacA = 0
   INTEGER                 :: nOpacB = 0
   INTEGER                 :: nMomB  = 0
   INTEGER                 :: nOpacB_NES = 2
   INTEGER                 :: nMomB_NES  = 1
   INTEGER                 :: nOpacB_TP  = 1
   INTEGER                 :: nMomB_TP   = 1
   INTEGER                 :: nOpacC = 0
   INTEGER                 :: nMomC  = 0
!---------------------------------------------------------------------
! Set E grid limits
!---------------------------------------------------------------------
   INTEGER,  PARAMETER     :: nPointsE = 40 
   REAL(dp), PARAMETER     :: Emin = 1.0d-1  
   REAL(dp), PARAMETER     :: Emax = 3.0d02

!---------------------------------------------------------------------
! Set Eta grid limits
!---------------------------------------------------------------------
   INTEGER                 :: nPointsEta = 20!100
   REAL(dp), PARAMETER     :: Etamin = 1.0d-8
   REAL(dp), PARAMETER     :: Etamax = 2.5d03

   ! --- other inner variables
   INTEGER                 :: i_r, i_rb, i_e, j_rho, k_t, l_ye, &
                              t_m, i_quad, i_eta, i_ep
   REAL(dp)                :: energy, rho, T, TMeV, ye, Z, A, &
                              chem_e, chem_n, chem_p, xheavy, xn, &
                              xp, xhe, bb, eta, minvar 
   REAL(dp), DIMENSION(nPointsE) :: absor, emit, cok 
   INTEGER, PARAMETER      :: nquad = 30
   REAL(dp), DIMENSION(nPointsE, nPointsE) :: NESK
   REAL(dp), DIMENSION(nPointsE, nPointsE) :: TPK
   REAL(dp)                :: j0i, j0ii, j1i, j1ii, bufferTP

   CALL idate(today)
   CALL itime(now)
   WRITE ( *, 10000 )  today(2), today(1), today(3), now

   10000 format ( ' Date ', i2.2, '/', i2.2, '/', i4.4, '; Time ',&
     &         i2.2, ':', i2.2, ':', i2.2 )

   CALL InitializeHDF( ) 
   CALL AllocateOpacityTable &
            ( OpacityTable, nOpacA, nOpacB, nMomB, &
              nOpacB_NES, nMomB_NES, nOpacB_TP, nMomB_TP, &
              nOpacC, nMomC, nPointsE, nPointsEta ) 
   CALL FinalizeHDF( )

! -- Set OpacityTableTypeA thermEmAb  

   OpacityTable % thermEmAb % nOpacities = nOpacA

   OpacityTable % thermEmAb % nPoints(1) = nPointsE

   OpacityTable % thermEmAb % nPoints(2:4) = OpacityTable % nPointsTS

   OpacityTable % thermEmAb % Names = &
                                (/'Electron Capture  Chimera    ', &
                                  'AntiElec Capture  Chimera    '/)  

   OpacityTable % thermEmAb % Species = &
                                (/'Electron Neutrino           ', &  
                                  'Electron Antineutrino       '/)

   OpacityTable % thermEmAb % Units = &
                                (/'Per Centimeter              ', &
                                  'Per Centimeter              '/) 

! -- Set OpacityTableTypeB scatt_Iso

   OpacityTable % scatt_Iso % nOpacities = nOpacB
   
   OpacityTable % scatt_Iso % nMoments = nMomB
   
   OpacityTable % scatt_Iso % nPoints(1) = nPointsE

   OpacityTable % scatt_Iso % nPoints(2:4) = OpacityTable % nPointsTS

   OpacityTable % scatt_Iso % Names = &
                                (/'Electron Iso-scattering Chime', &
                                  'AntiElec Iso-scattering Chime'/)

   OpacityTable % scatt_Iso % Species = &
                                (/'Electron Neutrino           ', &
                                  'Electron AntiNeutrino       '/)

   OpacityTable % scatt_Iso % Units = &
                                (/'Per Centimeter              ', &
                                  'Per Centimeter              '/)

! -- Set OpacityTableTypeB scatt_NES

   OpacityTable % scatt_NES % nOpacities = nOpacB_NES

   OpacityTable % scatt_NES % nMoments   = nMomB_NES

   OpacityTable % scatt_NES % nPoints(1) = nPointsE
   OpacityTable % scatt_NES % nPoints(2) = nPointsE
   OpacityTable % scatt_NES % nPoints(3) = OpacityTable % nPointsTS(2)
   OpacityTable % scatt_NES % nPoints(4) = nPointsEta

   OpacityTable % scatt_NES % Names = &
                                (/'Electron NES Kernel Moment Ch', &
                                  'Antielec NES Kernel Moment Ch'/)

   OpacityTable % scatt_NES % Species = &
                                (/'Electron Neutrino           ', &
                                  'Electron AntiNeutrino       '/)

   OpacityTable % scatt_NES % Units = &
                                (/'Per Centimeter Per MeV^3    ', &
                                  'Per Centimeter Per MeV^3    '/)

! -- Set OpacityTableTypeB scatt_TP

   OpacityTable % scatt_TP % nOpacities = nOpacB_TP

   OpacityTable % scatt_TP % nMoments = nMomB_TP

   OpacityTable % scatt_TP % nPoints(1) = nPointsE
   OpacityTable % scatt_TP % nPoints(2) = nPointsE
   OpacityTable % scatt_TP % nPoints(3) = OpacityTable % nPointsTS(2)
   OpacityTable % scatt_TP % nPoints(4) = nPointsEta

   OpacityTable % scatt_TP % Names = &
                                (/'Pair P. Kernel Momen Chimera'/)

   OpacityTable % scatt_TP % Species = &
                                (/'Electron Neutrino-Antineutri'/)

   OpacityTable % scatt_TP % Units = &
                                (/'Per Centimeter Per MeV^3    '/)

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
       iRho    => OpacityTable % EOSTable % TS % Indices % iRho , &
       iT      => OpacityTable % EOSTable % TS % Indices % iT   , &
       iYe     => OpacityTable % EOSTable % TS % Indices % iYe  , &
       Indices => OpacityTable % EOSTable % DV % Indices        , &
       DVOffs  => OpacityTable % EOSTable % DV % Offsets        , &
       DVar    => OpacityTable % EOSTable % DV % Variables  )

!-----------------  ECAPEM -------------------- 
   PRINT*, 'Calculating thermEmAb and Elastic Scattering Kernel ...'

   DO l_ye = 1, OpacityTable % nPointsTS(iYe)
   
     ye = OpacityTable % EOSTable % TS % States (iYe) % Values (l_ye)

     DO k_t = 1, OpacityTable % nPointsTS(iT)

       T = OpacityTable % EOSTable % TS % States (iT) % Values (k_t)

       DO j_rho = 1, OpacityTable % nPointsTS(iRho)

         rho = OpacityTable % EOSTable % TS % States (iRho) % &
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
       
         DO i_r = 1, nOpacA
          
             iaefnp = 1
!  iaefnp   : emission and absorption of e-neutrinos and e-antineutrinos
!   on free neutrons and protons switch.
!
!     iaefnp = 0: emission and absorption of e-neutrinos and
!      e-antineutrinos on free neutrons and protons turned off.
!     iaefnp = 1: emission and absorption of e-neutrinos and
!      e-antineutrinos on free neutrons and protons turned on.
             i_aeps = 0
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
             rhoaefnp = HUGE(1.d0) ! (?) 
!  rhoaefnp : density above which emission and absorption of e-neutrinos
!   and e-antineutrinos on free neutrons and protons is turned off.
             iaence = 1
!  iaence   : emission and absorption of e-neutrinos on nuclei as given
!   by the FFN prescription  switch
!
!     iaence = 0: emission and absorption of e-neutrinos on nuclei as
!      given by the FFN prescription turned off.
!     iaence = 1: emission and absorption of e-neutrinos on nuclei as
!      given by the FFN prescription turned on.
             edmpe = 3.d0
!  edmpe    : difference between the mean excited state energy of the
!   daughter nucleus and that of the parent nucleus [MeV].
             iaenca = 1 
!  iaenca   : emission and absorption of e-antineutrinos on nuclei as
!   given by the FFN prescription  switch
!
!     iaenca = 0: emission and absorption of e-antineutrinos on nuclei as
!      given by the FFN prescription turned off.
!     iaenca = 1: emission and absorption of e-antineutrinos on nuclei as
!      given by the FFN prescription turned on.
             edmpa = 3.d0
!  edmpa    : difference between the mean excited state energy of the
!   daughter nucleus and that of  the parent nucleus [MeV].
             iaenct = 0
!  iaenct: emission and absorption of e-neutrinos on nuclei as given by
!   the NSE-folded tabular data for Hix et al (2003)
!
!     iaenct  = 0: emission and absorption of e-neutrinos on nuclei as
!      given by the tabular data is turned off.
!     iaenct  = 1: emission and absorption of e-neutrinos on nuclei as
!      given by the tabular data is turned on.
             roaenct = TINY(1.d0)
!  roaenct: density above which emission and absorption of e-neutrinos on
!   nuclei as given by Hix et al (2003) is turned off.

             CALL abemrgn_weaklib &
                  ( i_r, OpacityTable % EnergyGrid % Values, &
                    rho, T, xn, xp, xheavy, &
                    A, Z, chem_n, chem_p, chem_e, & 
                    absor, emit, &
                    ye, nPointsE )

             DO i_e = 1, OpacityTable % nPointsE
                OpacityTable % thermEmAb % Absorptivity(i_r) % &
                        Values (i_e, j_rho, k_t, l_ye) &
                = absor(i_e) + emit(i_e)
             END DO  !i_e

          ! i_r = 1 = iNue for Eletron-type neutrino Chimera 
          ! i_r = 2 = iNuebar for Eletron-type antineutrino Chimera 
          
         END DO !i_r

!----------------  Scatt_Iso -----------------------
         DO i_rb = 1, nOpacB
           DO t_m = 1, nMomB

             CALL scatical_weaklib &
             ( i_rb, t_m-1, OpacityTable % EnergyGrid % Values, &
               nPointsE, rho, T, xn, xp, xhe, xheavy, A, Z, cok )

             OpacityTable % scatt_Iso % Kernel(i_rb) % Values &
             ( :, j_rho, k_t, l_ye, t_m )  = cok 

           END DO !t_m         
         END DO !i_rb

       END DO  !j_rho
     END DO  !k_t
   END DO  !l_ye

!------- thermEmAb % Offsets
   DO i_r = 1, nOpacA  
     minvar = MINVAL( OpacityTable % thermEmAb % Absorptivity(i_r) % Values )
     OpacityTable % thermEmAb % Offsets(i_r) = -2.d0 * MIN( 0.d0, minvar ) 
   END DO

!------- scatt_Iso % Offsets
   DO i_r = 1, nOpacB
     DO t_m = 1, nMomB
       minvar = MINVAL( OpacityTable % scatt_Iso % Kernel(i_r) &
                       % Values(:,:,:,:,t_m ) )
       OpacityTable % scatt_Iso % Offsets(i_r, t_m) =           &
                      -2.d0 * MIN( 0.d0, minvar )
     END DO
   END DO 

!----------------  Scatt_NES -----------------------

  PRINT*, 'Calculating Scatt_NES Kernel ... '

  DO i_rb = 1, nOpacB_NES
    DO t_m = 1, nMomB_NES 
      DO i_eta = 1, nPointsEta
       
        eta = OpacityTable % EtaGrid % Values(i_eta)
      
        DO k_t = 1, OpacityTable % nPointsTS(iT)

          T = OpacityTable % EOSTable % TS % States (iT) % Values (k_t)
          TMeV = T * kMeV
          chem_e = TMeV * eta

          CALL scatergn_weaklib&
              ( i_rb, t_m-1, nPointsE, OpacityTable % EnergyGrid % Values, TMeV, chem_e, NESK )

          DO i_e = 1, nPointsE
          
            DO i_ep = 1, nPointsE
             
              OpacityTable % scatt_NES % Kernel(i_rb) % Values &
                           ( i_ep, i_e, k_t, i_eta, t_m)       &
               = NESK(i_e, i_ep) ! NESK was saved as NESK(e,ep)
 
            END DO ! i_ep
          END DO ! i_e
        END DO  !k_t
      END DO !i_eta

!------- scatt_NES % Offsets
       minvar = MINVAL( OpacityTable % scatt_NES % Kernel(i_rb) &
                       % Values(:,:,:,:,t_m ) )
       OpacityTable % scatt_NES % Offsets(i_rb, t_m) =           &
                      -2.d0 * MIN( 0.d0, minvar ) 

    END DO ! t_m
  END DO ! i_rb
    
!----------------  Scatt_TP -----------------------

  PRINT*, 'Calculating Scatt_TP Kernel ... '

  DO i_rb = 1, nOpacB_TP 
    DO t_m = 1, nMomB_TP 
      DO i_eta = 1, nPointsEta
   
        eta = OpacityTable % EtaGrid % Values(i_eta)
      
        DO k_t = 1, OpacityTable % nPointsTS(iT)

          T = OpacityTable % EOSTable % TS % States (iT) % Values (k_t)
          TMeV = T * kMeV
          chem_e = TMeV * eta

          DO i_e = 1, nPointsE
          
            DO i_ep = 1, nPointsE
             
             CALL paircal_weaklib( OpacityTable % EnergyGrid % Values(i_e), &
                           OpacityTable % EnergyGrid % Values(i_ep), &
                           chem_e, T, j0i, j0ii, j1i, j1ii )

             bufferTP = j0i + j0ii

             OpacityTable % scatt_TP % Kernel(i_rb) % Values &
                          ( i_ep, i_e, k_t, i_eta, t_m)       &
              = bufferTP
 
            END DO ! i_ep
          END DO ! i_e
        END DO  !k_t
      END DO !i_eta

!------- scatt_TP % Offsets
       minvar = MINVAL( OpacityTable % scatt_TP % Kernel(i_rb) &
                       % Values(:,:,:,:,t_m ) )
       OpacityTable % scatt_TP % Offsets(i_rb, t_m) =          &
                      -2.d0 * MIN( 0.d0, minvar ) 

    END DO ! t_m

  END DO ! i_rb

   END ASSOCIATE ! rho-T-Ye
!---------------------------------------------------------------------
!      Describe the Table ( give the real physical value )            
!---------------------------------------------------------------------

  CALL DescribeOpacityTable( OpacityTable )

! --------------------------------------------------------------------
!          LOG the WHOLE table for storage
! --------------------------------------------------------------------

  WRITE(*,*) 'LOG the whole table with relevant offset for storage'

  DO i_r = 1, nOpacA
     OpacityTable % thermEmAb % Absorptivity(i_r) % Values&
     = LOG10( OpacityTable % thermEmAb % Absorptivity(i_r) % &
              Values + OpacityTable % thermEmAb % Offsets(i_r) + epsilon )
  END DO  !i_r

  DO i_rb = 1, nOpacB
    DO t_m = 1, nMomB
      OpacityTable % scatt_Iso % Kernel(i_rb) % Values(:,:,:,:,t_m) &
      = LOG10 ( OpacityTable % scatt_Iso % Kernel(i_rb) % Values(:,:,:,:,t_m) &
                + OpacityTable % scatt_Iso % Offsets(i_rb,t_m) + epsilon )
    END DO ! t_m
  END DO ! i_rb

  DO i_rb = 1, nOpacB_NES
    DO t_m = 1, nMomB_NES
      OpacityTable % scatt_NES % Kernel(i_rb) % Values(:,:,:,:,t_m) &
      = LOG10 ( OpacityTable % scatt_NES % Kernel(i_rb) % Values(:,:,:,:,t_m) &
                + OpacityTable % scatt_NES % Offsets(i_rb, t_m) + epsilon )
    END DO
  END DO !i_rb

  DO i_rb = 1, nOpacB_TP
    DO t_m = 1, nMomB_TP
      OpacityTable % scatt_TP % Kernel(i_rb) % Values(:,:,:,:,t_m) &
      = LOG10 ( OpacityTable % scatt_TP % Kernel(i_rb) % Values(:,:,:,:,t_m) &
                + OpacityTable % scatt_TP % Offsets(i_rb, t_m) + epsilon )
    END DO
  END DO !i_rb

! -- write into hdf5 file

  CALL InitializeHDF( )
  WRITE(*,*) 'Write data into file temp.h5 '
  CALL WriteOpacityTableHDF( OpacityTable, "temp.h5" )
  CALL FinalizeHDF( )
  
  WRITE (*,*) "HDF write successful"

  !=============================================================

  CALL itime(now)
  WRITE ( *, 10000 )  today(2), today(1), today(3), now

END PROGRAM wlCreateOpacityTable
