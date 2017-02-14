PROGRAM wlCreateOpacityTable
!-----------------------------------------------------------------------
!
!    File:         wlCreateOpacityTable.f90
!    Type:         Program 
!    Author:       R. Chu, Dept. Phys. & Astronomy
!                  U. Tennesee, Knoxville
!
!    Created:      2/29/16
!    WeakLib ver:  
!
!    Purpose:
!      Create table for opacity containing EoS table with BoltzTran.
!      Existing EoS table is readin. 
!      Function were created to fill OpacityTable
!!   
!    CONTAINS:
!    
!
!    Modules used:
!      wlKindModule, ONLY: dp
!      HDF5
!      wlGridModule, ONLY: MakeLogGrid
!      wlThermoStateModule
!      wlDependentVariablesModule
!      wlIOModuleHDF
!      wlOpacityTableModule
!      wlOpacityTableIOModuleHDF
!      wlExtPhysicalConstantsModule
!      B85
!
!-----------------------------------------------------------------------
!  NOTE: Only Type A interaction applied. Type B and Type C interaction 
!        needs to be added for future use.
!-----------------------------------------------------------------------
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
!-----------------------------------------------------------------------


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
  USE B85
  USE wlExtNumericalModule, ONLY: epsilon
  USE GreyVariables
 
IMPLICIT NONE

   INTEGER*4 today(3), now(3)    
   TYPE(OpacityTableType)  :: OpacityTable
   INTEGER                 :: nOpacA = 1
   INTEGER                 :: nOpacB = 1
   INTEGER                 :: nMomB  = 2
   INTEGER                 :: nOpacB_NES = 1
   INTEGER                 :: nMomB_NES  = 2
   INTEGER                 :: nOpacC = 0
   INTEGER                 :: nMomC  = 0
!-------------------------------------------------------------------------
! Set E grid limits
!-------------------------------------------------------------------------
   INTEGER,  PARAMETER     :: nPointsE = 40 
   REAL(dp), PARAMETER     :: Emin = 1.0d-1  
   REAL(dp), PARAMETER     :: Emax = 3.0d02

!-------------------------------------------------------------------------
! Set Eta grid limits
!-------------------------------------------------------------------------
   INTEGER                 :: nPointsEta = 100
   REAL(dp), PARAMETER     :: Etamin = 1.0d-8
   REAL(dp), PARAMETER     :: Etamax = 2.5d03

! --- inner variables
   INTEGER                 :: i_r, i_rb, i_e, j_rho, k_t, l_ye, t_m, i_quad, &
                              i_eta, i_ep
   REAL(dp)                :: energy, rho, T, TMeV, ye, Z, A, chem_e, chem_n, &
                              chem_p, xheavy, xn, xp, bb, eta, &
                              bufferquad1, bufferquad2, bufferquad3,&
                              bufferquad4, bufferquad21, bufferquad22, &
                              bufferquad23, bufferquad24
   INTEGER, PARAMETER      :: nquadGrey = 30
   INTEGER, PARAMETER      :: nquad = 30
   REAL(dp), DIMENSION(nPointsE, nPointsE) :: NESK

  CALL idate(today)
  CALL itime(now)
  write ( *, 1000 )  today(2), today(1), today(3), now
1000 format ( 'Date ', i2.2, '/', i2.2, '/', i4.4, '; time ',&
     &         i2.2, ':', i2.2, ':', i2.2 )

   CALL InitializeHDF( ) 
   CALL AllocateOpacityTable &
            ( OpacityTable, nOpacA, nOpacB, nMomB, &
              nOpacB_NES, nMomB_NES, nOpacC, nMomC, nPointsE, nPointsEta ) 
   CALL FinalizeHDF( )

! -- Set OpacityTableTypeA thermEmAb  

   OpacityTable % thermEmAb % nOpacities = nOpacA

   OpacityTable % thermEmAb % nPoints(1) = nPointsE

   OpacityTable % thermEmAb % nPoints(2:4) = OpacityTable % nPointsTS

   OpacityTable % thermEmAb % Names = &
                                (/'Electron Capture Absorptivity'/)  

   OpacityTable % thermEmAb % Species = &
                                (/'Electron Neutrino           '/)

   OpacityTable % thermEmAb % Units = &
                                (/'Per Centimeter              '/) 

   OpacityTable % thermEmAb % Offset = 1.0d-100

! -- Set OpacityTableTypeB scatt_Iso

   OpacityTable % scatt_Iso % nOpacities = nOpacB
   
   OpacityTable % scatt_Iso % nMoments = nMomB
   
   OpacityTable % scatt_Iso % nPoints(1) = nPointsE

   OpacityTable % scatt_Iso % nPoints(2:4) = OpacityTable % nPointsTS

   OpacityTable % scatt_Iso % Names = &
                                (/'Electron Iso-sca Kernel Momen'/)

   OpacityTable % scatt_Iso % Species = &
                                (/'Electron Neutrino           '/)

   OpacityTable % scatt_Iso % Units = &
                                (/'Per Centimeter              '/)

   OpacityTable % scatt_Iso % Offset = 1.0d-1

! -- Set OpacityTableTypeB scatt_NES

   OpacityTable % scatt_NES % nOpacities = nOpacB_NES

   OpacityTable % scatt_NES % nMoments = nMomB_NES

   OpacityTable % scatt_NES % nPoints(1) = nPointsE
   OpacityTable % scatt_NES % nPoints(2) = nPointsE
   OpacityTable % scatt_NES % nPoints(3) = OpacityTable % nPointsTS(2)
   OpacityTable % scatt_NES % nPoints(4) = nPointsEta

   OpacityTable % scatt_NES % Names = &
                                (/'Electron non-Iso Kernel Momen'/)

   OpacityTable % scatt_NES % Species = &
                                (/'Electron Neutrino           '/)

   OpacityTable % scatt_NES % Units = &
                                (/'Per Centimeter              '/)

   OpacityTable % scatt_NES % Offset = 1.0E-21 

!-----------------------------   
! Generate E grid from limits
!-----------------------------
PRINT*, "Making Energy Grid ... "

   ASSOCIATE( EnergyGrid => OpacityTable % EnergyGrid )

   EnergyGrid % Name &
     = 'Comoving Frame Neutrino Energy  '

   EnergyGrid % Unit &
     = 'MeV                             '

   EnergyGrid % MinValue = Emin
   EnergyGrid % MaxValue = Emax
   EnergyGrid % LogInterp = 1

   CALL MakeLogGrid &
          ( EnergyGrid % MinValue, EnergyGrid % MaxValue, &
            EnergyGrid % nPoints, EnergyGrid % Values )
 
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

!-------------------------------------------------------------------------
!              Fill OpacityTable
!-------------------------------------------------------------------------

PRINT*, 'Filling OpacityTable ...'
!-----------------  ECAPEM ----------------------- 
   ASSOCIATE( iRho    => OpacityTable % EOSTable % TS % Indices % iRho , &
              iT      => OpacityTable % EOSTable % TS % Indices % iT   , &
              iYe     => OpacityTable % EOSTable % TS % Indices % iYe )

PRINT*, 'Calculating thermEmAb and Elastic Scattering Kernel ...'
   DO l_ye = 1, OpacityTable % nPointsTS(3)
     
        ye = OpacityTable % EOSTable % TS % States (iYe) % Values (l_ye)

      DO k_t = 1, OpacityTable % nPointsTS(2)

           T = OpacityTable % EOSTable % TS % States (iT) % Values (k_t)

         DO j_rho = 1, OpacityTable % nPointsTS(1)

              rho = OpacityTable % EOSTable % TS % States (iRho) % Values (j_rho)

              chem_e = 10**OpacityTable % EOSTable % DV % Variables (4) %&
                       Values (j_rho, k_t, l_ye) - OpacityTable % EOSTable % &
                       DV % Offsets(4) - epsilon    !4 =Electron Chemical Potential             

              chem_p = 10**OpacityTable % EOSTable % DV % Variables (5) %&
                       Values (j_rho, k_t, l_ye) - OpacityTable % EOSTable % &
                       DV % Offsets(5) - epsilon    !5 =Proton Chemical Potential 

              chem_n = 10**OpacityTable % EOSTable % DV % Variables (6) %&
                       Values (j_rho, k_t, l_ye) - OpacityTable % EOSTable % &
                       DV % Offsets(6) - epsilon    !6 =Neutron Chemical Potential

                 xp  = 10**OpacityTable % EOSTable % DV % Variables (7) %&
                       Values (j_rho, k_t, l_ye) - OpacityTable % EOSTable % &
                       DV % Offsets(7) - epsilon    !7 =Proton Mass Fraction

                 xn  = 10**OpacityTable % EOSTable % DV % Variables (8) %&
                       Values (j_rho, k_t, l_ye) - OpacityTable % EOSTable % &
                       DV % Offsets(8) - epsilon    !8 =Neutron Mass Fraction

             xheavy  = 10**OpacityTable % EOSTable % DV % Variables (10) %&
                       Values (j_rho, k_t, l_ye) - OpacityTable % EOSTable % &
                       DV % Offsets(10) - epsilon   !10 =Heavy Mass Fraction

                 Z   = 10**OpacityTable % EOSTable % DV % Variables (11) %&
                       Values (j_rho, k_t, l_ye) - OpacityTable % EOSTable % &
                       DV % Offsets(11) - epsilon   !11 =Heavy Charge Number

                 A   = 10**OpacityTable % EOSTable % DV % Variables (12) %&
                       Values (j_rho, k_t, l_ye) - OpacityTable % EOSTable % &
                       DV % Offsets(12) - epsilon   !12 =Heavy Mass Number
 
                 bb  = (chem_e + chem_p - chem_n)/(T*kMev)
         
         Do i_r = 1, nOpacA
           DO i_e = 1, OpacityTable % nPointsE

              energy = OpacityTable % EnergyGrid % Values(i_e)

              OpacityTable % thermEmAb % Absorptivity(i_r) % Values (i_e, j_rho, k_t, l_ye) &
                     = totalECapEm(energy, rho, T, Z, A, chem_e, chem_n, chem_p, &
                       xheavy, xn, xp )
           END DO  !i_e
           
!
!           CALL GreyMomentWithGaussianQuadrature&
!                           ( nquadGrey, bb, &
!                             bufferquad1, "GreyMoment_Number ", .FALSE. )
!
!           CALL GreyMomentWithGaussianQuadrature&
!                           ( nquadGrey, bb, &
!                             bufferquad2, "GreyMoment_Energy ", .FALSE. )
!
!           CALL GreyOpacityWithGaussianQuadrature&
!                           ( nquadGrey, bb, &
!                             rho, T, Z, A, chem_e, chem_n,&
!                             chem_p, xheavy, xn, xp,&
!                             bufferquad3,"GreyOpacity_Number ", .FALSE. )
!
!           CALL GreyOpacityWithGaussianQuadrature&
!                           ( nquadGrey, bb, &
!                             rho, T, Z, A, chem_e, chem_n,&
!                             chem_p, xheavy, xn, xp,&
!                             bufferquad4,"GreyOpacity_Energy ", .FALSE. )
!
!           OpacityTable % thermEmAb % GreyMoment_Number_FD(i_r) % &
!                          Values ( j_rho, k_t, l_ye)  &
!              = bufferquad1 * (T*kMeV)**3 
!
!           OpacityTable % thermEmAb % GreyMoment_Energy_FD(i_r) % &
!                          Values ( j_rho, k_t, l_ye)  &
!              = bufferquad2 * (T*kMeV)**3
!
!           OpacityTable % thermEmAb % GreyOpacity_Number_FD(i_r) % &
!                         Values ( j_rho, k_t, l_ye)  &
!              = bufferquad3 * (T*kMeV)**3 
!
!           OpacityTable % thermEmAb % GreyOpacity_Energy_FD(i_r) % &
!                           Values ( j_rho, k_t, l_ye)  &
!              = bufferquad4 * (T*kMeV)**3
        END DO !i_r

!----------------  Scatt_Iso -----------------------
         DO i_rb = 1, nOpacB
           DO t_m = 1, nMomB
             DO i_e = 1, OpacityTable % nPointsE
 
               energy = OpacityTable % EnergyGrid % Values(i_e)

!              CALL GreyOpacityWithGaussianQuadrature_scattIso&
!                           ( nquadGrey, bb, &
!                            rho, T, xheavy, A, Z, xn, xp, t_m-1, &
!                            bufferquad23,"GreyOpacity_Number ", .FALSE. )
!
!              CALL GreyOpacityWithGaussianQuadrature_scattIso&
!                          ( nquadGrey, bb, &
!                            rho, T, xheavy, A, Z, xn, xp, t_m-1, &
!                            bufferquad24,"GreyOpacity_Energy ", .FALSE. )
! 
               OpacityTable % scatt_Iso % Kernel(i_rb) % Values &
                          ( i_e, j_rho, k_t, l_ye, t_m ) &
                = totalElasticScatteringKernel&
                  ( energy, rho, T, xheavy, A, Z, xn, xp, t_m-1 )
             END DO  !i_e

!             OpacityTable % scatt_Iso % GreyOpacity_Number_FD(i_rb) % &
!                            Values ( j_rho, k_t, l_ye, t_m)  &
!                = bufferquad23  * (T*kMeV)**3
!
!             OpacityTable % scatt_Iso % GreyOpacity_Energy_FD(i_rb) % &
!                             Values ( j_rho, k_t, l_ye, t_m)  &
!                = bufferquad24  * (T*kMeV)**3
!
           END DO !t_m 
         END DO !i_rb

       END DO  !j_rho
     END DO  !k_t
   END DO  !l_ye
  

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

          CALL TotalNESKernel&
              ( OpacityTable % EnergyGrid % Values, TMeV, chem_e, nquad, t_m-1, NESK)

          DO i_e = 1, nPointsE
          
            DO i_ep = 1, nPointsE
             
              OpacityTable % scatt_NES % Kernel(i_rb) % Values &
                           ( i_ep, i_e, k_t, i_eta, t_m)       &
               = NESK(i_e, i_ep) ! NESK was saved as NESK(e,ep)
 
            END DO ! i_ep
          END DO ! i_e
        END DO  !k_t
      END DO !i_eta
    END DO ! t_m
  END DO ! i_rb

   END ASSOCIATE ! EnergyGrid
!----------------------------------------------------------------------
!      Describe the Table ( give the real physical value )               
!----------------------------------------------------------------------

  CALL DescribeOpacityTable( OpacityTable )

! ---------------------------------------------------------------------
!          LOG the WHOLE table for storage
! ---------------------------------------------------------------------

  DO i_r = 1, nOpacA

     OpacityTable % thermEmAb % Absorptivity(i_r) % Values&
     = LOG10( OpacityTable % thermEmAb % Absorptivity(i_r) % &
              Values + OpacityTable % thermEmAb % Offset )

!     OpacityTable % thermEmAb % GreyMoment_Number_FD(i_r) % Values &
!     = LOG10 ( OpacityTable % thermEmAb % GreyMoment_Number_FD(i_r) % &
!              Values + OpacityTable % thermEmAb % Offset )
!
!     OpacityTable % thermEmAb % GreyMoment_Energy_FD(i_r) % Values &
!     = LOG10 ( OpacityTable % thermEmAb % GreyMoment_Energy_FD(i_r) % &
!              Values + OpacityTable % thermEmAb % Offset )
!
!     OpacityTable % thermEmAb % GreyOpacity_Number_FD(i_r) % Values &
!     = LOG10 ( OpacityTable % thermEmAb % GreyOpacity_Number_FD(i_r) % &
!              Values + OpacityTable % thermEmAb % Offset )
!
!     OpacityTable % thermEmAb % GreyOpacity_Energy_FD(i_r) % Values &
!     = LOG10 ( OpacityTable % thermEmAb % GreyOpacity_Energy_FD(i_r) % &
!              Values + OpacityTable % thermEmAb % Offset )
!
  END DO  !i_r

  DO i_rb = 1, nOpacB

    OpacityTable % scatt_Iso % Kernel(i_rb) % Values &
    = LOG10 ( OpacityTable % scatt_Iso % Kernel(i_rb) % Values &
              + OpacityTable % scatt_Iso % Offset )

!    OpacityTable % scatt_Iso % GreyOpacity_Number_FD(i_rb) % Values &
!    = LOG10 ( OpacityTable % scatt_Iso % GreyOpacity_Number_FD(i_rb) % Values &
!              + OpacityTable % scatt_Iso % Offset )
!
!    OpacityTable % scatt_Iso % GreyOpacity_Energy_FD(i_rb) % Values &
!    = LOG10 ( OpacityTable % scatt_Iso % GreyOpacity_Energy_FD(i_rb) % Values &
!              + OpacityTable % scatt_Iso % Offset )
!
  END DO !i_rb

  DO i_rb = 1, nOpacB_NES

    OpacityTable % scatt_NES % Kernel(i_rb) % Values &
    = LOG10 ( OpacityTable % scatt_NES % Kernel(i_rb) % Values &
              + OpacityTable % scatt_NES % Offset )

  END DO !i_rb

! -- write into hdf5 file

  CALL InitializeHDF( )
  CALL WriteOpacityTableHDF( OpacityTable, "wl-OP-LS220-20-40-100-Lower-T-nquad30-NoGrey.h5" )
  CALL FinalizeHDF( )
  
  WRITE (*,*) "HDF write successful"
!=============================================================

END PROGRAM wlCreateOpacityTable
