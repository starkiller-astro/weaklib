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
  USE HDF5
  USE wlGridModule, ONLY: MakeLogGrid  
  USE wlThermoStateModule
  USE wlDependentVariablesModule
  
  USE wlIOModuleHDF
  USE wlEnergyGridModule
  USE wlOpacityTableModule
  USE wlOpacityFieldsModule
  USE wlOpacityTableIOModuleHDF
  USE wlExtPhysicalConstantsModule
  USE B85
  USE wlExtNumericalModule, ONLY: epsilon
  USE GreyVariables
 
implicit none
    
   TYPE(OpacityTableType)  :: OpacityTable
   INTEGER                 :: nOpacA = 1
   INTEGER                 :: nOpacB = 1
   INTEGER                 :: nMomB  = 2
   INTEGER                 :: nOpacC = 0
   INTEGER                 :: nMomC  = 0
!-------------------------------------------------------------------------
! Set E grid limits
!-------------------------------------------------------------------------
   INTEGER                 :: nPointsE = 40 
   REAL(dp), PARAMETER     :: Emin = 1.0d-1  
   REAL(dp), PARAMETER     :: Emax = 3.0d02
   INTEGER, PARAMETER      :: nSpeciesA = 1
   

   INTEGER                 :: i_r, i_rb, i_e, j_rho, k_t, l_ye, t_m, i_quad
   REAL(dp)                :: energy, rho, T, ye, Z, A, chem_e, chem_n, &
                              chem_p, xheavy, xn, xp, bb, &
                              bufferquad1, bufferquad2, bufferquad3,&
                              bufferquad4, bufferquad21, bufferquad22, &
                              bufferquad23, bufferquad24
!   INTEGER, PARAMETER      :: nquad = 5
  INTEGER, PARAMETER      :: nquadGrey = 20

PRINT*, "Allocating OpacityTable"   

   CALL InitializeHDF( ) 
   CALL AllocateOpacityTable &
            ( OpacityTable, nOpacA, nOpacB, nMomB, &
              nOpacC, nMomC, nPointsE ) 
   CALL FinalizeHDF( )

! Set OpacityTableTypeA thermEmAb  
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

! Set OpacityTableTypeB scatt_Iso
   OpacityTable % scatt_Iso % nOpacities = nOpacB
   
   OpacityTable % scatt_Iso % nMoments = nMomB
   
   OpacityTable % scatt_Iso % nPoints(1) = nPointsE

   OpacityTable % scatt_Iso % nPoints(2:4) = OpacityTable % nPointsTS

   OpacityTable % scatt_Iso % Names = &
                                (/'Electron Iso-sca Absorptivity'/)
   OpacityTable % scatt_Iso % Species = &
                                (/'Electron Neutrino           '/)
   OpacityTable % scatt_Iso % Units = &
                                (/'Per Centimeter              '/)

   OpacityTable % scatt_Iso % Offset = 1.0d-1
!-----------------------------   
! Generate E grid from limits
!-----------------------------
PRINT*, "Making Energy Grid"

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

!-------------------------------------------------------------------------
!            Print The OpacityTable
!-------------------------------------------------------------------------

   CALL DescribeOpacityTable( OpacityTable )


!-------------------------------------------------------------------------
!              Fill OpacityTable
!-------------------------------------------------------------------------

!-----------------  ECAPEM ----------------------- 

   DO l_ye = 1, OpacityTable % nPointsTS(3)
     
        ye = OpacityTable % EOSTable % TS % States (3) % Values (l_ye)

      DO k_t = 1, OpacityTable % nPointsTS(2)

           T = OpacityTable % EOSTable % TS % States (2) % Values (k_t)

         DO j_rho = 1, OpacityTable % nPointsTS(1)

              rho = OpacityTable % EOSTable % TS % States (1) % Values (j_rho)

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
           

           CALL GreyMomentWithGaussianQuadrature&
                           ( nquadGrey, bb, &
                             bufferquad1, "GreyMoment_Number ", .FALSE. )

           CALL GreyMomentWithGaussianQuadrature&
                           ( nquadGrey, bb, &
                             bufferquad2, "GreyMoment_Energy ", .FALSE. )

           CALL GreyOpacityWithGaussianQuadrature&
                           ( nquadGrey, bb, &
                             rho, T, Z, A, chem_e, chem_n,&
                             chem_p, xheavy, xn, xp,&
                             bufferquad3,"GreyOpacity_Number ", .FALSE. )

           CALL GreyOpacityWithGaussianQuadrature&
                           ( nquadGrey, bb, &
                             rho, T, Z, A, chem_e, chem_n,&
                             chem_p, xheavy, xn, xp,&
                             bufferquad4,"GreyOpacity_Energy ", .FALSE. )

           OpacityTable % thermEmAb % GreyMoment_Number_FD(i_r) % &
                          Values ( j_rho, k_t, l_ye)  &
              = bufferquad1 * (T*kMeV)**3 

           OpacityTable % thermEmAb % GreyMoment_Energy_FD(i_r) % &
                          Values ( j_rho, k_t, l_ye)  &
              = bufferquad2 * (T*kMeV)**3

           OpacityTable % thermEmAb % GreyOpacity_Number_FD(i_r) % &
                         Values ( j_rho, k_t, l_ye)  &
              = bufferquad3 * (T*kMeV)**3 

           OpacityTable % thermEmAb % GreyOpacity_Energy_FD(i_r) % &
                           Values ( j_rho, k_t, l_ye)  &
              = bufferquad4 * (T*kMeV)**3
        END DO !i_r

!----------------  Scatt_Iso -----------------------
         DO i_rb = 1, nOpacB
           DO t_m = 1, nMomB

             DO i_e = 1, OpacityTable % nPointsE
 
               energy = OpacityTable % EnergyGrid % Values(i_e)

              CALL GreyOpacityWithGaussianQuadrature_scattIso&
                           ( nquadGrey, bb, &
                            rho, T, xheavy, A, Z, xn, xp, t_m-1, &
                            bufferquad23,"GreyOpacity_Number ", .FALSE. )

              CALL GreyOpacityWithGaussianQuadrature_scattIso&
                          ( nquadGrey, bb, &
                            rho, T, xheavy, A, Z, xn, xp, t_m-1, &
                            bufferquad24,"GreyOpacity_Energy ", .FALSE. )
 
               OpacityTable % scatt_Iso % Kernel(i_rb) % Values &
                          ( i_e, j_rho, k_t, l_ye, t_m ) &
                = totalElasticScatteringKernel&
                  ( energy, rho, T, xheavy, A, Z, xn, xp, t_m-1 )
             END DO  !i_e

             OpacityTable % scatt_Iso % GreyMoment_Number_FD(i_rb) % &
                            Values ( j_rho, k_t, l_ye, t_m)  &
                = bufferquad21 * (T*kMeV)**3
         
             OpacityTable % scatt_Iso % GreyMoment_Energy_FD(i_rb) % &
                            Values ( j_rho, k_t, l_ye, t_m)  &
                = bufferquad22 * (T*kMeV)**3
   
             OpacityTable % scatt_Iso % GreyOpacity_Number_FD(i_rb) % &
                            Values ( j_rho, k_t, l_ye, t_m)  &
                = bufferquad23  * (T*kMeV)**3

             OpacityTable % scatt_Iso % GreyOpacity_Energy_FD(i_rb) % &
                             Values ( j_rho, k_t, l_ye, t_m)  &
                = bufferquad24  * (T*kMeV)**3

           END DO !t_m 
         END DO !i_rb

       END DO  !j_rho
     END DO  !k_t
   END DO  !l_ye
  
!----------------------------------------------------------------------
!----------------------------------------------------------------------

  CALL DescribeOpacityTable( OpacityTable )

  DO i_r = 1, nOpacA

     OpacityTable % thermEmAb % Absorptivity(i_r) % Values&
     = LOG10( OpacityTable % thermEmAb % Absorptivity(i_r) % &
              Values + OpacityTable % thermEmAb % Offset )

     OpacityTable % thermEmAb % GreyMoment_Number_FD(i_r) % Values &
     = LOG10 ( OpacityTable % thermEmAb % GreyMoment_Number_FD(i_r) % &
              Values + OpacityTable % thermEmAb % Offset )

     OpacityTable % thermEmAb % GreyMoment_Energy_FD(i_r) % Values &
     = LOG10 ( OpacityTable % thermEmAb % GreyMoment_Energy_FD(i_r) % &
              Values + OpacityTable % thermEmAb % Offset )

     OpacityTable % thermEmAb % GreyOpacity_Number_FD(i_r) % Values &
     = LOG10 ( OpacityTable % thermEmAb % GreyOpacity_Number_FD(i_r) % &
              Values + OpacityTable % thermEmAb % Offset )

     OpacityTable % thermEmAb % GreyOpacity_Energy_FD(i_r) % Values &
     = LOG10 ( OpacityTable % thermEmAb % GreyOpacity_Energy_FD(i_r) % &
              Values + OpacityTable % thermEmAb % Offset )

  END DO  !i_r

  DO i_rb = 1, nOpacB

    OpacityTable % scatt_Iso % Kernel(i_rb) % Values &
    = LOG10 ( OpacityTable % scatt_Iso % Kernel(i_rb) % Values &
              + OpacityTable % scatt_Iso % Offset )

    OpacityTable % scatt_Iso % GreyMoment_Number_FD(i_rb) % Values &
    = LOG10 ( OpacityTable % scatt_Iso % GreyMoment_Number_FD(i_rb) % Values &
              + OpacityTable % scatt_Iso % Offset )

    OpacityTable % scatt_Iso % GreyMoment_Energy_FD(i_rb) % Values &
    = LOG10 ( OpacityTable % scatt_Iso % GreyMoment_Energy_FD(i_rb) % Values &
              + OpacityTable % scatt_Iso % Offset )

    OpacityTable % scatt_Iso % GreyOpacity_Number_FD(i_rb) % Values &
    = LOG10 ( OpacityTable % scatt_Iso % GreyOpacity_Number_FD(i_rb) % Values &
              + OpacityTable % scatt_Iso % Offset )

    OpacityTable % scatt_Iso % GreyOpacity_Energy_FD(i_rb) % Values &
    = LOG10 ( OpacityTable % scatt_Iso % GreyOpacity_Energy_FD(i_rb) % Values &
              + OpacityTable % scatt_Iso % Offset )

  END DO !i_rb

  CALL InitializeHDF( )
  CALL WriteOpacityTableHDF( OpacityTable, "wl-OP-LS220-20-40-100-Lower-T-nquad20-Grey.h5" )
  CALL FinalizeHDF( )
  
  WRITE (*,*) "HDF write successful"

!=============================================================

END PROGRAM wlCreateOpacityTable
