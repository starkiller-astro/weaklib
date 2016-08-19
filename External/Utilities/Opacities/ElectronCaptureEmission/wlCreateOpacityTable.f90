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
  USE GaussianQuadrature
 
implicit none
    
   TYPE(OpacityTableType)  :: OpacityTable
   INTEGER                 :: nOpacA = 1
   INTEGER                 :: nOpacB = 0
   INTEGER                 :: nMomB  = 0
   INTEGER                 :: nOpacC = 0
   INTEGER                 :: nMomC  = 0
!-------------------------------------------------------------------------
! Set E grid limits
!-------------------------------------------------------------------------
   INTEGER                 :: nPointsE = 40 
   REAL(dp), PARAMETER     :: Emin = 1.0d-1  
   REAL(dp), PARAMETER     :: Emax = 3.0d02
   INTEGER, PARAMETER      :: nSpeciesA = 1
   

   INTEGER                 :: i_r, i_e, j_rho, k_t, l_ye, i_quad 
   REAL(dp)                :: energy, rho, T, Z, A, chem_e, chem_n, &
                              chem_p, xheavy, xn, xp, bb, &
                              bufferquad1, bufferquad2, bufferquad3,&
                              bufferquad4
   INTEGER, PARAMETER      :: nquad = 6

PRINT*, "Allocating OpacityTable"   

   CALL InitializeHDF( ) 
   CALL AllocateOpacityTable &
            ( OpacityTable, nOpacA, nOpacB, nMomB, &
              nOpacC, nMomC, nPointsE ) 
   CALL FinalizeHDF( )
  
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

   DO i_r = 1, nOpacA
 
     DO l_ye = 1, 2!OpacityTable % nPointsTS(3)

       DO k_t = 1, 2!OpacityTable % nPointsTS(2)

             T = OpacityTable % EOSTable % TS % States (2) % Values (k_t)

         DO j_rho = 1, 2!OpacityTable % nPointsTS(1)

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

                 Z   = OpacityTable % EOSTable % DV % Variables (11) %&
                       Values (j_rho, k_t, l_ye) - OpacityTable % EOSTable % &
                       DV % Offsets(11) - epsilon   !11 =Heavy Charge Number

                 A   = OpacityTable % EOSTable % DV % Variables (12) %&
                       Values (j_rho, k_t, l_ye) - OpacityTable % EOSTable % &
                       DV % Offsets(12) - epsilon   !12 =Heavy Mass Number 

           DO i_e = 1, 2!OpacityTable % nPointsE

              energy = OpacityTable % EnergyGrid % Values(i_e)

              OpacityTable % thermEmAb % Absorptivity(i_r) % Values (i_e, j_rho, k_t, l_ye) &
               = totalECapEm(energy, rho, T, Z, A,&
                      chem_e, chem_n, chem_p, &
                      xheavy, xn, xp )
           END DO  !i_e
           
           bb = (chem_e + chem_p - chem_n)/(T*kMev)

           CALL GreyMomentWithGaussianQuadrature&
                           ( nquad, bb, &
                             bufferquad1, "GreyMoment_Number ", .FALSE. )

           CALL GreyMomentWithGaussianQuadrature&
                           ( nquad, bb, &
                             bufferquad2, "GreyMoment_Energy ", .FALSE. )

           CALL GreyOpacityWithGaussianQuadrature&
                           ( nquad, bb, &
                             rho, T, Z, A, chem_e, chem_n,&
                             chem_p, xheavy, xn, xp,&
                             bufferquad3,"GreyOpacity_Number ", .FALSE. )

           CALL GreyOpacityWithGaussianQuadrature&
                           ( nquad, bb, &
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

         END DO  !j_rho
       END DO  !k_t
     END DO  !l_ye
   END DO  !i_r
  
  
  CALL DescribeOpacityTable( OpacityTable )

   DO i_r = 1, nOpacA
     DO l_ye = 1, OpacityTable % nPointsTS(3)
       DO k_t = 1, OpacityTable % nPointsTS(2)
         DO j_rho = 1, OpacityTable % nPointsTS(1)
           DO i_e = 1, OpacityTable % nPointsE
             OpacityTable % thermEmAb % Absorptivity(i_r) % &
                          Values (i_e, j_rho, k_t, l_ye)  &
             = LOG10( OpacityTable % thermEmAb % Absorptivity(i_r) % &
                          Values (i_e, j_rho, k_t, l_ye) + &
                      OpacityTable % thermEmAb % Offset )
           END DO  !i_e
         END DO  !j_rho
       END DO  !k_t
     END DO  !l_ye
   END DO  !i_r

  CALL InitializeHDF( )
  CALL WriteOpacityTableHDF( OpacityTable, "wl-OP-LS220-20-40-100.h5" )
  CALL FinalizeHDF( )
  
  WRITE (*,*) "HDF write successful"

!=============================================================

END PROGRAM wlCreateOpacityTable
