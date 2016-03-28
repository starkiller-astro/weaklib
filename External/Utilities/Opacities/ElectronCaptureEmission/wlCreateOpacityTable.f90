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
!
!   
!
!    CONTAINS:
!    
!
!    Modules used:
!       wlKindModule
!       HDF5
!       wlThermoStateModule
!       wlEOSIOModuleHDF
!       wlOpacityTableModule
!       wlDependentVariablesModule
!       wlGridModule, ONLY: MakeLinearGrid, MakeLogGrid
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
  USE wlOpacityTableModule
  USE wlOpacityTableIOModuleHDF
  USE wlExtPhysicalConstantsModule

implicit none
    

   TYPE(OpacityTableType)  :: OpacityTable
   REAL                    :: totalECapEm ! Functions 
!-------------------------------------------------------------------------
! Set E grid limits
!-------------------------------------------------------------------------
   INTEGER                 :: nPointsE = 40 
   REAL(dp), Parameter     :: Emin = 1.0d00  
   REAL(dp), Parameter     :: Emax = 3.0d02
   INTEGER, Parameter      :: nSpeciesA = 1
   

   INTEGER                 :: i_r, i_e, j_rho, k_t, l_ye 
   REAL(dp)                :: energy, rho, T, Z, A, mue, mun, mup, xh, xn, xp


PRINT*, "Allocate OpacityTable"   
 
   CALL AllocateOpacityTable( OpacityTable, nSpeciesA, nPointsE ) 

   OpacityTable % ECAPEM     % Names = &
                                &(/'Neutrino Emissivity '/)  

PRINT*, "Allocate OpacityTable Units "
 
   OpacityTable % ECAPEM     % Units = &
                                &(/'NEEDED   '/)  ! == 0229 == the Unit?
PRINT*, "Set E grid limits "

   OpacityTable % EnergyGrid % minValue = Emin  
   OpacityTable % EnergyGrid % maxValue = Emax 
   OpacityTable % EnergyGrid % nPointsE = nPointsE

   
!-------------------------------------------------------------------------
! Generate E grid from limits
!-------------------------------------------------------------------------
PRINT*, "Make Energy Grid"

   CALL MakeLogGrid( Emin, Emax, nPointsE, &
                           & OpacityTable % EnergyGrid % Values)

!-------------------------------------------------------------------------
!            Printe The OpacityTable
!-------------------------------------------------------------------------
PRINT*, "Print The OpacityTable"

   CALL DescribeOpacityTable( OpacityTable )


!-------------------------------------------------------------------------
!              Fill OpacityTable
!-------------------------------------------------------------------------

!-----------------  ECAPEM ----------------------- 
 !  DO i_r = 1, SIZE(ECAPEM)
     i_r = 1;  ! i_r can be 1,2,3,4 
     DO l_ye = 1, OpacityTable % EOSTable % TS % nPoints(3)

       DO k_t = 1, OpacityTable % EOSTable % TS % nPoints(2)

              T = OpacityTable % EOSTable % TS % States (2) % Values (k_t)

         DO j_rho = 1, OpacityTable % EOSTable % TS % nPoints(1)

              rho = OpacityTable % EOSTable % TS % States (1) % Values (j_rho)

           DO i_e = 1, OpacityTable % EnergyGrid % nPointsE

              energy = OpacityTable % EnergyGrid % Values(i_e)
              mue = OpacityTable % EOSTable % DV % Variables (4) % Values (j_rho, k_t, l_ye)  !4 =Electron Chemical Potential             
              mup = OpacityTable % EOSTable % DV % Variables (5) % Values (j_rho, k_t, l_ye)  !5 =Proton Chemical Potential 
              mun = OpacityTable % EOSTable % DV % Variables (6) % Values (j_rho, k_t, l_ye)  !6 =Neutron Chemical Potential 
              xp  = OpacityTable % EOSTable % DV % Variables (7) % Values (j_rho, k_t, l_ye)  !7 =Proton Mass Fraction 
              xn  = OpacityTable % EOSTable % DV % Variables (8) % Values (j_rho, k_t, l_ye)  !8 =Neutron Mass Fraction 
              xh  = OpacityTable % EOSTable % DV % Variables (10) % Values (j_rho, k_t, l_ye) !10 =Heavy Mass Fraction
              Z   = OpacityTable % EOSTable % DV % Variables (11) % Values (j_rho, k_t, l_ye) !11 =Heavy Charge Number
              A   = OpacityTable % EOSTable % DV % Variables (12) % Values (j_rho, k_t, l_ye) !12 =Heavy Mass Number 

              OpacityTable % ECAPEM(i_r) % Values (i_e, j_rho, k_t, l_ye) &
                  = totalECapEm(energy, rho, T, Z, A, mue, mun, mup, xh, xn, xp )
          

           END DO
         END DO
       END DO
     END DO
  ! END DO


PRINT*, "Print The OpacityTable"

   CALL DescribeOpacityTable( OpacityTable )


!==================== Below Necessary ? ======================

!  CALL WriteEquationOfStateTableHDF( EOSTable )
!  CALL WriteEquationOfStateTableHDF( EOSTable, Description )

!  CALL FinalizeHDF( )

!  WRITE (*,*) "HDF write successful"

!=============================================================

END PROGRAM wlCreateOpacityTable





