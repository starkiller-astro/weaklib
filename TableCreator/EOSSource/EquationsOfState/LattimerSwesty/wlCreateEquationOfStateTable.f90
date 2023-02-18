PROGRAM wlCreateEquationOfStateTable
 
  USE wlKindModule, ONLY: dp
  USE HDF5
  USE wlExtEOSWrapperModule, ONLY: wlGetFullEOS
  USE wlGridModule, ONLY: MakeLinearGrid, MakeLogGrid
  USE wlEquationOfStateTableModule
  USE wlIOModuleHDF
  USE wlEOSIOModuleHDF
  USE wlExtPhysicalConstantsModule, ONLY: kMeV, cvel
                  
                 

  implicit none

  INTEGER                        :: i, j, k, l, kmax, count
  INTEGER                        :: count_LS_call_fails
  INTEGER                        :: count_unphysical
  INTEGER, DIMENSION(3)          :: nPoints
  INTEGER                        :: nVariables
  TYPE(EquationOfStateTableType) :: EOSTable
  INTEGER(HID_T)                 :: file_id
  INTEGER(HID_T)                 :: group_id

  CHARACTER(len=1)   :: EOSFlag     ! nuclear eos selection flag
  CHARACTER(len=3)   :: LScompress
  CHARACTER(len=128) :: LSFilePath

  LOGICAL            :: fail        ! did EoS fail to converge

  REAL(8)                                       :: Ye_tmp
  REAL(8), DIMENSION(:), ALLOCATABLE            :: Ye_save
 
  REAL(dp) :: total_fraction

  94 FORMAT ("rho=", es12.5,1x, "T=", es12.5,1x, "Ye=" , es12.5 )
  95 FORMAT ("Press=", es12.5,1x, "Entropy=", es12.5,1x, "Energy=" , es12.5 )
  96 FORMAT ("Elec Chem Pot=", es12.5,1x, "Prot Chem Pot=", es12.5,1x,         &
             "Neut Chem Pot=" , es12.5 )
  97 FORMAT ("Prot Mass Fract=", es12.5,1x, "Neut Mass Fract=", es12.5,1x,     &
             "Alpha Mass Fract=" , es12.5 )
  98 FORMAT ("Heavy Mass Fract=", es12.5,1x, "Heavy Charge # =", es12.5,1x,    &
             "Heavy Mass #=" , es12.5,1x, "Heavy Binding Energy=", es12.5 )

!  FOR FUTURE USE OF INPUT FILE
!  FileName = "input.d" 
!  OPEN( newunit = InputUnit, file = FileName )
!  READ( InputUnit, * ) Description, Resolution 

!  nPoints = (/81,24,24/) ! Low Res
!  nPoints = (/65,55,13/) ! Low Res, 8pts/dec, 20pts/dec, 25pts/100  
!  nPoints = (/151,47,49/) ! Low Res
!  nPoints = (/81,500,24/) ! High Res in T only
!  nPoints = (/161,93,25/) ! High Res in T only
!  nPoints = (/161,47,25/) ! Standard C Res
!  nPoints = (/161,47,49/) ! Hi Res in Ye
!  nPoints = (/161,108,49/) ! Standard D Res
  nPoints = (/201,109,1/)
!  nPoints = (/1,93,3/) ! one line
!  nPoints = (/321,47,25/) ! High Res in Rho
!  nPoints = (/321,93,49/) ! High Res
  nVariables = 15
  LScompress = '220'
  LSFilePath = '../../../../Old/LS/Data/'

  CALL wlExtInitializeEOS( LSFilePath, LScompress )

PRINT*, "Allocate EOS"
  CALL AllocateEquationOfStateTable( EOSTable, nPoints , nVariables )

  EOSTable % TS % Names(1:3) = (/'Density                         ',&
                                 'Temperature                     ',&
                                 'Electron Fraction               '/)

  EOSTable % TS % Indices % iRho = 1
  EOSTable % TS % Indices % iT   = 2
  EOSTable % TS % Indices % iYe  = 3

PRINT*, "Allocate Independent Variable Units " 

  EOSTable % TS % Units(1:3) = (/'Grams per cm^3                  ', &
                                 'K                               ', &
                                 '                                '/) 

 EOSTable % TS % minValues(1:3) =  (/1.0d07, 10.d0**9.3, 0.51d0/)
 !EOSTable % TS % maxValues(1:3) =  (/1.0d15, 1.0d12, 0.54d0/)
 EOSTable % TS % maxValues(1:3) =  (/1.0d15, 1.0d12, 0.51d0/)


!------------------------------------------------------------------------------
! Generate rho, T, Ye grid from limits
!------------------------------------------------------------------------------
PRINT*, "Make Grids"
  CALL MakeLogGrid( EOSTable % TS % minValues(1), EOSTable % TS % maxValues(1), &
         EOSTable % TS % nPoints(1), EOSTable % TS % States(1) % Values)
  CALL MakeLogGrid( EOSTable % TS % minValues(2), EOSTable % TS % maxValues(2), &
         EOSTable % TS % nPoints(2), EOSTable % TS % States(2) % Values)
  CALL MakeLinearGrid( EOSTable % TS % minValues(3), EOSTable % TS % maxValues(3), &
         EOSTable % TS % nPoints(3), EOSTable % TS % States(3) % Values)

PRINT*, "Label Grid Type"
  EOSTable % TS % LogInterp(1:3) =  (/1, 1, 0/)
  
PRINT*, "Allocate Names " 
  EOSTable % DV % Names(1:15) = (/'Pressure                        ', &
                                  'Entropy Per Baryon              ', &
                                  'Internal Energy Density         ', &
                                  'Electron Chemical Potential     ', &
                                  'Proton Chemical Potential       ', &
                                  'Neutron Chemical Potential      ', &
                                  'Proton Mass Fraction            ', &
                                  'Neutron Mass Fraction           ', &
                                  'Alpha Mass Fraction             ', &
                                  'Heavy Mass Fraction             ', &
                                  'Heavy Charge Number             ', &
                                  'Heavy Mass Number               ', &
                                  'Heavy Binding Energy            ', &
                                  'Thermal Energy                  ', &
                                  'Gamma1                          '/)

PRINT*, "Set Dependent Variable Identifier Indicies " 

    EOSTable % DV % Indices % iPressure = 1
    EOSTable % DV % Indices % iEntropyPerBaryon = 2
    EOSTable % DV % Indices % iInternalEnergyDensity = 3
    EOSTable % DV % Indices % iElectronChemicalPotential = 4
    EOSTable % DV % Indices % iProtonChemicalPotential = 5
    EOSTable % DV % Indices % iNeutronChemicalPotential = 6
    EOSTable % DV % Indices % iProtonMassFraction = 7
    EOSTable % DV % Indices % iNeutronMassFraction = 8
    EOSTable % DV % Indices % iAlphaMassFraction = 9
    EOSTable % DV % Indices % iHeavyMassFraction = 10
    EOSTable % DV % Indices % iHeavyChargeNumber = 11
    EOSTable % DV % Indices % iHeavyMassNumber = 12
    EOSTable % DV % Indices % iHeavyBindingEnergy = 13
    EOSTable % DV % Indices % iThermalEnergy = 14
    EOSTable % DV % Indices % iGamma1 = 15

      WRITE (*,*) "iRho, iGamma1", EOSTable % TS % Indices % iRho, EOSTable % DV % Indices % iGamma1

PRINT*, "Allocate Dependent Variable Units " 
  EOSTable % DV % Units(1:15) = (/'Dynes per cm^2                  ', &
                                          'k_b per baryon                  ', &
                                          'erg per gram                    ', &
                                          'MeV                             ', &
                                          'MeV                             ', &
                                          'MeV                             ', &
                                          '                                ', &
                                          '                                ', &
                                          '                                ', &
                                          '                                ', &
                                          '                                ', &
                                          '                                ', &
                                          'MeV                             ', &
                                          'MeV                             ', &
                                          '                                '/)

PRINT*, "Allocate Dependent Variable Logical " 
  EOSTable % DV % Repaired(:,:,:) = 0
  
PRINT*, "Begin Associate" 
  ASSOCIATE(&
    Density     => EOSTable % TS % States(1) % Values, &
    Temperature => EOSTable % TS % States(2) % Values, &
    Ye          => EOSTable % TS % States(3) % Values, &
    press       => EOSTable % DV % Variables(1) % Values(:,:,:), &
    entrop      => EOSTable % DV % Variables(2) % Values(:,:,:), &
    energ       => EOSTable % DV % Variables(3) % Values(:,:,:), &
    chem_e      => EOSTable % DV % Variables(4) % Values(:,:,:), &
    chem_p      => EOSTable % DV % Variables(5) % Values(:,:,:), &
    chem_n      => EOSTable % DV % Variables(6) % Values(:,:,:), &
    xn_prot     => EOSTable % DV % Variables(7) % Values(:,:,:), &
    xn_neut     => EOSTable % DV % Variables(8) % Values(:,:,:), &
    xn_alpha    => EOSTable % DV % Variables(9) % Values(:,:,:), &
    xn_heavy    => EOSTable % DV % Variables(10) % Values(:,:,:), &
    z_heavy     => EOSTable % DV % Variables(11) % Values(:,:,:), &
    a_heavy     => EOSTable % DV % Variables(12) % Values(:,:,:), &
    be_heavy    => EOSTable % DV % Variables(13) % Values(:,:,:), &
    thermalenergy    => EOSTable % DV % Variables(14) % Values(:,:,:), &  
    gamma1    => EOSTable % DV % Variables(15) % Values(:,:,:) ) 

  EOSFlag = "L" 

  count = 0
  count_LS_call_fails = 0
  count_unphysical    = 0

  EOSTable % DV % Repaired(:,:,:) = 0

  ALLOCATE( Ye_save( EOSTable % nPoints(3) ) ) 


      DO i = 1, EOSTable % nPoints(1) 
       DO k = EOSTable % nPoints(3), 1, -1
        DO j = EOSTable % nPoints(2), 1, -1
          Ye_tmp = Ye(k)
          CALL wlGetFullEOS( Density(i), Temperature(j), Ye_tmp, EOSFlag, fail,      &
                       press(i,j,k), energ(i,j,k), entrop(i,j,k), chem_n(i,j,k),    &
                       chem_p(i,j,k), chem_e(i,j,k), xn_neut(i,j,k), xn_prot(i,j,k),&
                       xn_alpha(i,j,k), xn_heavy(i,j,k), a_heavy(i,j,k),            &
                       z_heavy(i,j,k), be_heavy(i,j,k), thermalenergy(i,j,k), gamma1(i,j,k), i, j, k )   

          total_fraction = xn_prot(i,j,k) + xn_neut(i,j,k) + xn_alpha(i,j,k) + xn_heavy(i,j,k)

          IF ( energ(i,j,k) < 0.0d0 .or. entrop(i,j,k) < 0.0d0 .or. press(i,j,k) < 0.0d0   &
            .or. xn_prot(i,j,k) < 0.0d0 .or. xn_prot(i,j,k) > 1.0d0 &
            .or. xn_neut(i,j,k) < 0.0d0 .or. xn_neut(i,j,k) > 1.0d0 &
            .or. xn_alpha(i,j,k) < 0.0d0 .or. xn_alpha(i,j,k) > 1.0d0 &
            .or. xn_heavy(i,j,k) < 0.0d0 .or. xn_heavy(i,j,k) > 1.0d0 &
            .or. a_heavy(i,j,k) < 0.0d0 .or. z_heavy(i,j,k) < 0.0d0 &
            .or. ABS(total_fraction-1.0d0) > 1d-10 &
            .or. gamma1(i,j,k) < 0.0d0 ) THEN
            write(*,'(A15,f5.2,A8,es12.3,A6,f5.2,A7,f5.2)') 'Total fraction ', &
            total_fraction, ' rho = ', Density(i), ' T = ', Temperature(j)*kMeV, &
            ' Ye = ', Ye_tmp
            write(*,'(A8,es12.3)') 'energ = ', energ(i,j,k) 
            write(*,'(A6,es12.3)') 'ent = ', entrop(i,j,k) 
            write(*,'(A8,es12.3)') 'press = ', press(i,j,k) 
            write(*,'(A7,es12.3)') 'xn_n = ', xn_neut(i,j,k) 
            write(*,'(A7,es12.3)') 'xn_p = ', xn_prot(i,j,k) 
            write(*,'(A7,es12.3)') 'xn_a = ', xn_alpha(i,j,k) 
            write(*,'(A7,es12.3)') 'xn_h = ', xn_heavy(i,j,k) 
            write(*,'(A5,es12.3)') 'ah = ', a_heavy(i,j,k) 
            write(*,'(A5,es12.3)') 'zh = ', z_heavy(i,j,k) 
            write(*,'(A9,es12.3)') 'Gamma1 = ', gamma1(i,j,k) 
            count_unphysical = count_unphysical + 1
            EOSTable % DV % Repaired(i,j,k) = -1
            !count = count + 1
          END IF
          IF ( fail ) THEN 
            EOSTable % DV % Repaired(i,j,k) = -1
            count_LS_call_fails = count_LS_call_fails + 1
          ENDIF
               
      END DO
    END DO
  END DO

WRITE (*,*) count_LS_call_fails, " fails from LS call out of " , nPoints(1)*nPoints(2)*nPoints(3) 
WRITE (*,*) count_unphysical, " unphysical output out of " , nPoints(1)*nPoints(2)*nPoints(3) 

  END ASSOCIATE

  CALL InitializeHDF( )

  WRITE (*,*) "Starting HDF write "
  CALL WriteEquationOfStateTableHDF( EOSTable )

  CALL FinalizeHDF( )

  WRITE (*,*) "HDF write successful"

  CALL DeAllocateEquationOfStateTable( EOSTable )

END PROGRAM wlCreateEquationOfStateTable
