PROGRAM wlWriteEquationOfStateTest
 
  USE wlKindModule, ONLY: dp
  USE HDF5
  USE wlExtEOSWrapperModule, ONLY: wlGetFullEOS
  USE wlGridModule, ONLY: MakeLinearGrid, MakeLogGrid
  USE wlEquationOfStateTableModule
  USE wlIOModuleHDF, ONLY: InitializeHDF, OpenFileHDF, OpenGroupHDF,          &
                           WriteDependentVariablesHDF, CloseGroupHDF,         & 
                           CloseFileHDF, FinalizeHDF, WriteThermoStateHDF
                           !Replace all of the above with WriteEquationOfStateTableHDF

  implicit none

  INTEGER                        :: i, j, k, l
  INTEGER, DIMENSION(3)          :: nPoints
  INTEGER                        :: nVariables
  TYPE(EquationOfStateTableType) :: EOSTable
  INTEGER(HID_T)                 :: file_id
  INTEGER(HID_T)                 :: group_id

  CHARACTER(len=1)   :: EOSFlag     ! nuclear eos selection flag
  CHARACTER(len=3)   :: LScompress
  CHARACTER(len=128) :: LSFilePath

  LOGICAL            :: fail        ! did EoS fail to converge

  nPoints = (/81,31,47/)
  nVariables = 13
  LScompress = '220'
  LSFilePath = '../../../External/LS/Data'

  CALL wlExtInitializeEOS( LSFilePath, LScompress )

PRINT*, "Allocate EOS"
  CALL AllocateEquationOfStateTable( EOSTable, nPoints , nVariables )

  EOSTable % TS % Names(1:3) = (/'Density                         ',&
                                 'Temperature                     ',&
                                 'Electron Fraction               '/)

PRINT*, "Allocate Independent Variable Units " 

  EOSTable % TS % Units(1:3) = (/'Grams per cm^3                  ', &
                                 'K                               ', &
                                 '                                '/) 

  EOSTable % TS % minValues(1:3) =  (/1.0d07, 1.0d09, 0.05d0/)
  EOSTable % TS % maxValues(1:3) =  (/1.0d15, 1.0d12, 0.51d0/)

!------------------------------------------------------------------------------
! Generate rho, T, Ye grid from limits
!------------------------------------------------------------------------------
PRINT*, "Make Grids"
  CALL MakeLogGrid( EOSTable % TS % minValues(1), EOSTable % TS % maxValues(1),&
         EOSTable % TS % nPoints(1), EOSTable % TS % States(1) % Values)
  CALL MakeLogGrid( EOSTable % TS % minValues(2), EOSTable % TS % maxValues(2),&
         EOSTable % TS % nPoints(2), EOSTable % TS % States(2) % Values)
  CALL MakeLinearGrid( EOSTable % TS % minValues(3), EOSTable % TS % maxValues(3),&
         EOSTable % TS % nPoints(3), EOSTable % TS % States(3) % Values)

PRINT*, "Allocate Names " 
  EOSTable % DV % Names(1:13) = (/'Pressure                        ', &
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
                                  'Heavy Binding Energy            '/)


PRINT*, "Allocate Dependent Variable Units " 
  EOSTable % DV % Units(1:13) = (/'Dynes per cm^2                  ', &
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
                                  'MeV                             '/)
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
    be_heavy    => EOSTable % DV % Variables(13) % Values(:,:,:) ) 

  EOSFlag = "L" 

  DO k = 1, EOSTable % nPoints(3) 
    DO j = 1, EOSTable % nPoints(2)
      DO i = 1, EOSTable % nPoints(1) 
         ! PRINT*, "i,j,k = ", i,j,k
          CALL wlGetFullEOS( Density(i), Temperature(j), Ye(k), EOSFlag, fail,      &
                       press(i,j,k), entrop(i,j,k), energ(i,j,k), chem_e(i,j,k),    &
                       chem_p(i,j,k), chem_n(i,j,k), xn_prot(i,j,k), xn_neut(i,j,k),&
                       xn_alpha(i,j,k), xn_heavy(i,j,k), z_heavy(i,j,k),            &
                       a_heavy(i,j,k), be_heavy(i,j,k) )   
      END DO
    END DO
  END DO

  END ASSOCIATE


  DO l = 1, nVariables 
    EOSTable % DV % Offsets(l) & 
      = MAX( 0.0_dp, MINVAL( EOSTable % DV % Variables(l) % Values(:,:,:) ) ) 
  END DO 

! Turn what follows into a WriteEquationOfStateTableHDF( EOSTable, Name of table ) subroutine in wlIOModul
  CALL InitializeHDF( )

  CALL OpenFileHDF( "EquationOfStateTable.h5", .true., file_id )

  CALL OpenGroupHDF( "ThermoState", .true., file_id, group_id )
  CALL WriteThermoStateHDF( EOSTable % TS, group_id )
  CALL CloseGroupHDF( group_id )

  CALL OpenGroupHDF( "DependentVariables", .true., file_id, group_id )
  CALL WriteDependentVariablesHDF( EOSTable % DV, group_id )
  CALL CloseGroupHDF( group_id )

  CALL CloseFileHDF( file_id )

  CALL FinalizeHDF( )

  WRITE (*,*) "HDF write successful"

  CALL DeAllocateEquationOfStateTable( EOSTable )

END PROGRAM wlWriteEquationOfStateTest
