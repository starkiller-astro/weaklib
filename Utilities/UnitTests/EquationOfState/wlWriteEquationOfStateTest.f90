PROGRAM wlWriteEquationOfStateTest
 
  USE wlDependentVariablesModule
  USE HDF5
  USE wlThermoStateModule
  USE wlEquationOfStateTableModule
  USE wlIOModuleHDF, ONLY: InitializeHDF, OpenFileHDF, OpenGroupHDF,         &
                           WriteDependentVariablesHDF, CloseGroupHDF,        & 
                           ReadDependentVariablesHDF, CloseFileHDF,          &
                           FinalizeHDF, ReadDimensionsHDF,                   &
                           LoadDependentVariablesHDF, ReadNumberVariablesHDF
  implicit none

  INTEGER :: i
  INTEGER, DIMENSION(3) :: nPoints
  INTEGER :: nVariables
  INTEGER :: j
  TYPE(EquationOfStateTableType) :: EOSTable
  TYPE(DependentVariablesType) :: DV, DV2 
  INTEGER(HID_T) :: file_id
  INTEGER(HID_T) :: group_id


  nPoints = (/10,10,10/)
  nVariables = 12
  CALL AllocateEquationOfStateTable( EOSTable, nPoints , nVariables )

STOP

print*,"2"

! Insert stuff from ThermoStateTest
  DV % Names(1:12) = (/'Pressure                        ', &
                       'Entropy Per Baryon              ', &
                       'Internal Energy Density         ', &
                       'Neutron Chemical Potential      ', &
                       'Electron Chemical Potential     ', &
                       'Proton Chemical Potential       ', &
                       'Neutron Mass Fraction           ', &
                       'Proton Mass Fraction            ', &
                       'Helium Mass Fraction            ', &
                       'Heavy Mass Fraction             ', &
                       'Heavy Mass Number               ', &
                       'Heavy Charge Number             '/)

  DV % Units(1:12) = (/'Dynes per cm^2                  ', &
                       'Entropy Per Baryon Units        ', &
                       'Internal Energy Density         ', &
                       'Neutron Chemical Potential      ', &
                       'Electron Chemical Potential     ', &
                       'Proton Chemical Potential       ', &
                       'Neutron Mass Fraction           ', &
                       'Proton Mass Fraction            ', &
                       'Helium Mass Fraction            ', &
                       'Heavy Mass Fraction             ', &
                       'Heavy Mass Number               ', &
                       'Heavy Charge Number             '/)


print*,"3"

  DO i = 1, SIZE( DV % Variables )
    WRITE (*,*) SHAPE( DV % Variables(i) % Values )
  END DO

  DO i = 1, SIZE( DV % Variables )
    DV % Variables(i) % Values = i 
  END DO


  !DO j = 1,12
  !  WRITE(*,*)
  !  WRITE(*,*) TRIM( DV % Names(j) ) , j
  !  WRITE(*,*)
    !WRITE(*,*) DV2 % nValues(j)
  !  WRITE(*,*) DV % Variables(j) % Values(:,:,:)
  !END DO

print*,"4"

  CALL InitializeHDF( )

  CALL OpenFileHDF( "EquationOfStateTable.h5", .true., file_id )

! Add next three lines, but for thermostate

  CALL OpenGroupHDF( "DependentVariables", .true., file_id, group_id )
  CALL WriteDependentVariablesHDF( DV, group_id )
  CALL CloseGroupHDF( group_id )

  CALL CloseFileHDF( file_id )
  CALL DeAllocateDependentVariables( DV )

  CALL FinalizeHDF( )

END PROGRAM wlWriteEquationOfStateTest
