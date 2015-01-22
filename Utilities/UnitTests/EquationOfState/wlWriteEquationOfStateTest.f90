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

  EOSTable % TS % nPoints(1:3) = nPoints(1:3)
  EOSTable % TS % Names(1:3) = (/'Density                         ',&
                                 'Temperature                     ',&
                                 'Electron Fraction               '/)

  EOSTable % TS % minValues(1:3) =  (/1.0d06,0.1d00,1.0d-02/)
  EOSTable % TS % maxValues(1:3) =  (/1.0d15,1.0d02,6.1d-01/)

  CALL MakeLogGrid( EOSTable % TS % minValues(1), EOSTable % TS % maxValues(1),&
         EOSTable % TS % nValues(1), EOSTable % TS % States(1) % Values)
  CALL MakeLogGrid( EOSTable % TS % minValues(2), EOSTable % TS % maxValues(2),&
         EOSTable % TS % nValues(2), EOSTable % TS % States(2) % Values)
  CALL MakeLinearGrid( EOSTable % TS % minValues(3), EOSTable % TS % maxValues(3),&
         EOSTable % TS % nValues(3), EOSTable % TS % States(3) % Values)

  EOSTable % DV % Names(1:12) = (/'Pressure                        ', &
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

  EOSTable % DV % Units(1:12) = (/'Dynes per cm^2                  ', &
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

  DO i = 1, SIZE( EOSTable % DV % Variables )
    WRITE (*,*) SHAPE( EOSTable % DV % Variables(i) % Values )
  END DO

  DO i = 1, SIZE( EOSTable % DV % Variables )
    EOSTable % DV % Variables(i) % Values = i 
  END DO


  !DO j = 1,12
  !  WRITE(*,*)
  !  WRITE(*,*) TRIM( DV % Names(j) ) , j
  !  WRITE(*,*)
  !  WRITE(*,*) DV % Variables(j) % Values(:,:,:)
  !END DO

print*,"4"

  CALL InitializeHDF( )

  CALL OpenFileHDF( "EquationOfStateTable.h5", .true., file_id )

  CALL OpenGroupHDF( "ThermoState", .true., file_id, group_id )
  CALL WriteThermoStateHDF( EOSTable % TS, group_id )
  CALL CloseGroupHDF( group_id )

  CALL OpenGroupHDF( "DependentVariables", .true., file_id, group_id )
  CALL WriteDependentVariablesHDF( EOSTable % DV, group_id )
  CALL CloseGroupHDF( group_id )

  CALL CloseFileHDF( file_id )
  CALL DeAllocateThermoState( EOSTable % TS )
  CALL DeAllocateDependentVariables( EOSTable % DV )

  CALL FinalizeHDF( )

END PROGRAM wlWriteEquationOfStateTest
