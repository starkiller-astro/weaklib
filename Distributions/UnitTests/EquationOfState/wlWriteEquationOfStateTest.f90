PROGRAM wlWriteEquationOfStateTest
 
  USE wlDependentVariablesModule
  USE HDF5
  USE wlThermoStateModule
  USE wlGridModule, ONLY: MakeLinearGrid, MakeLogGrid
  USE wlEquationOfStateTableModule
  USE wlIOModuleHDF
  USE wlEOSIOModuleHDF

  implicit none

  INTEGER                        :: i
  INTEGER, DIMENSION(3)          :: nPoints
  INTEGER                        :: nVariables
  INTEGER                        :: j
  TYPE(EquationOfStateTableType) :: EOSTable
  INTEGER(HID_T)                 :: file_id
  INTEGER(HID_T)                 :: group_id


  nPoints = (/10,10,10/)
  nVariables = 13

  CALL AllocateEquationOfStateTable( EOSTable, nPoints , nVariables )

  EOSTable % TS % Names(1:3) = (/'Density                         ',&
                                 'Temperature                     ',&
                                 'Electron Fraction               '/)

  EOSTable % TS % minValues(1:3) =  (/1.0d06,0.1d00,1.0d-02/)
  EOSTable % TS % maxValues(1:3) =  (/1.0d15,1.0d02,6.1d-01/)

  CALL MakeLogGrid( EOSTable % TS % minValues(1), EOSTable % TS % maxValues(1),&
         EOSTable % TS % nPoints(1), EOSTable % TS % States(1) % Values)
  CALL MakeLogGrid( EOSTable % TS % minValues(2), EOSTable % TS % maxValues(2),&
         EOSTable % TS % nPoints(2), EOSTable % TS % States(2) % Values)
  CALL MakeLinearGrid( EOSTable % TS % minValues(3), EOSTable % TS % maxValues(3),&
         EOSTable % TS % nPoints(3), EOSTable % TS % States(3) % Values)

  EOSTable % DV % Names(1:13) = (/'Pressure                        ', &
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
                                  'Heavy Charge Number             ', &
                                  'Heavy Binding Energy            '/)

  EOSTable % DV % Units(1:13) = (/'Dynes per cm^2                  ', &
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
                                  'Heavy Mass Number               ', &
                                  'MeV                             '/)

  WRITE (*,*) "Dimensions of rho, T, Ye"
  
  DO i = 1, SIZE( EOSTable % DV % Variables )
    WRITE (*,*) SHAPE( EOSTable % DV % Variables(i) % Values )
  END DO

  DO i = 1, SIZE( EOSTable % DV % Variables )
    EOSTable % DV % Variables(i) % Values = i 
  END DO

  EOSTable % MD % IDTag = 'wl-EOS-LS220-20-40-100, Today'
  EOSTable % MD % TableResolution = '20 pts/dec rho, 40 pts/dec, delta ye = .01'
  EOSTable % MD % NucEOSLink = 'Nuc EOS Paper Link'
  EOSTable % MD % LeptonEOSLink = 'Lepton EOS Paper Link'
  EOSTable % MD % SourceLink = 'Table Source Link'
  EOSTable % MD % WLRevision = 'Rev'
  EOSTable % MD % TableLink = 'Table Storage Link'


  DO j = 1,3
    WRITE(*,*) TRIM( EOSTable % TS % Names(j) )
    WRITE(*,*) EOSTable % TS % nPoints(j)
    WRITE(*,*) EOSTable % TS % minValues(j), EOSTable % TS % maxValues(j)
    WRITE(*,*) EOSTable % TS % States(j) % Values(:)
  END DO

  DO j = 1,13
    WRITE(*,*)
    WRITE(*,*) TRIM( EOSTable % DV % Names(j) ) , j
    WRITE(*,*)
    WRITE(*,*) EOSTable % DV % Variables(j) % Values(:,:,:)
  END DO

  WRITE(*,*), "Table IDTag", EOSTable % MD % IDTag 
  WRITE(*,*), "Table Rez", EOSTable % MD % TableResolution 
  WRITE(*,*), "Table Nuc EOS Paper", EOSTable % MD % NucEOSLink 
  WRITE(*,*), "Table Lepton EOS Paper", EOSTable % MD % LeptonEOSLink 
  WRITE(*,*), "Table Source Link", EOSTable % MD % SourceLink 
  WRITE(*,*), "WeakLib Revision", EOSTable % MD % WLRevision
  WRITE(*,*), "WeakLib Table Link", EOSTable % MD % TableLink


  CALL InitializeHDF( )

  CALL WriteEquationOfStateTableHDF( EOSTable )

  CALL DeAllocateEquationOfStateTable( EOSTable )
  !CALL FinalizeHDF( )
  
  WRITE (*,*) "HDF write successful"

END PROGRAM wlWriteEquationOfStateTest
