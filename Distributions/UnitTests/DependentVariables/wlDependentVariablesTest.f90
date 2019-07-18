PROGRAM wlDependentVariablesTest
 
  USE wlDependentVariablesModule, ONLY: &
    DependentVariablesType, &
    AllocateDependentVariables, &
    DeAllocateDependentVariables
  USE wlIOModuleHDF, ONLY: &
    InitializeHDF, OpenFileHDF, OpenGroupHDF, &
    WriteDependentVariablesHDF, CloseGroupHDF, & 
    ReadDependentVariablesHDF, CloseFileHDF, &
    FinalizeHDF, ReadDimensionsHDF, &
    ReadDependentVariablesHDF, ReadNumberVariablesHDF

  USE HDF5

  IMPLICIT NONE

  INTEGER                      :: iVar, nVariables
  INTEGER, DIMENSION(3)        :: nPoints
  TYPE(DependentVariablesType) :: DV
  INTEGER(HID_T)               :: file_id, group_id

  ! --- Create Table Structures for Dependent Variables ---

  CALL AllocateDependentVariables &
         ( DV, nPoints = (/ 10, 10, 10 /) , nVariables = 13 )

  DV % Names(1:DV % nVariables) &
    = (/ 'Pressure                        ', &
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
         'Heavy Binding Energy            ' /)

  DV % Units(1:DV % nVariables) &
    = (/ 'Dynes per cm^2                  ', &
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
         'MeV                             ' /)

  WRITE(*,*)
  DO iVar = 1, SIZE( DV % Variables )
    WRITE (*,'(A4,A32,A10,3I5.3)') &
      ' ', TRIM( DV % Names(iVar) ), ', Shape = ', SHAPE( DV % Variables(iVar) % Values )
  END DO
  WRITE(*,*)

  ! --- Initialize Table Entries with Mock Values --- 

  DO iVar = 1, SIZE( DV % Variables )
    DV % Variables(iVar) % Values = REAL( iVar )
  END DO

  ! --- Write Table to HDF5 File --- 

  CALL InitializeHDF( )
  CALL OpenFileHDF( "DependentVariablesFile.h5", .true., file_id )
  CALL OpenGroupHDF( "DependentVariables", .true., file_id, group_id )
  CALL WriteDependentVariablesHDF( DV, group_id )
  CALL CloseGroupHDF( group_id )
  CALL CloseFileHDF( file_id )

  CALL DeAllocateDependentVariables( DV )

  IF( ALLOCATED( DV % Variables ) )THEN
    WRITE(*,*)
    WRITE(*,'(A4,A)') &
      ' ', 'ERROR: Dependent Variables Allocated After Deallocation'
    WRITE(*,*)
    STOP
  END IF

  nPoints = 0; nVariables = 0

  ! --- Read Dependent Variables from HDF5 File ---

  CALL OpenFileHDF( "DependentVariablesFile.h5", .false., file_id )
  CALL OpenGroupHDF( "DependentVariables", .false., file_id, group_id )

  CALL ReadDimensionsHDF( nPoints, group_id )
  CALL ReadNumberVariablesHDF( nVariables, group_id )

  CALL AllocateDependentVariables &
         ( DV, nPoints = nPoints, nVariables = nVariables )
  CALL ReadDependentVariablesHDF( DV, file_id )

  CALL CloseFileHDF( file_id )

  ! --- Write Contents of Dependent Variables ---

  DO iVar = 1, DV % nVariables
    WRITE(*,*)
    WRITE(*,'(A4,A19,I2.2,A3,A)') &
      ' ', 'Dependent Variable ', iVar, ' = ', TRIM( DV % Names(iVar) )
    WRITE(*,'(A6,A9,A)') &
      ' ', 'Units  = ', TRIM( DV % Units(iVar) )
    WRITE(*,'(A6,A8,3I4.3)') &
      ' ', 'Shape  =', SHAPE( DV % Variables(iVar) % Values )
    WRITE(*,'(A6,A19,ES10.4,A3,ES10.4)') &
      ' ', 'Min / Max Values = ', MINVAL( DV % Variables(iVar) % Values ), &
                           ' / ', MAXVAL( DV % Variables(iVar) % Values )
  END DO
  WRITE(*,*)

  CALL DeAllocateDependentVariables( DV )

  CALL FinalizeHDF( )
 
END PROGRAM wlDependentVariablesTest
