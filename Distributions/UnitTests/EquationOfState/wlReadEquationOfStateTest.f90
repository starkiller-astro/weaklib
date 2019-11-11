PROGRAM wlReadEquationOfStateTest
 
  USE HDF5
  USE wlEquationOfStateTableModule
  USE wlIOModuleHDF
  USE wlEOSIOModuleHDF

  implicit none

  INTEGER                        :: i
  TYPE(EquationOfStateTableType) :: EOSTable

  CALL InitializeHDF( )

  CALL ReadEquationOfStateTableHDF( EOSTable, "wl-EOS-SFHo-15-25-50.h5" )

  WRITE (*,*) 'EOSTable % nPoints:', EOSTable % nPoints
  WRITE (*,*) 'EOSTable % nVariables:', EOSTable % nVariables

  DO i = 1, SIZE( EOSTable % DV % Variables )
    EOSTable % DV % Variables(i) % Values = i
  END DO

  DO i = 1,3
    WRITE(*,*) TRIM( EOSTable % TS % Names(i) )
    WRITE(*,*) EOSTable % TS % nPoints(i)
    WRITE(*,*) EOSTable % TS % minValues(i), EOSTable % TS % maxValues(i)
!!$    WRITE(*,*) EOSTable % TS % States(i) % Values(:)
  END DO

  DO i = 1, SIZE( EOSTable % DV % Variables )
    WRITE(*,*)
    WRITE(*,*) TRIM( EOSTable % DV % Names(i) ) , i
    WRITE(*,*)
!!$    WRITE(*,*) EOSTable % DV % Variables(i) % Values(:,:,:)
  END DO

  WRITE(*,*) "Metadata!" 
  WRITE(*,*) "Table IDTag: ", EOSTable % MD % IDTag
  WRITE(*,*) "Table Rez: ", EOSTable % MD % TableResolution
  WRITE(*,*) "Table Nuc EOS Paper: ", EOSTable % MD % NucEOSLink
  WRITE(*,*) "Table Lepton EOS Paper: ", EOSTable % MD % LeptonEOSLink
  WRITE(*,*) "Table Source Link: ", EOSTable % MD % SourceLink
  WRITE(*,*) "WeakLib Revision: ", EOSTable % MD % WLRevision
  WRITE(*,*) "WeakLib Table Link: ", EOSTable % MD % TableLink

  CALL DeAllocateEquationOfStateTable( EOSTable )

  CALL FinalizeHDF( )

END PROGRAM wlReadEquationOfStateTest
