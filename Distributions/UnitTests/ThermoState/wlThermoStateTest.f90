PROGRAM wlThermoStateTest

  USE wlKindModule, ONLY: dp
  USE wlGridModule, ONLY: &
    MakeLinearGrid, MakeLogGrid
  USE wlThermoStateModule, ONLY: &
    ThermoStateType, AllocateThermoState, DeAllocateThermoState
  USE wlIOModuleHDF, ONLY: &
    InitializeHDF, OpenFileHDF, OpenGroupHDF, &
    WriteThermoStateHDF, ReadThermoStateHDF, &
    CloseGroupHDF, CloseFileHDF, FinalizeHDF, &
    ReadDimensionsHDF, ReadThermoStateHDF

  USE HDF5

  IMPLICIT NONE

  INTEGER               :: iTS, iPt
  INTEGER,    PARAMETER :: iD = 1, iT = 2, iY = 3
  INTEGER, DIMENSION(3) :: nPoints
  TYPE(ThermoStateType) :: TS
  INTEGER(HID_T)        :: file_id, group_id

  ! --- Create Grid of Thermodynamic States ---

  nPoints = (/ 7, 8, 9 /)

  CALL AllocateThermoState( TS, nPoints )

  TS % nPoints(1:3) = nPoints(1:3)

  TS % Names(1:3) &
    = (/ 'Density                         ',&
         'Temperature                     ',&
         'Electron Fraction               ' /)

  TS % Units(1:3) &
    = (/ 'Grams per cm^3                  ', &
         'K                               ', &
         '                                ' /)

  TS % minValues(1:3) &
    =  (/ 1.0d06, 0.1d00, 1.0d-02 /)
  TS % maxValues(1:3) &
    =  (/ 1.0d15, 1.0d02, 6.1d-01 /)

  CALL MakeLogGrid &
         ( TS % minValues(iD), TS % maxValues(iD), TS % nPoints(iD), &
           TS % States(iD) % Values )
  CALL MakeLogGrid &
         ( TS % minValues(iT), TS % maxValues(iT), TS % nPoints(iT), &
           TS % States(iT) % Values )
  CALL MakeLinearGrid &
         ( TS % minValues(iY), TS % maxValues(iY), TS % nPoints(iY), &
           TS % States(iY) % Values )

  ! --- Write Thermodynamic States to HDF5 File ---

  CALL InitializeHDF( )
  CALL OpenFileHDF( "ThermoStateFile.h5", .true., file_id )
  CALL OpenGroupHDF( "ThermoState", .true., file_id, group_id )
  CALL WriteThermoStateHDF( TS, group_id )
  CALL CloseGroupHDF( group_id )
  CALL CloseFileHDF( file_id )
 
  CALL DeAllocateThermoState( TS )

  ! --- Read Thermodynamic States from File ---

  nPoints = 0

  CALL OpenFileHDF( "ThermoStateFile.h5", .false., file_id )
  CALL OpenGroupHDF( "ThermoState", .false., file_id, group_id )
  CALL ReadDimensionsHDF( nPoints, group_id )

  CALL AllocateThermoState( TS, nPoints )
  CALL ReadThermoStateHDF( TS, file_id )

  CALL CloseFileHDF( file_id )

  ! --- Write Contents of Thermodynamic States ---

  DO iTS = 1, 3
    WRITE(*,*)
    WRITE(*,'(A4,A21,I1,A3,A)') &
      ' ', 'Independent Variable ', iTS, ' = ', TRIM( TS % Names(iTS) )
    WRITE(*,'(A6,A8,I1,A4,I4.4)') &
      ' ', 'nPoints(', iTS, ') = ', TS % nPoints(iTS)
    WRITE(*,'(A6,A12,ES10.4E2,A3,ES10.4E2)') &
      ' ', 'Min / Max = ', TS % minValues(iTS),' / ', TS % maxValues(iTS)
    WRITE(*,'(A6,A7)') ' ', 'Values:'
    DO iPt = 1, nPoints(iTS)
      WRITE(*,'(A8,A6,I1,A4,ES10.4E2,A1,A)') &
        ' ', 'Value(', iPt,') = ', TS % States(iTS) % Values(iPt), ' ', TRIM( TS % Units(iTS) )
    END DO
  END DO
  WRITE(*,*)

  CALL DeAllocateThermoState( TS )

  CALL FinalizeHDF( )

END PROGRAM wlThermoStateTest
