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

  INTEGER :: i
  INTEGER, DIMENSION(3) :: npts
  INTEGER, DIMENSION(3) :: npts2
  TYPE(ThermoStateType) :: ThermoState
  TYPE(ThermoStateType) :: ThermoState2
  INTEGER(HID_T) :: file_id
  INTEGER(HID_T) :: group_id

  npts = (/10,31,61/)

  CALL AllocateThermoState( ThermoState, npts )

  ThermoState % nPoints(1:3) = npts(1:3)
  ThermoState % Names(1:3) = (/ 'Density                         ',&
                                'Temperature                     ',&
                                'Electron Fraction               ' /)

  ThermoState % Units(1:3) = (/ 'Grams per cm^3                  ', &
                                'K                               ', &
                                '                                ' /)
    
  ThermoState % minValues(1:3) =  (/1.0d06,0.1d00,1.0d-02/)
  ThermoState % maxValues(1:3) =  (/1.0d15,1.0d02,6.1d-01/)

  CALL MakeLogGrid( ThermoState % minValues(1), ThermoState % maxValues(1),&
         ThermoState % nPoints(1), ThermoState % States(1) % Values)
  CALL MakeLogGrid( ThermoState % minValues(2), ThermoState % maxValues(2),&
         ThermoState % nPoints(2), ThermoState % States(2) % Values)
  CALL MakeLinearGrid( ThermoState % minValues(3), ThermoState % maxValues(3),&
         ThermoState % nPoints(3), ThermoState % States(3) % Values)

  CALL InitializeHDF( )
  CALL OpenFileHDF( "ThermoStateFile.h5", .true., file_id )
  CALL OpenGroupHDF( "ThermoState", .true., file_id, group_id )
  CALL WriteThermoStateHDF( ThermoState, group_id )
  CALL CloseGroupHDF( group_id )
  CALL CloseFileHDF( file_id )
 
  CALL DeAllocateThermoState( ThermoState )

  CALL OpenFileHDF( "ThermoStateFile.h5", .false., file_id )
  CALL OpenGroupHDF( "ThermoState", .false., file_id, group_id )
  CALL ReadDimensionsHDF( npts2, group_id )
  CALL AllocateThermoState( ThermoState2, npts2 )
  ThermoState2 % nPoints(1:3) = npts(1:3)
  ThermoState2 % Names(1:3) = (/'Density                         ',&
                               'Temperature                     ',&
                               'Electron Fraction               '/)
    
  ThermoState2 % minValues(1:3) =  (/1.0d06,0.1d00,1.0d-02/)
  ThermoState2 % maxValues(1:3) =  (/1.0d15,1.0d02,6.1d-01/)
  ThermoState2 % nPoints(1:3) = npts2(1:3)
  CALL ReadThermoStateHDF( ThermoState2, file_id )
  CALL CloseFileHDF( file_id )

  DO i = 1,3
    WRITE(*,*) TRIM( ThermoState2 % Names(i) )
    WRITE(*,*) ThermoState2 % nPoints(i)
    WRITE(*,*) ThermoState2 % minValues(i), ThermoState2 % maxValues(i)
    WRITE(*,*) ThermoState2 % States(i) % Values(:)
  END DO

  CALL DeAllocateThermoState( ThermoState2 )

  CALL FinalizeHDF( )

END PROGRAM wlThermoStateTest
