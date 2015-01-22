PROGRAM wlThermoStateTest

  USE wlThermoStateModule
  USE wlKindModule, ONLY: dp
  USE wlGridModule, ONLY: MakeLinearGrid, MakeLogGrid
  USE HDF5
  USE wlIOModuleHDF, ONLY: InitializeHDF, OpenFileHDF, OpenGroupHDF,         &
                           WriteThermoStateHDF, ReadThermoStateHDF,          & 
                           CloseGroupHDF, CloseFileHDF, FinalizeHDF,         &
                           ReadDimensionsHDF, ReadThermoStateHDF

  implicit none


  INTEGER, DIMENSION(3) :: npts
  INTEGER, DIMENSION(3) :: npts2
  INTEGER :: nvar
  INTEGER :: j   
  TYPE(ThermoStateType) :: ThermoState
  TYPE(ThermoStateType) :: ThermoState2
  INTEGER(HID_T) :: file_id
  INTEGER(HID_T) :: group_id

  npts = (/10,31,61/)


  CALL AllocateThermoState( ThermoState, npts )

  ThermoState % nPoints(1:3) = npts(1:3)
  ThermoState % Names(1:3) = (/'Density                         ',&
                               'Temperature                     ',&
                               'Electron Fraction               '/)
    
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
  CALL ReadThermoStateHDF( ThermoState2, file_id )
  !CALL OpenGroupHDF( "ThermoState", .false., file_id, group_id )
  !CALL ReadDimensionsHDF( npts, group_id )
  !CALL AllocateThermoState( ThermoState2, npts )

  !ThermoState2 % nPoints(1:3) = npts(1:3)

  !CALL ReadThermoStateHDF( ThermoState2, group_id )
  !CALL CloseGroupHDF( group_id )
  CALL CloseFileHDF( file_id )

  DO j = 1,3
    WRITE(*,*) TRIM( ThermoState2 % Names(j) )
    WRITE(*,*) ThermoState2 % nPoints(j)
    WRITE(*,*) ThermoState2 % minValues(j), ThermoState2 % maxValues(j)
    WRITE(*,*) ThermoState2 % States(j) % Values(:)
  END DO

  CALL DeAllocateThermoState( ThermoState2 )

  CALL FinalizeHDF( )

END PROGRAM wlThermoStateTest
