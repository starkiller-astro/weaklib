MODULE wlOpacityTableIOModuleHDF
!-----------------------------------------------------------------------
!
!    File:         wlOpacityTableIOModuleHDF.f90
!    Module:       wlOpacityTableIOModuleHDF
!    Type:         Module w/ Subroutines
!    Author:       R. Chu, Dept. Phys. & Astronomy
!                  U. Tennesee, Knoxville
!
!    Created:      3/22/16
!    WeakLib ver:  
!
!    Purpose:
!      Subroutines needed for reading, printing OpacityTable
!
!    CONTAINS:
!       DescribeOpacityTable
!
!    Modules used:
!       wlOpacityTableModule, ONLY: OpacityTableType, OpacityTypeA
!       wlEOSIOModuleHDF, ONLY: DescribeEquationOfStateTable
!       wlGridModule, ONLY: EnergyGridType
!       wlKindModule, ONLY:dp
!       HDF5
!       wlIOModuleHDF
!       wlEquationOfStateTableModule
!
!-----------------------------------------------------------------------
!  NOTE: Only Type A interaction applied. Type B and Type C interaction 
!        needs to be added for future use.
!-----------------------------------------------------------------------

  USE wlKindModule, ONLY:dp
  USE wlEnergyGridModule, ONLY: &
    EnergyGridType
  USE wlOpacityTableModule, ONLY: OpacityTableType, AllocateOpacityTable
  USE wlOpacityFieldsModule, ONLY: OpacityTypeA
  USE wlIOModuleHDF
  USE wlEquationOfStateTableModule

  USE HDF5

  IMPLICIT NONE
  PRIVATE

  INTEGER :: hdferr

  PUBLIC WriteOpacityTableHDF
  PUBLIC ReadOpacityTableHDF

CONTAINS
 
  SUBROUTINE WriteOpacityTableHDF( OpacityTable )
 
    TYPE(OpacityTableType), INTENT(inout) :: OpacityTable

    INTEGER(HID_T)                                :: file_id
    INTEGER(HID_T)                                :: group_id

    CHARACTER(LEN=32), DIMENSION(1)             :: tempString
    INTEGER, DIMENSION(1)                       :: tempInteger
    INTEGER(HSIZE_T), DIMENSION(1)              :: datasize1d
   
    CALL OpenFileHDF( "OpacityTable_NS.h5", .true., file_id )

    datasize1d(1) = 1
    tempInteger(1) = OpacityTable % nOpacitiesA
    CALL Write1dHDF_integer&
         ( "nOpacitiesA", tempInteger, file_id, datasize1d )

    tempInteger(1) = OpacityTable % nOpacitiesB 
    CALL Write1dHDF_integer&
         ( "nOpacitiesB", tempInteger, file_id, datasize1d )
  
    tempInteger(1) = OpacityTable % nMomentsB     
    CALL Write1dHDF_integer&
         ( "nMomentsB", tempInteger, file_id, datasize1d )

    tempInteger(1) = OpacityTable % nOpacitiesC     
    CALL Write1dHDF_integer&
         ( "nOpacitiesC", tempInteger, file_id, datasize1d )

    tempInteger(1) = OpacityTable % nMomentsC   
    CALL Write1dHDF_integer&
         ( "nMomentsC", tempInteger, file_id, datasize1d )

    tempInteger(1) = OpacityTable % nPointsE  
    CALL Write1dHDF_integer&
         ( "nPointsE", tempInteger, file_id, datasize1d )

    datasize1d = 3
    CALL Write1dHDF_integer&
         ( "nPointsTS", OpacityTable % nPointsTS, file_id, datasize1d )

    CALL OpenGroupHDF( "thermEmAb", .true., file_id, group_id )
    CALL WriteOpacityTableTypeAHDF( OpacityTable % thermEmAb, group_id )
    CALL CloseGroupHDF( group_id )

    CALL OpenGroupHDF( "EnergyGrid", .true., file_id, group_id )
    CALL WriteEnergyGridHDF( OpacityTable % EnergyGrid, group_id )
    CALL CloseGroupHDF( group_id )

    CALL OpenGroupHDF( "EOSTable", .true., file_id, group_id )
    CALL WriteEOSTableHDF( OpacityTable % EOSTable, file_id, group_id )
    CALL CloseGroupHDF( group_id )
     
    CALL CloseFileHDF( file_id )

  END SUBROUTINE WriteOpacityTableHDF


  SUBROUTINE WriteEOSTableHDF( EOSTable, file_id, group_id )

    TYPE(EquationOfStateTableType), INTENT(in)    :: EOSTable
    INTEGER(HID_T), INTENT(in)                    :: file_id
    INTEGER(HID_T), INTENT(in)                    :: group_id

    INTEGER(HID_T)                                :: subgroup_id

    CALL OpenGroupHDF( "ThermoState", .true., group_id, subgroup_id )
    CALL WriteThermoStateHDF( EOSTable % TS, subgroup_id )
    CALL CloseGroupHDF( subgroup_id )

    CALL OpenGroupHDF( "DependentVariables", .true., group_id, subgroup_id )
    CALL WriteDependentVariablesHDF( EOSTable % DV, subgroup_id )
    CALL CloseGroupHDF( subgroup_id )

  END SUBROUTINE WriteEOSTableHDF


  SUBROUTINE WriteEnergyGridHDF( EnergyGrid, group_id )

    TYPE(EnergyGridType), INTENT(in)           :: EnergyGrid
    INTEGER(HID_T), INTENT(in)                 :: group_id

    INTEGER(HSIZE_T), DIMENSION(1)              :: datasize1d
    INTEGER                                     :: i
    INTEGER, DIMENSION(1)                       :: buffer

    CHARACTER(LEN=32), DIMENSION(1)             :: tempString
    INTEGER, DIMENSION(1)                       :: tempInteger          

    datasize1d(1) = 1

    tempString(1) = EnergyGrid % Name
    CALL Write1dHDF_string( "Name", tempString, &
                             group_id, datasize1d )
    
    tempString(1) = EnergyGrid % Unit
    CALL Write1dHDF_string( "Unit", tempString, &
                            group_id, datasize1d )

    tempInteger(1) = EnergyGrid % nPoints  
    CALL Write1dHDF_integer( "nPoints", tempInteger, &
                            group_id, datasize1d )
   
    tempInteger(1) = EnergyGrid % LogInterp 
    CALL Write1dHDF_integer( "LogInterp", tempInteger, &
                             group_id, datasize1d )
   
    datasize1d(1) = EnergyGrid % nPoints
    CALL Write1dHDF_double( "Values", EnergyGrid % Values(:), &
                              group_id, datasize1d )
  
  END SUBROUTINE WriteEnergyGridHDF


  SUBROUTINE WriteOpacityTableTypeAHDF( thermEmAb, group_id )

    TYPE(OpacityTypeA), INTENT(in)              :: thermEmAb
    INTEGER(HID_T), INTENT(in)                  :: group_id

    INTEGER(HSIZE_T)                            :: datasize1d
    INTEGER(HSIZE_T), DIMENSION(4)              :: datasize4d
    INTEGER                                     :: i
    INTEGER, DIMENSION(1)                       :: buffer


    CHARACTER(LEN=32), DIMENSION(1)             :: tempString
    INTEGER, DIMENSION(1)                       :: tempInteger
    INTEGER(HSIZE_T), DIMENSION(1)              :: datasize1dtemp

    INTEGER(HID_T)                                :: subgroup_id


    datasize1dtemp(1) = 1
    tempInteger(1) = thermEmAb % nOpacities
    CALL Write1dHDF_integer&
         ( "nOpacities", tempInteger, group_id, datasize1dtemp )

    datasize1dtemp(1) = 4
    CALL Write1dHDF_integer&
         ( "nPoints", thermEmAb % nPoints, group_id, datasize1dtemp )

    datasize1dtemp(1) = thermEmAb % nOpacities
    CALL Write1dHDF_string&
         ( "Names", thermEmAb % Names, group_id, datasize1dtemp ) 

    CALL Write1dHDF_string&
         ( "Species", thermEmAb % Species, group_id, datasize1dtemp ) 

    CALL Write1dHDF_string&
         ( "Units", thermEmAb % Units, group_id, datasize1dtemp ) 

    datasize1d = thermEmAb % nOpacities 
    Write(*,*) 'datasize1d before write Ab in writing loop is ', datasize1d
    CALL OpenGroupHDF( "Absorptivity", .true., group_id, subgroup_id )
      DO i = 1, datasize1d
        datasize4d = SHAPE( thermEmAb % Absorptivity(i) % Values )
       CALL Write4dHDF_double( "Values", thermEmAb % Absorptivity(i) % Values(:,:,:,:), &
                              subgroup_id, datasize4d )
      END DO
    CALL CloseGroupHDF( subgroup_id )
 
  END SUBROUTINE WriteOpacityTableTypeAHDF


  SUBROUTINE Write4dHDF_double &
              ( name, values, group_id, datasize, desc_option, unit_option )

    CHARACTER(*), INTENT(in)                    :: name
    CHARACTER(*), INTENT(in), OPTIONAL          :: unit_option
    CHARACTER(*), INTENT(in), OPTIONAL          :: desc_option
    INTEGER(HID_T)                              :: group_id
    INTEGER(HSIZE_T), DIMENSION(4), INTENT(in)  :: datasize
    REAL(dp), DIMENSION(:,:,:,:), INTENT(in)      :: values
   
    INTEGER(HID_T)                              :: dataset_id
    INTEGER(HID_T)                              :: dataspace_id
    INTEGER(HID_T)                              :: atype_id
    INTEGER(HID_T)                              :: attr_id
    INTEGER(SIZE_T)                             :: attr_len
    INTEGER(HSIZE_T), DIMENSION(1)              :: adims = (/1/)
  
    
    CALL h5screate_simple_f( 4, datasize, dataspace_id, hdferr )

    CALL h5dcreate_f( group_id, name, H5T_NATIVE_DOUBLE, &
           dataspace_id, dataset_id, hdferr )

    CALL h5dwrite_f( dataset_id, H5T_NATIVE_DOUBLE, &
           values, datasize, hdferr )

    CALL h5sclose_f( dataspace_id, hdferr ) 

    CALL h5dclose_f( dataset_id, hdferr )

  END SUBROUTINE Write4dHDF_double


  SUBROUTINE ReadOpacityTableHDF( OpacityTable, FileName )
 
    TYPE(OpacityTableType), INTENT(inout)   :: OpacityTable
    CHARACTER(len=*), INTENT(in)            :: FileName

!    INTEGER, DIMENSION(3)                         :: nPoints
    INTEGER                                       :: nPointsE
    INTEGER                                       :: nOpacA
    INTEGER                                       :: nOpacB, nMomB
    INTEGER                                       :: nOpacC, nMomC
!    INTEGER                                       :: nVariables
    INTEGER(HID_T)                                :: file_id
    INTEGER(HID_T)                                :: group_id
    INTEGER(HID_T)                                :: subgroup_id
    INTEGER(HSIZE_T), DIMENSION(1)                :: datasize1d
    INTEGER, DIMENSION(1)                         :: buffer

    CALL OpenFileHDF( FileName, .false., file_id )

    datasize1d(1) = 1
    CALL Read1dHDF_integer( "nOpacitiesA", buffer, file_id, datasize1d )
    nOpacA = buffer(1)   

    CALL Read1dHDF_integer( "nOpacitiesB", buffer, file_id, datasize1d )
    nOpacB = buffer(1)

    CALL Read1dHDF_integer( "nMomentsB", buffer, file_id, datasize1d )
    nMomB = buffer(1)

    CALL Read1dHDF_integer( "nOpacitiesC", buffer, file_id, datasize1d )
    nOpacC = buffer(1)

    CALL Read1dHDF_integer( "nMomentsC", buffer, file_id, datasize1d )
    nMomC = buffer(1)

    CALL Read1dHDF_integer( "nPointsE", buffer, file_id, datasize1d )
    nPointsE = buffer(1)

     
 !   CALL AllocateOpacityTable &
!               ( OpacityTable, nOpacA, nOpacB, nMomB, nOpacC, nMomC, nPointsE )
!    CALL OpenGroupHDF( "EOSTable", .false., file_id, group_id )
   
!    Call OpenGroupHDF( "DependentVariables", .false., group_id, subgroup_id )     
!    CALL ReadDimensionsHDF( nPoints, subgroup_id )
!!$    Write(*,*),'EOS nPoints', nPoints
!!$    CALL ReadNumberVariablesHDF( nVariables, subgroup_id )
!!$    Write(*,*),'EOS nVariables', nVariables
!!$    CALL CloseGroupHDF( subgroup_id )
!!$   
!!$    CALL CloseGroupHDF( group_id )
!!$  
!!$    CALL OpenGroupHDF( "EnergyGrid", .false., file_id, group_id )
!!$    CALL ReadEnergyPointsHDF( nPointsE, group_id)
!!$    Write(*,*),'EngeryGrid nPointsE', nPointsE
!!$    CALL CloseGroupHDF( group_id )
!!$
!!$    CALL OpenGroupHDF( "EcapEm", .false., file_id, group_id )
!!$  !  CALL ReadTableDimensionsHDF( nSpeciesA, group_id )
!!$    WRITE(*,*),'EcapEm number', nSpeciesA
!!$    CALL CloseGroupHDF( group_id )
!!$    
!!$    CALL AllocateEmptyOpacityTable &
!!$           ( OpacityTable, nSpeciesA, nPointsE, nPoints, nVariables )
!!$   
!!$    CALL OpenGroupHDF( "EOSTable", .false., file_id, group_id )
!!$    CALL ReadThermoStateHDF( OpacityTable % EOSTable % TS, group_id )
!!$    CALL ReadDependentVariablesHDF( OpacityTable % EOSTable % DV, group_id )
!!$    CALL CloseGroupHDF( group_id )
!!$ 
!!$    CALL ReadOpacityTypeAHDF( OpacityTable % ECAPEM , file_id )
!!$
!!$    CALL ReadEnergyGridHDF( OpacityTable % EnergyGrid, file_id  )
!!$
!!$    CALL DescribeOpacityTable( OpacityTable )

    CALL CloseFileHDF( file_id )

  END SUBROUTINE ReadOpacityTableHDF


  SUBROUTINE ReadOpacityTypeAHDF( ECAPEM, file_id )

    TYPE(OpacityTypeA), DIMENSION(:),INTENT(inout)   :: ECAPEM
    INTEGER(HID_T), INTENT(in)                       :: file_id

!!$    INTEGER(HID_T)                              :: group_id
!!$    INTEGER(HSIZE_T), DIMENSION(4)              :: datasize4d
!!$    INTEGER                                     :: i
!!$    INTEGER, DIMENSION(1)                       :: buffer
!!$
!!$    CALL OpenGroupHDF( "EcapEm", .false., file_id, group_id )  
!!$      
!!$    DO i = 1, 1
!!$    CALL Read4dHDF_double( "Valuse ", ECAPEM(i) % Values, &
!!$                              group_id, datasize4d )
!!$    END DO

  END SUBROUTINE ReadOpacityTypeAHDF


  SUBROUTINE Read4dHDF_double( name, values, group_id, datasize )

    CHARACTER(*), INTENT(in)                    :: name
    INTEGER(HID_T)                              :: group_id
    INTEGER(HSIZE_T), DIMENSION(4), INTENT(in)  :: datasize
    REAL(dp), DIMENSION(:,:,:,:), INTENT(out)   :: values
   
!!$    INTEGER(HID_T)                               :: dataset_id
!!$ 
!!$    CALL h5dopen_f( group_id, name, dataset_id, hdferr )
!!$    CALL h5dread_f( dataset_id, H5T_NATIVE_DOUBLE, values, datasize, hdferr )
!!$    CALL h5dclose_f( dataset_id, hdferr )

  END SUBROUTINE Read4dHDF_double


  SUBROUTINE ReadEnergyGridHDF( EnergyGrid, file_id )

    TYPE(EnergyGridType), INTENT(inout)         :: EnergyGrid
    INTEGER(HID_T), INTENT(in)                  :: file_id

!!$    INTEGER(HID_T)                              :: group_id
!!$    INTEGER(HSIZE_T), DIMENSION(1)              :: datasize1d
!!$
!!$    CALL OpenGroupHDF( "EnergyGrid", .false., file_id, group_id )
!!$
!!$ 
!!$!    CALL Read1HDF_string( "Names", EnergyGrid % Names, &
!!$!                              group_id )
!!$!    CALL Read1HDF_string( "Units", EnergyGrid % Units, &
!!$!                              group_id )
!!$    CALL Read1dHDF_double( "Values", EnergyGrid % Values(:), &
!!$                              group_id, datasize1d )
!!$
!!$    CALL CloseGroupHDF( group_id )

  END SUBROUTINE ReadEnergyGridHDF


  SUBROUTINE ReadEnergyPointsHDF ( nPointsE, group_id )

    INTEGER(HID_T), INTENT(in)                  :: group_id
    INTEGER, INTENT(inout)                      :: nPointsE

!!$    INTEGER, DIMENSION(1)                       :: nVarTemp
!!$    INTEGER(HSIZE_T), DIMENSION(1)              :: datasize1d
!!$
!!$    datasize1d(1) = 1
!!$    CALL Read1dHDF_integer( "nPointsE", nVarTemp(:), group_id, datasize1d )
!!$    nPointsE = nVarTemp(1)

  END SUBROUTINE ReadEnergyPointsHDF


  SUBROUTINE Read1HDF_string( name, values, group_id )
  
    CHARACTER(*), INTENT(in)                      :: name
    INTEGER(HID_T)                                :: group_id
    CHARACTER(len=*), INTENT(inout)               :: values

!!$    INTEGER(HSIZE_T), DIMENSION(1)                :: datasize
!!$
!!$!    CALL Read1dHDF_string( name, tempvalues, group_id, datasize )

  END SUBROUTINE Read1HDF_string


END MODULE wlOpacityTableIOModuleHDF
