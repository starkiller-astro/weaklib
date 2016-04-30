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
  USE wlOpacityTableModule, ONLY: OpacityTableType
  USE wlOpacityFieldsModule, ONLY: OpacityTypeA
  USE wlEOSIOModuleHDF, ONLY: DescribeEquationOfStateTable
  USE wlIOModuleHDF
  USE wlEquationOfStateTableModule

  USE HDF5

  implicit none
  PRIVATE
  INTEGER                                     :: hdferr

  PUBLIC WriteOpacityTableHDF
  PUBLIC ReadOpacityTableHDF
  PUBLIC DescribeOpacityTable
 ! PUBLIC BroadcastOpacityTableParallel
  PUBLIC DescribeOpacityTypeA
  PUBLIC WriteOpacityTableTypeAHDF
  PUBLIC WriteEnergyGridHDF

CONTAINS
 
  SUBROUTINE WriteOpacityTableHDF( OpacityTable )
 
    TYPE(OpacityTableType), INTENT(inout) :: OpacityTable

!!$    INTEGER(HID_T)                                :: file_id
!!$    INTEGER(HID_T)                                :: group_id
!!$    
!!$    CALL OpenFileHDF( "OpacityTable.h5", .true., file_id )
!!$
!!$    CALL OpenGroupHDF( "EcapEm", .true., file_id, group_id )
!!$    CALL WriteOpacityTableTypeAHDF( OpacityTable % ECAPEM, group_id )
!!$    CALL CloseGroupHDF( group_id )
!!$
!!$    CALL OpenGroupHDF( "EnergyGrid", .true., file_id, group_id )
!!$    CALL WriteEnergyGridHDF( OpacityTable % EnergyGrid, group_id )
!!$    CALL CloseGroupHDF( group_id )
!!$
!!$    CALL OpenGroupHDF( "EOSTable", .true., file_id, group_id )
!!$    CALL WriteEOSTableHDF( OpacityTable % EOSTable, file_id, group_id )
!!$    CALL CloseGroupHDF( group_id )
!!$     
!!$    CALL CloseFileHDF( file_id )

  END SUBROUTINE WriteOpacityTableHDF


  SUBROUTINE DescribeOpacityTable( OpacityTable )
  
    TYPE(OpacityTableType), INTENT(inout) :: OpacityTable

!!$    CALL DescribeEquationOfStateTable( OpacityTable % EOSTable )
!!$    CALL DescribeEnergyGrid( OpacityTable % EnergyGrid )
!!$    CALL DescribeOpacityTypeA( OpacityTable % ECAPEM )

  END SUBROUTINE DescribeOpacityTable


  SUBROUTINE DescribeOpacityTypeA( ECAPEM )

    TYPE(OpacityTypeA), DIMENSION(:), INTENT(inout) :: ECAPEM

!!$    INTEGER :: i, j, k, l, ii, jj, kk, ll, iii
!!$    
!!$    DO iii = 1, SIZE(ECAPEM)
!!$
!!$      WRITE (UNIT=6,FMT=*) "EcapEm Names =", ECAPEM(iii) % Names
!!$      WRITE (UNIT=6,FMT=*) "EcapEm SIZE =", SIZE(ECAPEM)
!!$
!!$      i=SIZE(ECAPEM(iii) % Values,1)  ! E
!!$      j=SIZE(ECAPEM(iii) % Values,2)  ! rho
!!$      k=SIZE(ECAPEM(iii) % Values,3)  ! T
!!$      l=SIZE(ECAPEM(iii) % Values,4)  ! Ye
!!$  
!!$      WRITE (*,*) "EcapEm % SIZE",i, j, k, l
!!$      DO ii = 20,20
!!$         DO jj = 81, 81
!!$            DO kk = 11, 12
!!$               DO ll = 1, l
!!$                  WRITE (*,*) "EcapEm(", iii, ") % Values",ii, jj,&
!!$                  kk, ll, "=",  ECAPEM(iii) % Values(ii,jj,kk,ll),&
!!$                   "Units:", ECAPEM(iii) % Units
!!$               END DO
!!$            END DO
!!$         END DO
!!$      END DO
!!$
!!$    END DO

  END SUBROUTINE DescribeOpacityTypeA

  SUBROUTINE WriteEnergyGridHDF( EnergyGrid, group_id )

    TYPE(EnergyGridType), INTENT(in)           :: EnergyGrid
    INTEGER(HID_T), INTENT(in)                 :: group_id

!!$    INTEGER(HSIZE_T), DIMENSION(1)              :: datasize1d
!!$    INTEGER                                     :: i
!!$    INTEGER, DIMENSION(1)                       :: buffer
!!$
!!$  WRITE (*,*) "Starting HDF EnergyGrid write "
!!$
!!$    datasize1d(1) = 1
!!$
!!$    CALL Write1HDF_integer( "nPointsE", EnergyGrid % nPointsE, &
!!$                             group_id, datasize1d )
!!$    
!!$    CALL Write1HDF_string( "Names", EnergyGrid % Names, &
!!$                             group_id, datasize1d )
!!$
!!$    CALL Write1HDF_string( "Units", EnergyGrid % Units, &
!!$                             group_id, datasize1d )
!!$   
!!$    datasize1d(1) = EnergyGrid % nPointsE
!!$    CALL Write1dHDF_double( "Values", EnergyGrid % Values(:), &
!!$                              group_id, datasize1d )
  
  END SUBROUTINE WriteEnergyGridHDF

  SUBROUTINE WriteOpacityTableTypeAHDF( ECAPEM, group_id )

    TYPE(OpacityTypeA), DIMENSION(:), INTENT(in)           :: ECAPEM
    INTEGER(HID_T), INTENT(in)                 :: group_id

!!$    INTEGER(HSIZE_T)                            :: datasize1d
!!$    INTEGER(HSIZE_T), DIMENSION(4)              :: datasize4d
!!$    INTEGER                                     :: i
!!$    INTEGER, DIMENSION(1)                       :: buffer
!!$
!!$  WRITE (*,*) "Starting HDF OpacityTypeA write "
!!$    datasize1d = SIZE( ECAPEM )
!!$    DO i = 1, datasize1d
!!$     datasize4d = SHAPE( ECAPEM(i) % Values )
!!$     CALL Write4dHDF_double( "Values", ECAPEM(i) % Values(:,:,:,:), &
!!$                              group_id, datasize4d )
!!$    END DO
  
  END SUBROUTINE WriteOpacityTableTypeAHDF

  SUBROUTINE Write4dHDF_double &
              ( name, values, group_id, datasize, desc_option, unit_option )

    CHARACTER(*), INTENT(in)                    :: name
    CHARACTER(*), INTENT(in), OPTIONAL          :: unit_option
    CHARACTER(*), INTENT(in), OPTIONAL          :: desc_option
    INTEGER(HID_T)                              :: group_id
    INTEGER(HSIZE_T), DIMENSION(4), INTENT(in)  :: datasize
    REAL(dp), DIMENSION(:,:,:,:), INTENT(in)      :: values
   
!!$    INTEGER(HID_T)                              :: dataset_id
!!$    INTEGER(HID_T)                              :: dataspace_id
!!$    INTEGER(HID_T)                              :: atype_id
!!$    INTEGER(HID_T)                              :: attr_id
!!$    INTEGER(SIZE_T)                             :: attr_len
!!$    INTEGER(HSIZE_T), DIMENSION(1)              :: adims = (/1/)
!!$  
!!$    
!!$    CALL h5screate_simple_f( 4, datasize, dataspace_id, hdferr )
!!$
!!$    CALL h5dcreate_f( group_id, name, H5T_NATIVE_DOUBLE, &
!!$           dataspace_id, dataset_id, hdferr )
!!$
!!$    CALL h5dwrite_f( dataset_id, H5T_NATIVE_DOUBLE, &
!!$           values, datasize, hdferr )
!!$
!!$    CALL h5sclose_f( dataspace_id, hdferr ) 
!!$
!!$    CALL h5dclose_f( dataset_id, hdferr )

  END SUBROUTINE Write4dHDF_double


  SUBROUTINE Write1HDF_integer &
               ( name, values, group_id, datasize, desc_option, unit_option)

    CHARACTER(*), INTENT(in)                    :: name
    CHARACTER(*), INTENT(in), OPTIONAL          :: unit_option
    CHARACTER(*), INTENT(in), OPTIONAL          :: desc_option
    INTEGER(HID_T)                              :: group_id
    INTEGER(HSIZE_T), DIMENSION(1), INTENT(in)  :: datasize
    INTEGER, INTENT(in)                         :: values
   
!!$    INTEGER(HID_T)                              :: dataset_id
!!$    INTEGER(HID_T)                              :: dataspace_id
!!$    INTEGER(HID_T)                              :: atype_id
!!$    INTEGER(HID_T)                              :: attr_id
!!$    INTEGER(SIZE_T)                             :: attr_len
!!$    INTEGER(HSIZE_T), DIMENSION(1)              :: adims = (/1/)
!!$  
!!$    
!!$    CALL h5screate_simple_f( 1, datasize, dataspace_id, hdferr )
!!$
!!$    CALL h5dcreate_f( group_id, name, H5T_NATIVE_INTEGER, &
!!$           dataspace_id, dataset_id, hdferr )
!!$
!!$    CALL h5dwrite_f( dataset_id, H5T_NATIVE_INTEGER, &
!!$           values, datasize, hdferr )
!!$
!!$    CALL h5sclose_f( dataspace_id, hdferr ) 
!!$
!!$    CALL h5dclose_f( dataset_id, hdferr )

  END SUBROUTINE Write1HDF_integer


  SUBROUTINE Write1HDF_string &
               ( name, values, group_id, datasize, desc_option, unit_option )

    CHARACTER(*), INTENT(in)                    :: name
    CHARACTER(*), INTENT(in), OPTIONAL          :: unit_option
    CHARACTER(*), INTENT(in), OPTIONAL          :: desc_option
    INTEGER(HID_T)                              :: group_id
    INTEGER(HSIZE_T), DIMENSION(1), INTENT(in)  :: datasize
    CHARACTER(len=*), INTENT(in)                :: values
   
!!$    INTEGER(HSIZE_T)                            :: sizechar
!!$    INTEGER(HID_T)                              :: dataset_id
!!$    INTEGER(HID_T)                              :: dataspace_id
!!$    INTEGER(HID_T)                              :: atype_id
!!$    INTEGER(HID_T)                              :: attr_id
!!$    INTEGER(SIZE_T)                             :: attr_len
!!$    INTEGER(HSIZE_T), DIMENSION(1)              :: adims = (/1/)
!!$  
!!$    
!!$    CALL h5screate_simple_f( 1, datasize, dataspace_id, hdferr )
!!$    sizechar = LEN( values )
!!$
!!$    CALL h5tset_size_f( H5T_NATIVE_CHARACTER, sizechar, hdferr )
!!$    CALL h5dcreate_f( group_id, name, H5T_NATIVE_CHARACTER, &
!!$           dataspace_id, dataset_id, hdferr )
!!$
!!$    CALL h5dwrite_f( dataset_id, H5T_NATIVE_CHARACTER, &
!!$           values, datasize, hdferr )
!!$
!!$    CALL h5sclose_f( dataspace_id, hdferr ) 
!!$
!!$    CALL h5dclose_f( dataset_id, hdferr )

  END SUBROUTINE Write1HDF_string


  SUBROUTINE WriteEOSTableHDF( EOSTable, file_id, group_id )

    TYPE(EquationOfStateTableType), INTENT(in)    :: EOSTable
    INTEGER(HID_T), INTENT(in)                    :: file_id
    INTEGER(HID_T), INTENT(in)                    :: group_id

!!$    INTEGER(HID_T)                                :: subgroup_id
!!$
!!$    CALL OpenGroupHDF( "ThermoState", .true., group_id, subgroup_id )
!!$    CALL WriteThermoStateHDF( EOSTable % TS, subgroup_id )
!!$    CALL CloseGroupHDF( subgroup_id )
!!$
!!$    CALL OpenGroupHDF( "DependentVariables", .true., group_id, subgroup_id )
!!$    CALL WriteDependentVariablesHDF( EOSTable % DV, subgroup_id )
!!$    CALL CloseGroupHDF( subgroup_id )

  END SUBROUTINE WriteEOSTableHDF
  

  SUBROUTINE ReadOpacityTableHDF( OpacityTable, FileName )
 
    TYPE(OpacityTableType), INTENT(inout)   :: OpacityTable
    CHARACTER(len=*), INTENT(in)            :: FileName

!!$    INTEGER, DIMENSION(3)                         :: nPoints
!!$    INTEGER                                       :: nPointsE
!!$    INTEGER                                       :: nSpeciesA
!!$    INTEGER                                       :: nVariables
!!$    INTEGER(HID_T)                                :: file_id
!!$    INTEGER(HID_T)                                :: group_id
!!$    INTEGER(HID_T)                                :: subgroup_id
!!$
!!$
!!$    CALL OpenFileHDF( FileName, .false., file_id )
!!$
!!$    CALL OpenGroupHDF( "EOSTable", .false., file_id, group_id )
!!$   
!!$    Call OpenGroupHDF( "DependentVariables", .false., group_id, subgroup_id )     
!!$    CALL ReadDimensionsHDF( nPoints, subgroup_id )
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
!!$
!!$    CALL CloseFileHDF( file_id )

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
