MODULE wlHelmIOModuleHDF
    
    USE wlKindModule, ONLY: dp
    USE wlLeptonEOSTableModule, ONLY: &
        HelmTableType, &
        AllocateHelmholtzTable, &
        DeallocateHelmholtzTable
    USE wlIOModuleHDF
    USE HDF5
    
    IMPLICIT NONE
    PRIVATE
    
    PUBLIC :: WriteHelmholtzTableHDF
    PUBLIC :: ReadHelmholtzTableHDF
    
    INTEGER :: hdferr
    
CONTAINS
  
  SUBROUTINE WriteHelmholtzTableHDF( HelmTable, FileName, WhichTableName, NewFile )
    
    TYPE(HelmTableType), INTENT(INOUT)  :: HelmTable
    CHARACTER(len=*), INTENT(IN) :: FileName
    CHARACTER(len=*), INTENT(IN) :: WhichTableName
    LOGICAL, INTENT(IN) :: NewFile
    
    INTEGER(HID_T) :: file_id
    INTEGER(HID_T) :: group_id
    INTEGER(HSIZE_T) :: datasize2d(2)
    INTEGER(HSIZE_T) :: datasize1d(1)
    
    REAL(dp), DIMENSION(1) :: buffer
    
    CALL OpenFileHDF( FileName, NewFile, file_id, ReadWrite_Option = .TRUE. )
    
    datasize2d = (/ HelmTable % nPointsDen, HelmTable % nPointsTemp /) 
    
    CALL OpenGroupHDF( TRIM(ADJUSTL(WhichTableName)), .true., file_id, group_id )
    datasize1d = 2
    CALL WriteHDF( "DimensionsHelmTable", &
      (/ HelmTable % nPointsDen, HelmTable % nPointsTemp /), &
      group_id, datasize1d )
    
    datasize1d = 1
    buffer(1) = HelmTable % lepton_mass
    CALL WriteHDF( "LeptonMass", buffer, group_id, datasize1d )

    datasize1d = HelmTable % nPointsDen
    CALL WriteHDF( "Density", HelmTable % d, group_id, datasize1d )
    datasize1d = HelmTable % nPointsTemp
    CALL WriteHDF( "Temperature", HelmTable % t, group_id, datasize1d )
    CALL CloseGroupHDF( group_id )

    CALL OpenGroupHDF( TRIM(ADJUSTL(WhichTableName)) // "/FreeEnergyTable", .true., file_id, group_id )
    CALL WriteHDF( "f", HelmTable % f(:,:),     group_id, datasize2d )
    CALL WriteHDF( "fd", HelmTable % fd(:,:),    group_id, datasize2d )
    CALL WriteHDF( "ft", HelmTable % ft(:,:),    group_id, datasize2d )
    CALL WriteHDF( "fdd", HelmTable % fdd(:,:),   group_id, datasize2d )
    CALL WriteHDF( "ftt", HelmTable % ftt(:,:),   group_id, datasize2d )
    CALL WriteHDF( "fdt", HelmTable % fdt(:,:),   group_id, datasize2d )
    CALL WriteHDF( "fddt", HelmTable % fddt(:,:),  group_id, datasize2d )
    CALL WriteHDF( "fdtt", HelmTable % fdtt(:,:),  group_id, datasize2d )
    CALL WriteHDF( "fddtt", HelmTable % fddtt(:,:), group_id, datasize2d )
    CALL CloseGroupHDF( group_id )
    
    CALL OpenGroupHDF( TRIM(ADJUSTL(WhichTableName)) // "/dPdRhoTable", .true., file_id, group_id )
    CALL WriteHDF( "dpdf", HelmTable % dpdf(:,:),   group_id, datasize2d )
    CALL WriteHDF( "dpdfd", HelmTable % dpdfd(:,:),  group_id, datasize2d )
    CALL WriteHDF( "dpdft", HelmTable % dpdft(:,:),  group_id, datasize2d )
    CALL WriteHDF( "dpdfdt", HelmTable % dpdfdt(:,:), group_id, datasize2d )
    CALL CloseGroupHDF( group_id )

    CALL OpenGroupHDF( TRIM(ADJUSTL(WhichTableName)) // "/EleChemPotTable", .true., file_id, group_id )
    CALL WriteHDF( "ef", HelmTable % ef(:,:),   group_id, datasize2d )
    CALL WriteHDF( "efd", HelmTable % efd(:,:),  group_id, datasize2d )
    CALL WriteHDF( "eft", HelmTable % eft(:,:),  group_id, datasize2d )
    CALL WriteHDF( "efdt", HelmTable % efdt(:,:), group_id, datasize2d )
    CALL CloseGroupHDF( group_id )
    
    CALL OpenGroupHDF( TRIM(ADJUSTL(WhichTableName)) // "/NumberDensityTable", .true., file_id, group_id )
    CALL WriteHDF( "xf", HelmTable % xf(:,:),   group_id, datasize2d )
    CALL WriteHDF( "xfd", HelmTable % xfd(:,:),  group_id, datasize2d )
    CALL WriteHDF( "xft", HelmTable % xft(:,:),  group_id, datasize2d )
    CALL WriteHDF( "xfdt", HelmTable % xfdt(:,:), group_id, datasize2d )
    CALL CloseGroupHDF( group_id )
    
    datasize1d = HelmTable % nPointsTemp
    CALL OpenGroupHDF( TRIM(ADJUSTL(WhichTableName)) // "/DeltasTempTable", .true., file_id, group_id )
    CALL WriteHDF( "dt", HelmTable % dt(:),   group_id, datasize1d )
    CALL WriteHDF( "dt2", HelmTable % dt2(:),  group_id, datasize1d )
    CALL WriteHDF( "dti", HelmTable % dti(:),  group_id, datasize1d )
    CALL WriteHDF( "dt2i", HelmTable % dt2i(:), group_id, datasize1d )
    CALL CloseGroupHDF( group_id )
    
    datasize1d = HelmTable % nPointsDen
    CALL OpenGroupHDF( TRIM(ADJUSTL(WhichTableName)) // "/DeltasDenTable", .true., file_id, group_id )
    CALL WriteHDF( "dd", HelmTable % dd(:),   group_id, datasize1d )
    CALL WriteHDF( "dd2", HelmTable % dd2(:),  group_id, datasize1d )
    CALL WriteHDF( "ddi", HelmTable % ddi(:),  group_id, datasize1d )
    CALL WriteHDF( "dd2i", HelmTable % dd2i(:), group_id, datasize1d )
    CALL CloseGroupHDF( group_id )
    
    datasize1d = 1
    CALL OpenGroupHDF( TRIM(ADJUSTL(WhichTableName)) // "/LimitsTable", .true., file_id, group_id )
    buffer(1) = HelmTable % mintemp
    CALL WriteHDF( "mintemp", buffer, group_id, datasize1d )
    buffer(1) = HelmTable % maxtemp
    CALL WriteHDF( "maxtemp", buffer, group_id, datasize1d )
    buffer(1) = HelmTable % mindens
    CALL WriteHDF( "mindens", buffer, group_id, datasize1d )
    buffer(1) = HelmTable % maxdens
    CALL WriteHDF( "maxdens", buffer, group_id, datasize1d )
    CALL CloseGroupHDF( group_id )
    
    !CALL OpenGroupHDF( "Metadata", .true., file_id, group_id )
    !CALL WriteEOSMetadataHDF( HelmTable % MD, group_id )
    !CALL CloseGroupHDF( group_id )
    
    CALL CloseFileHDF( file_id )
        
  END SUBROUTINE WriteHelmholtzTableHDF
  
  SUBROUTINE ReadHelmholtzTableHDF( HelmTable, FileName, WhichTableName, eos_minD )
    
    TYPE(HelmTableType), INTENT(INOUT)  :: HelmTable
    CHARACTER(len=*), INTENT(IN) :: FileName
    CHARACTER(len=*), INTENT(IN) :: WhichTableName
    REAL(DP), INTENT(in), OPTIONAL :: eos_minD

    INTEGER(HID_T) :: file_id
    INTEGER(HID_T) :: group_id
    INTEGER(HSIZE_T) :: datasize2d(2)
    INTEGER(HSIZE_T) :: datasize1d(1)
    
    REAL(dp), DIMENSION(1) :: buffer
    INTEGER,  DIMENSION(2) :: nPoints
    
    CALL OpenFileHDF( FileName, .false., file_id )
    
    datasize1d(1) = 2
    CALL OpenGroupHDF( TRIM(ADJUSTL(WhichTableName)) // "", .false., file_id, group_id )
    CALL ReadHDF( "DimensionsHelmTable", nPoints(:), group_id, datasize1d )
    CALL CloseGroupHDF( group_id )
    
    ! Allocate Helmholtz EOS
    CALL AllocateHelmholtzTable( HelmTable, nPoints, eos_minD = eos_minD )

    CALL OpenGroupHDF( TRIM(ADJUSTL(WhichTableName)) // "", .false., file_id, group_id )
    datasize1d = HelmTable % nPointsDen
    CALL ReadHDF( "Density", HelmTable % d, group_id, datasize1d )
    datasize1d = HelmTable % nPointsTemp
    CALL ReadHDF( "Temperature", HelmTable % t, group_id, datasize1d )
    datasize1d = 1
    CALL ReadHDF( "LeptonMass", buffer, group_id, datasize1d )
    HelmTable % lepton_mass = buffer(1)
    CALL CloseGroupHDF( group_id )

    datasize2d = (/ HelmTable % nPointsDen, HelmTable % nPointsTemp /) 
    CALL OpenGroupHDF( TRIM(ADJUSTL(WhichTableName)) // "/FreeEnergyTable", .false., file_id, group_id )
    CALL ReadHDF( "f", HelmTable % f(:,:),     group_id, datasize2d )
    CALL ReadHDF( "fd", HelmTable % fd(:,:),    group_id, datasize2d )
    CALL ReadHDF( "ft", HelmTable % ft(:,:),    group_id, datasize2d )
    CALL ReadHDF( "fdd", HelmTable % fdd(:,:),   group_id, datasize2d )
    CALL ReadHDF( "ftt", HelmTable % ftt(:,:),   group_id, datasize2d )
    CALL ReadHDF( "fdt", HelmTable % fdt(:,:),   group_id, datasize2d )
    CALL ReadHDF( "fddt", HelmTable % fddt(:,:),  group_id, datasize2d )
    CALL ReadHDF( "fdtt", HelmTable % fdtt(:,:),  group_id, datasize2d )
    CALL ReadHDF( "fddtt", HelmTable % fddtt(:,:), group_id, datasize2d )
    CALL CloseGroupHDF( group_id )
    
    CALL OpenGroupHDF( TRIM(ADJUSTL(WhichTableName)) // "/dPdRhoTable", .false., file_id, group_id )
    CALL ReadHDF( "dpdf", HelmTable % dpdf(:,:),   group_id, datasize2d )
    CALL ReadHDF( "dpdfd", HelmTable % dpdfd(:,:),  group_id, datasize2d )
    CALL ReadHDF( "dpdft", HelmTable % dpdft(:,:),  group_id, datasize2d )
    CALL ReadHDF( "dpdfdt", HelmTable % dpdfdt(:,:), group_id, datasize2d )
    CALL CloseGroupHDF( group_id )
    
    CALL OpenGroupHDF( TRIM(ADJUSTL(WhichTableName)) // "/EleChemPotTable", .false., file_id, group_id )
    CALL ReadHDF( "ef", HelmTable % ef(:,:),   group_id, datasize2d )
    CALL ReadHDF( "efd", HelmTable % efd(:,:),  group_id, datasize2d )
    CALL ReadHDF( "eft", HelmTable % eft(:,:),  group_id, datasize2d )
    CALL ReadHDF( "efdt", HelmTable % efdt(:,:), group_id, datasize2d )
    CALL CloseGroupHDF( group_id )
    
    CALL OpenGroupHDF( TRIM(ADJUSTL(WhichTableName)) // "/NumberDensityTable", .false., file_id, group_id )
    CALL ReadHDF( "xf", HelmTable % xf(:,:),   group_id, datasize2d )
    CALL ReadHDF( "xfd", HelmTable % xfd(:,:),  group_id, datasize2d )
    CALL ReadHDF( "xft", HelmTable % xft(:,:),  group_id, datasize2d )
    CALL ReadHDF( "xfdt", HelmTable % xfdt(:,:), group_id, datasize2d )
    CALL CloseGroupHDF( group_id )
    
    datasize1d = HelmTable % nPointsTemp
    CALL OpenGroupHDF( TRIM(ADJUSTL(WhichTableName)) // "/DeltasTempTable", .false., file_id, group_id )
    CALL ReadHDF( "dt", HelmTable % dt(:),   group_id, datasize1d )
    CALL ReadHDF( "dt2", HelmTable % dt2(:),  group_id, datasize1d )
    CALL ReadHDF( "dti", HelmTable % dti(:),  group_id, datasize1d )
    CALL ReadHDF( "dt2i", HelmTable % dt2i(:), group_id, datasize1d )
    CALL CloseGroupHDF( group_id )
    
    datasize1d = HelmTable % nPointsDen
    CALL OpenGroupHDF( TRIM(ADJUSTL(WhichTableName)) // "/DeltasDenTable", .false., file_id, group_id )
    CALL ReadHDF( "dd", HelmTable % dd(:),   group_id, datasize1d )
    CALL ReadHDF( "dd2", HelmTable % dd2(:),  group_id, datasize1d )
    CALL ReadHDF( "ddi", HelmTable % ddi(:),  group_id, datasize1d )
    CALL ReadHDF( "dd2i", HelmTable % dd2i(:), group_id, datasize1d )
    CALL CloseGroupHDF( group_id )
    
    datasize1d = 1
    CALL OpenGroupHDF( TRIM(ADJUSTL(WhichTableName)) // "/LimitsTable", .false., file_id, group_id )
    CALL ReadHDF( "mintemp", buffer, group_id, datasize1d )
    HelmTable % mintemp = buffer(1)
    CALL ReadHDF( "maxtemp", buffer, group_id, datasize1d )
    HelmTable % maxtemp = buffer(1)
    CALL ReadHDF( "mindens", buffer, group_id, datasize1d )
    HelmTable % mindens = buffer(1)
    CALL ReadHDF( "maxdens", buffer, group_id, datasize1d )
    HelmTable % maxdens = buffer(1)
    CALL CloseGroupHDF( group_id )
        
    !CALL OpenGroupHDF( "Metadata", .true., file_id, group_id )
    !CALL WriteEOSMetadataHDF( HelmTable % MD, group_id )
    !CALL CloseGroupHDF( group_id )
    
    CALL CloseFileHDF( file_id )
    
  END SUBROUTINE ReadHelmholtzTableHDF

END MODULE wlHelmIOModuleHDF
        
