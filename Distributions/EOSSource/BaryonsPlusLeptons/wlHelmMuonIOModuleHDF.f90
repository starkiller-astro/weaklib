MODULE wlHelmMuonIOModuleHDF
    
    USE wlKindModule, ONLY: dp
    USE wlLeptonEOSModule, ONLY: &
        HelmholtzTableType, &
        AllocateHelmholtzTable, DeallocateHelmholtzTable, &
        MuonTableType, &
        AllocateMuonEOS, DeAllocateMuonEOS
    USE wlIOModuleHDF
    USE HDF5
    
    IMPLICIT NONE
    PRIVATE
    
    PUBLIC :: WriteHelmholtzTableHDF
    PUBLIC :: ReadHelmholtzTableHDF
    
    PUBLIC :: WriteMuonTableHDF
    PUBLIC :: ReadMuonTableHDF
    
    INTEGER :: hdferr
    
CONTAINS
    
    SUBROUTINE WriteHelmholtzTableHDF( HelmholtzTable, FileName, NewFile )
        
        TYPE(HelmholtzTableType), INTENT(INOUT)  :: HelmholtzTable
        CHARACTER(len=*), INTENT(IN) :: FileName
        LOGICAL, INTENT(IN) :: NewFile
        
        INTEGER(HID_T) :: file_id
        INTEGER(HID_T) :: group_id
        INTEGER(HSIZE_T) :: datasize2d(2)
        INTEGER(HSIZE_T) :: datasize1d(1)
        
        REAL(dp), DIMENSION(1) :: buffer
        
        CALL OpenFileHDF( FileName, NewFile, file_id, ReadWrite_Option = .TRUE. )
        
        datasize2d = (/ HelmholtzTable % nPointsDen, HelmholtzTable % nPointsTemp /) 
        
        CALL OpenGroupHDF( "HelmholtzTable", .true., file_id, group_id )
        datasize1d = 2
        CALL WriteHDF( "DimensionsHelmTable", &
        (/ HelmholtzTable % nPointsDen, HelmholtzTable % nPointsTemp /), group_id, datasize1d )
        
        datasize1d = HelmholtzTable % nPointsDen
        CALL WriteHDF( "Density", HelmholtzTable % d, group_id, datasize1d )
        datasize1d = HelmholtzTable % nPointsTemp
        CALL WriteHDF( "Temperature", HelmholtzTable % t, group_id, datasize1d )
        CALL CloseGroupHDF( group_id )

        CALL OpenGroupHDF( "HelmholtzTable/FreeEnergyTable", .true., file_id, group_id )
        CALL WriteHDF( "f", HelmholtzTable % f(:,:),     group_id, datasize2d )
        CALL WriteHDF( "fd", HelmholtzTable % fd(:,:),    group_id, datasize2d )
        CALL WriteHDF( "ft", HelmholtzTable % ft(:,:),    group_id, datasize2d )
        CALL WriteHDF( "fdd", HelmholtzTable % fdd(:,:),   group_id, datasize2d )
        CALL WriteHDF( "ftt", HelmholtzTable % ftt(:,:),   group_id, datasize2d )
        CALL WriteHDF( "fdt", HelmholtzTable % fdt(:,:),   group_id, datasize2d )
        CALL WriteHDF( "fddt", HelmholtzTable % fddt(:,:),  group_id, datasize2d )
        CALL WriteHDF( "fdtt", HelmholtzTable % fdtt(:,:),  group_id, datasize2d )
        CALL WriteHDF( "fddtt", HelmholtzTable % fddtt(:,:), group_id, datasize2d )
        CALL CloseGroupHDF( group_id )
        
        CALL OpenGroupHDF( "HelmholtzTable/dPdRhoTable", .true., file_id, group_id )
        CALL WriteHDF( "dpdf", HelmholtzTable % dpdf(:,:),   group_id, datasize2d )
        CALL WriteHDF( "dpdfd", HelmholtzTable % dpdfd(:,:),  group_id, datasize2d )
        CALL WriteHDF( "dpdft", HelmholtzTable % dpdft(:,:),  group_id, datasize2d )
        CALL WriteHDF( "dpdfdt", HelmholtzTable % dpdfdt(:,:), group_id, datasize2d )
        CALL CloseGroupHDF( group_id )
    
        CALL OpenGroupHDF( "HelmholtzTable/EleChemPotTable", .true., file_id, group_id )
        CALL WriteHDF( "ef", HelmholtzTable % ef(:,:),   group_id, datasize2d )
        CALL WriteHDF( "efd", HelmholtzTable % efd(:,:),  group_id, datasize2d )
        CALL WriteHDF( "eft", HelmholtzTable % eft(:,:),  group_id, datasize2d )
        CALL WriteHDF( "efdt", HelmholtzTable % efdt(:,:), group_id, datasize2d )
        CALL CloseGroupHDF( group_id )
        
        CALL OpenGroupHDF( "HelmholtzTable/NumberDensityTable", .true., file_id, group_id )
        CALL WriteHDF( "xf", HelmholtzTable % xf(:,:),   group_id, datasize2d )
        CALL WriteHDF( "xfd", HelmholtzTable % xfd(:,:),  group_id, datasize2d )
        CALL WriteHDF( "xft", HelmholtzTable % xft(:,:),  group_id, datasize2d )
        CALL WriteHDF( "xfdt", HelmholtzTable % xfdt(:,:), group_id, datasize2d )
        CALL CloseGroupHDF( group_id )
        
        datasize1d = HelmholtzTable % nPointsTemp
        CALL OpenGroupHDF( "HelmholtzTable/DeltasTempTable", .true., file_id, group_id )
        CALL WriteHDF( "dt", HelmholtzTable % dt(:),   group_id, datasize1d )
        CALL WriteHDF( "dt2", HelmholtzTable % dt2(:),  group_id, datasize1d )
        CALL WriteHDF( "dti", HelmholtzTable % dti(:),  group_id, datasize1d )
        CALL WriteHDF( "dt2i", HelmholtzTable % dt2i(:), group_id, datasize1d )
        CALL CloseGroupHDF( group_id )
        
        datasize1d = HelmholtzTable % nPointsDen
        CALL OpenGroupHDF( "HelmholtzTable/DeltasDenTable", .true., file_id, group_id )
        CALL WriteHDF( "dd", HelmholtzTable % dd(:),   group_id, datasize1d )
        CALL WriteHDF( "dd2", HelmholtzTable % dd2(:),  group_id, datasize1d )
        CALL WriteHDF( "ddi", HelmholtzTable % ddi(:),  group_id, datasize1d )
        CALL WriteHDF( "dd2i", HelmholtzTable % dd2i(:), group_id, datasize1d )
        CALL CloseGroupHDF( group_id )
        
        datasize1d = 1
        CALL OpenGroupHDF( "HelmholtzTable/LimitsTable", .true., file_id, group_id )
        buffer(1) = HelmholtzTable % mintemp
        CALL WriteHDF( "mintemp", buffer, group_id, datasize1d )
        buffer(1) = HelmholtzTable % maxtemp
        CALL WriteHDF( "maxtemp", buffer, group_id, datasize1d )
        buffer(1) = HelmholtzTable % mindens
        CALL WriteHDF( "mindens", buffer, group_id, datasize1d )
        buffer(1) = HelmholtzTable % maxdens
        CALL WriteHDF( "maxdens", buffer, group_id, datasize1d )
        CALL CloseGroupHDF( group_id )
        
        !CALL OpenGroupHDF( "Metadata", .true., file_id, group_id )
        !CALL WriteEOSMetadataHDF( HelmholtzTable % MD, group_id )
        !CALL CloseGroupHDF( group_id )
        
        CALL CloseFileHDF( file_id )
            
    END SUBROUTINE WriteHelmholtzTableHDF
    
    SUBROUTINE ReadHelmholtzTableHDF( HelmholtzTable, FileName )
    
    TYPE(HelmholtzTableType), INTENT(INOUT)  :: HelmholtzTable
    CHARACTER(len=*), INTENT(IN) :: FileName
    
    INTEGER(HID_T) :: file_id
    INTEGER(HID_T) :: group_id
    INTEGER(HSIZE_T) :: datasize2d(2)
    INTEGER(HSIZE_T) :: datasize1d(1)
    
    REAL(dp), DIMENSION(1) :: buffer
    INTEGER, DIMENSION(2) :: nPoints
    
    CALL OpenFileHDF( FileName, .false., file_id )
    
    datasize1d(1) = 2
    CALL OpenGroupHDF( "HelmholtzTable", .false., file_id, group_id )
    CALL ReadHDF( "DimensionsHelmTable", nPoints(:), group_id, datasize1d )
    CALL CloseGroupHDF( group_id )
    
    ! Allocate Helmholtz EOS
    CALL AllocateHelmholtzTable( HelmholtzTable, nPoints )

    CALL OpenGroupHDF( "HelmholtzTable", .false., file_id, group_id )
    datasize1d = HelmholtzTable % nPointsDen
    CALL ReadHDF( "Density", HelmholtzTable % d, group_id, datasize1d )
    datasize1d = HelmholtzTable % nPointsTemp
    CALL ReadHDF( "Temperature", HelmholtzTable % t, group_id, datasize1d )
    CALL CloseGroupHDF( group_id )

    datasize2d = (/ HelmholtzTable % nPointsDen, HelmholtzTable % nPointsTemp /) 
    CALL OpenGroupHDF( "HelmholtzTable/FreeEnergyTable", .false., file_id, group_id )
    CALL ReadHDF( "f", HelmholtzTable % f(:,:),     group_id, datasize2d )
    CALL ReadHDF( "fd", HelmholtzTable % fd(:,:),    group_id, datasize2d )
    CALL ReadHDF( "ft", HelmholtzTable % ft(:,:),    group_id, datasize2d )
    CALL ReadHDF( "fdd", HelmholtzTable % fdd(:,:),   group_id, datasize2d )
    CALL ReadHDF( "ftt", HelmholtzTable % ftt(:,:),   group_id, datasize2d )
    CALL ReadHDF( "fdt", HelmholtzTable % fdt(:,:),   group_id, datasize2d )
    CALL ReadHDF( "fddt", HelmholtzTable % fddt(:,:),  group_id, datasize2d )
    CALL ReadHDF( "fdtt", HelmholtzTable % fdtt(:,:),  group_id, datasize2d )
    CALL ReadHDF( "fddtt", HelmholtzTable % fddtt(:,:), group_id, datasize2d )
    CALL CloseGroupHDF( group_id )
    
    CALL OpenGroupHDF( "HelmholtzTable/dPdRhoTable", .false., file_id, group_id )
    CALL ReadHDF( "dpdf", HelmholtzTable % dpdf(:,:),   group_id, datasize2d )
    CALL ReadHDF( "dpdfd", HelmholtzTable % dpdfd(:,:),  group_id, datasize2d )
    CALL ReadHDF( "dpdft", HelmholtzTable % dpdft(:,:),  group_id, datasize2d )
    CALL ReadHDF( "dpdfdt", HelmholtzTable % dpdfdt(:,:), group_id, datasize2d )
    CALL CloseGroupHDF( group_id )
    
    CALL OpenGroupHDF( "HelmholtzTable/EleChemPotTable", .false., file_id, group_id )
    CALL ReadHDF( "ef", HelmholtzTable % ef(:,:),   group_id, datasize2d )
    CALL ReadHDF( "efd", HelmholtzTable % efd(:,:),  group_id, datasize2d )
    CALL ReadHDF( "eft", HelmholtzTable % eft(:,:),  group_id, datasize2d )
    CALL ReadHDF( "efdt", HelmholtzTable % efdt(:,:), group_id, datasize2d )
    CALL CloseGroupHDF( group_id )
    
    CALL OpenGroupHDF( "HelmholtzTable/NumberDensityTable", .false., file_id, group_id )
    CALL ReadHDF( "xf", HelmholtzTable % xf(:,:),   group_id, datasize2d )
    CALL ReadHDF( "xfd", HelmholtzTable % xfd(:,:),  group_id, datasize2d )
    CALL ReadHDF( "xft", HelmholtzTable % xft(:,:),  group_id, datasize2d )
    CALL ReadHDF( "xfdt", HelmholtzTable % xfdt(:,:), group_id, datasize2d )
    CALL CloseGroupHDF( group_id )
    
    datasize1d = HelmholtzTable % nPointsTemp
    CALL OpenGroupHDF( "HelmholtzTable/DeltasTempTable", .false., file_id, group_id )
    CALL ReadHDF( "dt", HelmholtzTable % dt(:),   group_id, datasize1d )
    CALL ReadHDF( "dt2", HelmholtzTable % dt2(:),  group_id, datasize1d )
    CALL ReadHDF( "dti", HelmholtzTable % dti(:),  group_id, datasize1d )
    CALL ReadHDF( "dt2i", HelmholtzTable % dt2i(:), group_id, datasize1d )
    CALL CloseGroupHDF( group_id )
    
    datasize1d = HelmholtzTable % nPointsDen
    CALL OpenGroupHDF( "HelmholtzTable/DeltasDenTable", .false., file_id, group_id )
    CALL ReadHDF( "dd", HelmholtzTable % dd(:),   group_id, datasize1d )
    CALL ReadHDF( "dd2", HelmholtzTable % dd2(:),  group_id, datasize1d )
    CALL ReadHDF( "ddi", HelmholtzTable % ddi(:),  group_id, datasize1d )
    CALL ReadHDF( "dd2i", HelmholtzTable % dd2i(:), group_id, datasize1d )
    CALL CloseGroupHDF( group_id )
    
    datasize1d = 1
    CALL OpenGroupHDF( "HelmholtzTable/LimitsTable", .false., file_id, group_id )
    CALL ReadHDF( "mintemp", buffer, group_id, datasize1d )
    HelmholtzTable % mintemp = buffer(1)
    CALL ReadHDF( "maxtemp", buffer, group_id, datasize1d )
    HelmholtzTable % maxtemp = buffer(1)
    CALL ReadHDF( "mindens", buffer, group_id, datasize1d )
    HelmholtzTable % mindens = buffer(1)
    CALL ReadHDF( "maxdens", buffer, group_id, datasize1d )
    HelmholtzTable % maxdens = buffer(1)
    CALL CloseGroupHDF( group_id )
        
    !CALL OpenGroupHDF( "Metadata", .true., file_id, group_id )
    !CALL WriteEOSMetadataHDF( HelmholtzTable % MD, group_id )
    !CALL CloseGroupHDF( group_id )
    
    CALL CloseFileHDF( file_id )
    
    END SUBROUTINE ReadHelmholtzTableHDF
        
    SUBROUTINE WriteMuonTableHDF( MuonTable, FileName, NewFile )
        
        TYPE(MuonTableType), INTENT(INOUT)  :: MuonTable
        CHARACTER(len=*), INTENT(IN) :: FileName
        LOGICAL, INTENT(IN) :: NewFile
        
        INTEGER(HID_T) :: file_id
        INTEGER(HID_T) :: group_id
        INTEGER(HSIZE_T) :: datasize2d(2)
        INTEGER(HSIZE_T) :: datasize1d(1)
                
        CALL OpenFileHDF( FileName, NewFile, file_id, ReadWrite_Option = .TRUE. )
        
        datasize2d = (/ MuonTable % nPointsTemp, MuonTable % nPointsDen /) 
        
        CALL OpenGroupHDF( "MuonTable", .true., file_id, group_id )
        datasize1d = 2
        CALL WriteHDF( "DimensionsMuonTable", &
        (/ MuonTable % nPointsTemp, MuonTable % nPointsDen /), group_id, datasize1d )

        datasize1d = MuonTable % nPointsTemp
        CALL WriteHDF( "Temperature", MuonTable % t, group_id, datasize1d )
        
        datasize1d = MuonTable % nPointsDen
        CALL WriteHDF( "Density", MuonTable % rhoym, group_id, datasize1d )

        datasize2d = (/ MuonTable % nPointsTemp, MuonTable % nPointsDen /) 
        CALL WriteHDF( "Pressure", MuonTable % p(:,:),     group_id, datasize2d )
        CALL WriteHDF( "InternalEnergy", MuonTable % e(:,:),    group_id, datasize2d )
        CALL WriteHDF( "Entropy", MuonTable % s(:,:),    group_id, datasize2d )
        CALL WriteHDF( "Mu", MuonTable % mu, group_id, datasize2d )
        CALL CloseGroupHDF( group_id )
        CALL CloseFileHDF( file_id )
  
    END SUBROUTINE WriteMuonTableHDF
        
        
    SUBROUTINE ReadMuonTableHDF( MuonTable, FileName )
    
        TYPE(MuonTableType), INTENT(INOUT)  :: MuonTable
        CHARACTER(len=*), INTENT(IN) :: FileName
        
        INTEGER(HID_T) :: file_id
        INTEGER(HID_T) :: group_id
        INTEGER(HSIZE_T) :: datasize2d(2)
        INTEGER(HSIZE_T) :: datasize1d(1)
        
        REAL(dp), DIMENSION(1) :: buffer
        INTEGER, DIMENSION(2) :: nPoints
        
        CALL OpenFileHDF( FileName, .false., file_id )
        
        datasize1d(1) = 2
        CALL OpenGroupHDF( "MuonTable", .false., file_id, group_id )
        CALL ReadHDF( "DimensionsMuonTable", nPoints(:), group_id, datasize1d )
        CALL CloseGroupHDF( group_id )
        
        ! Allocate Muon EOS
        CALL AllocateMuonEOS( MuonTable, nPoints )

        CALL OpenGroupHDF( "MuonTable", .false., file_id, group_id )

        datasize1d = MuonTable % nPointsTemp
        CALL ReadHDF( "Temperature", MuonTable % t, group_id, datasize1d )
        CALL ReadHDF( "Density", MuonTable % rhoym, group_id, datasize1d )

        datasize2d = (/ MuonTable % nPointsTemp, MuonTable % nPointsDen /) 
        CALL ReadHDF( "Pressure", MuonTable % p(:,:),     group_id, datasize2d )
        CALL ReadHDF( "InternalEnergy", MuonTable % e(:,:),    group_id, datasize2d )
        CALL ReadHDF( "Entropy", MuonTable % s(:,:),    group_id, datasize2d )
        CALL ReadHDF( "Mu", MuonTable % mu, group_id, datasize2d )

        CALL CloseGroupHDF( group_id )
        CALL CloseFileHDF( file_id )
    
    END SUBROUTINE ReadMuonTableHDF

END MODULE wlHelmMuonIOModuleHDF
        