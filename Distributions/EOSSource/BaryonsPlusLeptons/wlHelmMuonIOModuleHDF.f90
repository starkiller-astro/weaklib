MODULE wlHelmMuonIOModuleHDF
    
    USE wlKindModule, ONLY: dp
    USE wlLeptonEOSModule, ONLY: &
        HelmTableType, &
        AllocateHelmholtzTable, DeallocateHelmholtzTable, &
        MuonTableType, &
        AllocateMuonTable, DeAllocateMuonTable
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
    
    SUBROUTINE WriteHelmholtzTableHDF( HelmTable, FileName, NewFile )
        
        TYPE(HelmTableType), INTENT(INOUT)  :: HelmTable
        CHARACTER(len=*), INTENT(IN) :: FileName
        LOGICAL, INTENT(IN) :: NewFile
        
        INTEGER(HID_T) :: file_id
        INTEGER(HID_T) :: group_id
        INTEGER(HSIZE_T) :: datasize2d(2)
        INTEGER(HSIZE_T) :: datasize1d(1)
        
        REAL(dp), DIMENSION(1) :: buffer
        
        CALL OpenFileHDF( FileName, NewFile, file_id, ReadWrite_Option = .TRUE. )
        
        datasize2d = (/ HelmTable % nPointsDen, HelmTable % nPointsTemp /) 
        
        CALL OpenGroupHDF( "HelmTable", .true., file_id, group_id )
        datasize1d = 2
        CALL WriteHDF( "DimensionsHelmTable", &
        (/ HelmTable % nPointsDen, HelmTable % nPointsTemp /), group_id, datasize1d )
        
        datasize1d = HelmTable % nPointsDen
        CALL WriteHDF( "Density", HelmTable % d, group_id, datasize1d )
        datasize1d = HelmTable % nPointsTemp
        CALL WriteHDF( "Temperature", HelmTable % t, group_id, datasize1d )
        CALL CloseGroupHDF( group_id )

        CALL OpenGroupHDF( "HelmTable/FreeEnergyTable", .true., file_id, group_id )
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
        
        CALL OpenGroupHDF( "HelmTable/dPdRhoTable", .true., file_id, group_id )
        CALL WriteHDF( "dpdf", HelmTable % dpdf(:,:),   group_id, datasize2d )
        CALL WriteHDF( "dpdfd", HelmTable % dpdfd(:,:),  group_id, datasize2d )
        CALL WriteHDF( "dpdft", HelmTable % dpdft(:,:),  group_id, datasize2d )
        CALL WriteHDF( "dpdfdt", HelmTable % dpdfdt(:,:), group_id, datasize2d )
        CALL CloseGroupHDF( group_id )
    
        CALL OpenGroupHDF( "HelmTable/EleChemPotTable", .true., file_id, group_id )
        CALL WriteHDF( "ef", HelmTable % ef(:,:),   group_id, datasize2d )
        CALL WriteHDF( "efd", HelmTable % efd(:,:),  group_id, datasize2d )
        CALL WriteHDF( "eft", HelmTable % eft(:,:),  group_id, datasize2d )
        CALL WriteHDF( "efdt", HelmTable % efdt(:,:), group_id, datasize2d )
        CALL CloseGroupHDF( group_id )
        
        CALL OpenGroupHDF( "HelmTable/NumberDensityTable", .true., file_id, group_id )
        CALL WriteHDF( "xf", HelmTable % xf(:,:),   group_id, datasize2d )
        CALL WriteHDF( "xfd", HelmTable % xfd(:,:),  group_id, datasize2d )
        CALL WriteHDF( "xft", HelmTable % xft(:,:),  group_id, datasize2d )
        CALL WriteHDF( "xfdt", HelmTable % xfdt(:,:), group_id, datasize2d )
        CALL CloseGroupHDF( group_id )
        
        datasize1d = HelmTable % nPointsTemp
        CALL OpenGroupHDF( "HelmTable/DeltasTempTable", .true., file_id, group_id )
        CALL WriteHDF( "dt", HelmTable % dt(:),   group_id, datasize1d )
        CALL WriteHDF( "dt2", HelmTable % dt2(:),  group_id, datasize1d )
        CALL WriteHDF( "dti", HelmTable % dti(:),  group_id, datasize1d )
        CALL WriteHDF( "dt2i", HelmTable % dt2i(:), group_id, datasize1d )
        CALL CloseGroupHDF( group_id )
        
        datasize1d = HelmTable % nPointsDen
        CALL OpenGroupHDF( "HelmTable/DeltasDenTable", .true., file_id, group_id )
        CALL WriteHDF( "dd", HelmTable % dd(:),   group_id, datasize1d )
        CALL WriteHDF( "dd2", HelmTable % dd2(:),  group_id, datasize1d )
        CALL WriteHDF( "ddi", HelmTable % ddi(:),  group_id, datasize1d )
        CALL WriteHDF( "dd2i", HelmTable % dd2i(:), group_id, datasize1d )
        CALL CloseGroupHDF( group_id )
        
        datasize1d = 1
        CALL OpenGroupHDF( "HelmTable/LimitsTable", .true., file_id, group_id )
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
    
    SUBROUTINE ReadHelmholtzTableHDF( HelmTable, FileName )
    
    TYPE(HelmTableType), INTENT(INOUT)  :: HelmTable
    CHARACTER(len=*), INTENT(IN) :: FileName
    
    INTEGER(HID_T) :: file_id
    INTEGER(HID_T) :: group_id
    INTEGER(HSIZE_T) :: datasize2d(2)
    INTEGER(HSIZE_T) :: datasize1d(1)
    
    REAL(dp), DIMENSION(1) :: buffer
    INTEGER, DIMENSION(2) :: nPoints
    
    CALL OpenFileHDF( FileName, .false., file_id )
    
    datasize1d(1) = 2
    CALL OpenGroupHDF( "HelmTable", .false., file_id, group_id )
    CALL ReadHDF( "DimensionsHelmTable", nPoints(:), group_id, datasize1d )
    CALL CloseGroupHDF( group_id )
    
    ! Allocate Helmholtz EOS
    CALL AllocateHelmholtzTable( HelmTable, nPoints )

    CALL OpenGroupHDF( "HelmTable", .false., file_id, group_id )
    datasize1d = HelmTable % nPointsDen
    CALL ReadHDF( "Density", HelmTable % d, group_id, datasize1d )
    datasize1d = HelmTable % nPointsTemp
    CALL ReadHDF( "Temperature", HelmTable % t, group_id, datasize1d )
    CALL CloseGroupHDF( group_id )

    datasize2d = (/ HelmTable % nPointsDen, HelmTable % nPointsTemp /) 
    CALL OpenGroupHDF( "HelmTable/FreeEnergyTable", .false., file_id, group_id )
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
    
    CALL OpenGroupHDF( "HelmTable/dPdRhoTable", .false., file_id, group_id )
    CALL ReadHDF( "dpdf", HelmTable % dpdf(:,:),   group_id, datasize2d )
    CALL ReadHDF( "dpdfd", HelmTable % dpdfd(:,:),  group_id, datasize2d )
    CALL ReadHDF( "dpdft", HelmTable % dpdft(:,:),  group_id, datasize2d )
    CALL ReadHDF( "dpdfdt", HelmTable % dpdfdt(:,:), group_id, datasize2d )
    CALL CloseGroupHDF( group_id )
    
    CALL OpenGroupHDF( "HelmTable/EleChemPotTable", .false., file_id, group_id )
    CALL ReadHDF( "ef", HelmTable % ef(:,:),   group_id, datasize2d )
    CALL ReadHDF( "efd", HelmTable % efd(:,:),  group_id, datasize2d )
    CALL ReadHDF( "eft", HelmTable % eft(:,:),  group_id, datasize2d )
    CALL ReadHDF( "efdt", HelmTable % efdt(:,:), group_id, datasize2d )
    CALL CloseGroupHDF( group_id )
    
    CALL OpenGroupHDF( "HelmTable/NumberDensityTable", .false., file_id, group_id )
    CALL ReadHDF( "xf", HelmTable % xf(:,:),   group_id, datasize2d )
    CALL ReadHDF( "xfd", HelmTable % xfd(:,:),  group_id, datasize2d )
    CALL ReadHDF( "xft", HelmTable % xft(:,:),  group_id, datasize2d )
    CALL ReadHDF( "xfdt", HelmTable % xfdt(:,:), group_id, datasize2d )
    CALL CloseGroupHDF( group_id )
    
    datasize1d = HelmTable % nPointsTemp
    CALL OpenGroupHDF( "HelmTable/DeltasTempTable", .false., file_id, group_id )
    CALL ReadHDF( "dt", HelmTable % dt(:),   group_id, datasize1d )
    CALL ReadHDF( "dt2", HelmTable % dt2(:),  group_id, datasize1d )
    CALL ReadHDF( "dti", HelmTable % dti(:),  group_id, datasize1d )
    CALL ReadHDF( "dt2i", HelmTable % dt2i(:), group_id, datasize1d )
    CALL CloseGroupHDF( group_id )
    
    datasize1d = HelmTable % nPointsDen
    CALL OpenGroupHDF( "HelmTable/DeltasDenTable", .false., file_id, group_id )
    CALL ReadHDF( "dd", HelmTable % dd(:),   group_id, datasize1d )
    CALL ReadHDF( "dd2", HelmTable % dd2(:),  group_id, datasize1d )
    CALL ReadHDF( "ddi", HelmTable % ddi(:),  group_id, datasize1d )
    CALL ReadHDF( "dd2i", HelmTable % dd2i(:), group_id, datasize1d )
    CALL CloseGroupHDF( group_id )
    
    datasize1d = 1
    CALL OpenGroupHDF( "HelmTable/LimitsTable", .false., file_id, group_id )
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
        
        datasize2d = (/ MuonTable % nPointsTemp, MuonTable % nPointsDen /) 
        CALL WriteHDF( "Mu", MuonTable % mu(:,:), group_id, datasize2d )
        CALL WriteHDF( "Density", MuonTable % rhoym(:,:), group_id, datasize2d )
        CALL WriteHDF( "Pressure", MuonTable % p(:,:),     group_id, datasize2d )
        CALL WriteHDF( "InternalEnergy", MuonTable % e(:,:),    group_id, datasize2d )
        CALL WriteHDF( "Entropy", MuonTable % s(:,:),    group_id, datasize2d )
        CALL CloseGroupHDF( group_id )
        CALL CloseFileHDF( file_id )
  
    END SUBROUTINE WriteMuonTableHDF
        
        
    SUBROUTINE ReadMuonTableHDF( MuonTable, FileName, eos_minD )
    
        TYPE(MuonTableType), INTENT(INOUT)  :: MuonTable
        CHARACTER(len=*), INTENT(IN) :: FileName
        REAL(DP), INTENT(in), OPTIONAL :: eos_minD
        
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
        CALL AllocateMuonTable( MuonTable, nPoints, eos_minD = eos_minD )

        CALL OpenGroupHDF( "MuonTable", .false., file_id, group_id )

        datasize1d = MuonTable % nPointsTemp
        CALL ReadHDF( "Temperature", MuonTable % t, group_id, datasize1d )

        datasize2d = (/ MuonTable % nPointsTemp, MuonTable % nPointsDen /) 
        CALL ReadHDF( "Mu", MuonTable % mu(:,:), group_id, datasize2d )
        CALL ReadHDF( "Density", MuonTable % rhoym(:,:), group_id, datasize2d )
        CALL ReadHDF( "Pressure", MuonTable % p(:,:),     group_id, datasize2d )
        CALL ReadHDF( "InternalEnergy", MuonTable % e(:,:),    group_id, datasize2d )
        CALL ReadHDF( "Entropy", MuonTable % s(:,:),    group_id, datasize2d )

        CALL CloseGroupHDF( group_id )
        CALL CloseFileHDF( file_id )
    
    END SUBROUTINE ReadMuonTableHDF

END MODULE wlHelmMuonIOModuleHDF
        
