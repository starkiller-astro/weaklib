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
!-----------------------------------------------------------------------
!  NOTE: Only Type A interaction applied. Type B and Type C interaction 
!        needs to be added for future use.
!-----------------------------------------------------------------------

  USE wlKindModule, ONLY:         &
    dp
  USE wlGridModule, ONLY:         &
    GridType
  USE wlOpacityTableModule, ONLY: &
    OpacityTableType,             &
    AllocateOpacityTable
  USE wlOpacityFieldsModule, ONLY:&
    OpacityTypeEmAb,              &
    OpacityTypeScat,                 &
    OpacityTypeC
  USE wlIOModuleHDF, ONLY:        &
    ReadHDF,                      &
    WriteHDF,                     &
    OpenFileHDF,                  &
    CloseFileHDF,                 &
    OpenGroupHDF,                 &
    CloseGroupHDF,                &
    WriteThermoStateHDF,          &
    ReadThermoStateHDF
  USE wlEquationOfStateTableModule
  USE HDF5

  IMPLICIT NONE
  PRIVATE

  INTEGER :: hdferr

  PUBLIC WriteOpacityTableHDF
  PUBLIC WriteOpacityTableHDF_New
  PUBLIC ReadOpacityTableHDF
  PUBLIC ReadOpacityTableHDF_New

CONTAINS

  SUBROUTINE WriteOpacityTableHDF_New &
    ( OpacityTable, FileName, WriteOpacity_EmAb_Option )
 
    TYPE(OpacityTableType), INTENT(inout)        :: OpacityTable
    CHARACTER(len=*),       INTENT(in)           :: FileName
    LOGICAL,                INTENT(in), OPTIONAL :: WriteOpacity_EmAb_Option

    LOGICAL           :: WriteOpacity_EmAb
    CHARACTER(LEN=32) :: tempString(1)
    INTEGER           :: tempInteger(1)
    INTEGER(HID_T)    :: file_id
    INTEGER(HID_T)    :: group_id
    INTEGER(HSIZE_T)  :: datasize1d(1)

    IF( PRESENT( WriteOpacity_EmAb_Option ) )THEN
      WriteOpacity_EmAb = WriteOpacity_EmAb_Option
    ELSE
      WriteOpacity_EmAb = .FALSE.
    END IF
   
    CALL OpenFileHDF( FileName, .true., file_id )

    datasize1d(1) = 1

    CALL OpenGroupHDF( "EnergyGrid", .true., file_id, group_id )
    CALL WriteGridHDF( OpacityTable % EnergyGrid, group_id )
    CALL CloseGroupHDF( group_id )
  
    CALL OpenGroupHDF( "ThermoState", .true., file_id, group_id )
    CALL WriteThermoStateHDF( OpacityTable % TS, group_id )
    CALL CloseGroupHDF( group_id )

    IF( WriteOpacity_EmAb )THEN

      IF( .NOT. ALLOCATED( OpacityTable % EmAb % Names ) )THEN

        ! --- Insert Appropriate Reaction ---
        WRITE(*,'(A4,A)') &
          '', 'OpacityTable % EmAb not allocated.  Write Skipped.'

      ELSE

        CALL OpenGroupHDF &
               ( "EmAb_CorrectedAbsorption", .true., file_id, group_id )
        CALL WriteOpacityTableHDF_EmAb( OpacityTable % EmAb, group_id )
        CALL CloseGroupHDF( group_id )

      END IF

    END IF

    CALL CloseFileHDF( file_id )

  END SUBROUTINE WriteOpacityTableHDF_New


  SUBROUTINE WriteOpacityTableHDF( OpacityTable, FileName )
 
    TYPE(OpacityTableType), INTENT(inout)       :: OpacityTable
    CHARACTER(len=*), INTENT(in)                :: FileName

    INTEGER(HID_T)                              :: file_id
    INTEGER(HID_T)                              :: group_id

    CHARACTER(LEN=32), DIMENSION(1)             :: tempString
    INTEGER, DIMENSION(1)                       :: tempInteger
    INTEGER(HSIZE_T), DIMENSION(1)              :: datasize1d
   
    CALL OpenFileHDF( FileName, .true., file_id )

    datasize1d(1) = 1

    tempInteger(1) = OpacityTable % nOpacitiesA
    CALL WriteHDF&
         ( "nOpacitiesA", tempInteger, file_id, datasize1d )

    tempInteger(1) = OpacityTable % nOpacitiesB 
    CALL WriteHDF&
         ( "nOpacitiesB", tempInteger, file_id, datasize1d )
  
    tempInteger(1) = OpacityTable % nMomentsB     
    CALL WriteHDF&
         ( "nMomentsB", tempInteger, file_id, datasize1d )

    tempInteger(1) = OpacityTable % nOpacitiesB_NES
    CALL WriteHDF&
         ( "nOpacitiesB_NES", tempInteger, file_id, datasize1d )

    tempInteger(1) = OpacityTable % nMomentsB_NES
    CALL WriteHDF&
         ( "nMomentsB_NES", tempInteger, file_id, datasize1d )

    tempInteger(1) = OpacityTable % nOpacitiesB_TP
    CALL WriteHDF&
         ( "nOpacitiesB_TP", tempInteger, file_id, datasize1d )

    tempInteger(1) = OpacityTable % nMomentsB_TP
    CALL WriteHDF&
         ( "nMomentsB_TP", tempInteger, file_id, datasize1d )
 
    tempInteger(1) = OpacityTable % nOpacitiesC     
    CALL WriteHDF&
         ( "nOpacitiesC", tempInteger, file_id, datasize1d )

    tempInteger(1) = OpacityTable % nMomentsC   
    CALL WriteHDF&
         ( "nMomentsC", tempInteger, file_id, datasize1d )

    tempInteger(1) = OpacityTable % nPointsE  
    CALL WriteHDF&
         ( "nPointsE", tempInteger, file_id, datasize1d )

    tempInteger(1) = OpacityTable % nPointsEta
    CALL WriteHDF&
         ( "nPointsEta", tempInteger, file_id, datasize1d )

    datasize1d = 3
    CALL WriteHDF&
         ( "nPointsTS", OpacityTable % nPointsTS, file_id, datasize1d )

    CALL OpenGroupHDF( "EnergyGrid", .true., file_id, group_id )
    CALL WriteGridHDF( OpacityTable % EnergyGrid, group_id )
    CALL CloseGroupHDF( group_id )

    CALL OpenGroupHDF( "EtaGrid", .true., file_id, group_id )
    CALL WriteGridHDF( OpacityTable % EtaGrid, group_id )
    CALL CloseGroupHDF( group_id )
  
    CALL OpenGroupHDF( "ThermoState", .true., file_id, group_id )
    CALL WriteThermoStateHDF( OpacityTable % TS, group_id )
    CALL CloseGroupHDF( group_id )

    CALL OpenGroupHDF( "EmAb", .true., file_id, group_id )
    CALL WriteOpacityTableHDF_EmAb( OpacityTable % EmAb, group_id )
    CALL CloseGroupHDF( group_id )

    CALL OpenGroupHDF( "Scat_Iso", .true., file_id, group_id )
    CALL WriteOpacityTableHDF_Scat( OpacityTable % Scat_Iso, group_id )
    CALL CloseGroupHDF( group_id )

    CALL OpenGroupHDF( "Scat_NES", .true., file_id, group_id )
    CALL WriteOpacityTableHDF_Scat( OpacityTable % Scat_NES, group_id )
    CALL CloseGroupHDF( group_id )

    CALL OpenGroupHDF( "Scat_Pair", .true., file_id, group_id )
    CALL WriteOpacityTableHDF_Scat( OpacityTable % Scat_Pair, group_id )
    CALL CloseGroupHDF( group_id )

    CALL OpenGroupHDF( "Scat_nIso", .true., file_id, group_id )
    CALL WriteOpacityTableTypeCHDF( OpacityTable % Scat_nIso, group_id )
    CALL CloseGroupHDF( group_id )   

    CALL CloseFileHDF( file_id )

  END SUBROUTINE WriteOpacityTableHDF


  SUBROUTINE WriteGridHDF( Grid, group_id )

    TYPE(GridType), INTENT(in)                  :: Grid
    INTEGER(HID_T), INTENT(in)                  :: group_id

    CHARACTER(LEN=32)                           :: tempString(1)
    INTEGER                                     :: tempInteger(1)
    INTEGER(HSIZE_T)                            :: datasize1d(1)

    datasize1d(1) = 1

    tempString(1) = Grid % Name
    CALL WriteHDF( "Name",      tempString,       group_id, datasize1d )
    
    tempString(1) = Grid % Unit
    CALL WriteHDF( "Unit",      tempString,       group_id, datasize1d )

    tempInteger(1) = Grid % nPoints  
    CALL WriteHDF( "nPoints",   tempInteger,      group_id, datasize1d )
   
    tempInteger(1) = Grid % LogInterp 
    CALL WriteHDF( "LogInterp", tempInteger,      group_id, datasize1d )
   
    datasize1d(1) = Grid % nPoints
    CALL WriteHDF( "Values",    Grid % Values(:), group_id, datasize1d )

  END SUBROUTINE WriteGridHDF


  SUBROUTINE WriteOpacityTableHDF_EmAb( EmAb, group_id )

    TYPE(OpacityTypeEmAb), INTENT(in) :: EmAb
    INTEGER(HID_T),        INTENT(in) :: group_id

    INTEGER(HSIZE_T) :: datasize1d(1)
    INTEGER(HSIZE_T) :: datasize4d(4)

    datasize1d = EmAb % nOpacities
    CALL WriteHDF &
           ( "Units", EmAb % Units, group_id, datasize1d ) 

    CALL WriteHDF &
           ( "Offsets", EmAb % Offsets, group_id, datasize1d )

    datasize4d = EmAb % nPoints

    CALL WriteHDF &
           ( TRIM( EmAb % Names(1) ), &
             EmAb % Absorptivity(1) % Values(:,:,:,:), group_id, datasize4d )

    CALL WriteHDF &
           ( TRIM( EmAb % Names(2) ), &
             EmAb % Absorptivity(2) % Values(:,:,:,:), group_id, datasize4d )
  
  END SUBROUTINE WriteOpacityTableHDF_EmAb


  SUBROUTINE WriteOpacityTableHDF_Scat( Scat_Iso , group_id )

    TYPE(OpacityTypeScat), INTENT(in)              :: Scat_Iso
    INTEGER(HID_T), INTENT(in)                  :: group_id

    INTEGER(HSIZE_T)                            :: datasize1d
    INTEGER(HSIZE_T), DIMENSION(2)              :: datasize2d
    INTEGER(HSIZE_T), DIMENSION(4)              :: datasize4d
    INTEGER(HSIZE_T), DIMENSION(5)              :: datasize5d
    INTEGER                                     :: i
    INTEGER, DIMENSION(1)                       :: buffer

    CHARACTER(LEN=32), DIMENSION(1)             :: tempString
    INTEGER, DIMENSION(1)                       :: tempInteger
    REAL(dp), DIMENSION(1)                      :: tempReal
    INTEGER(HSIZE_T), DIMENSION(1)              :: datasize1dtemp
    INTEGER(HID_T)                              :: subgroup_id

    datasize1dtemp(1) = 1
    tempInteger(1) = Scat_Iso % nOpacities
    CALL WriteHDF&
         ( "nOpacities", tempInteger, group_id, datasize1dtemp )

    tempInteger(1) = Scat_Iso % nMoments
    CALL WriteHDF&
         ( "nMoments", tempInteger, group_id, datasize1dtemp )

    datasize1dtemp(1) = 4
    CALL WriteHDF&
         ( "nPoints", Scat_Iso % nPoints, group_id, datasize1dtemp )

    datasize1dtemp(1) = Scat_Iso % nOpacities
    CALL WriteHDF&
         ( "Names", Scat_Iso % Names, group_id, datasize1dtemp )

    CALL WriteHDF&
         ( "Species", Scat_Iso % Species, group_id, datasize1dtemp )

    CALL WriteHDF&
         ( "Units", Scat_Iso % Units, group_id, datasize1dtemp )

    datasize2d = (/Scat_Iso % nOpacities, Scat_Iso % nMoments/)
    CALL WriteHDF&
         ( "Offsets", Scat_Iso % Offsets, group_id, datasize2d )

    datasize1d = Scat_Iso % nOpacities
    datasize5d(1:4) = Scat_Iso % nPoints
    datasize5d(5) = Scat_Iso % nMoments

    CALL OpenGroupHDF( "Kernel", .true., group_id, subgroup_id )
    DO i = 1, datasize1d
     CALL WriteHDF&
        ( Scat_Iso % Names(i), Scat_Iso % Kernel(i) % Values(:,:,:,:,:),&
                            subgroup_id, datasize5d )
    END DO
    CALL CloseGroupHDF( subgroup_id )

  END SUBROUTINE WriteOpacityTableHDF_Scat


  SUBROUTINE WriteOpacityTableTypeCHDF( scattn, group_id )

    TYPE(OpacityTypeC), INTENT(in)              :: scattn
    INTEGER(HID_T), INTENT(in)                  :: group_id

  END SUBROUTINE WriteOpacityTableTypeCHDF


  SUBROUTINE ReadOpacityTableHDF_New &
    ( OpacityTable, FileName, ReadOpacity_EmAb_Option )
 
    TYPE(OpacityTableType), INTENT(inout)        :: OpacityTable
    CHARACTER(len=*),       INTENT(in)           :: FileName
    LOGICAL,                INTENT(in), OPTIONAL :: ReadOpacity_EmAb_Option

    INTEGER, DIMENSION(3)                         :: nPointsTS
    INTEGER                                       :: nPointsE
    INTEGER                                       :: nPointsEta
    INTEGER                                       :: nOpacA
    INTEGER                                       :: nOpacB, nMomB
    INTEGER                                       :: nOpacB_NES, nMomB_NES
    INTEGER                                       :: nOpacB_TP, nMomB_TP
    INTEGER(HID_T)                                :: file_id
    INTEGER(HID_T)                                :: group_id
    INTEGER(HID_T)                                :: subgroup_id
    INTEGER(HSIZE_T), DIMENSION(1)                :: datasize1d
    INTEGER(HSIZE_T), DIMENSION(4)                :: datasize4d
    INTEGER, DIMENSION(1)                         :: buffer
    CHARACTER(LEN=32), DIMENSION(1)               :: buffer_string

    INTEGER                                       :: hdfreadErr

    hdfreadErr = 0

    WRITE(*,*) "           File in"
    WRITE(*,*) " Reading ", FileName, " hdf5 file ... "

    CALL OpenFileHDF( FileName, .false., file_id )

    CALL OpenGroupHDF( "EnergyGrid", .false., file_id, group_id )
    CALL ReadHDF( "nPoints", buffer, group_id, datasize1d )
    nPointsE = buffer(1)
    CALL CloseGroupHDF( group_id )

    CALL AllocateOpacityTable &
           ( OpacityTable, 2, 0, 0, 0, 0, &
             0, 0, 0, 0, nPointsE, 0 )
   
!    IF( hdfreadErr == 0 ) THEN

      WRITE(*,*) "Now read-in EmAb table"
      CALL OpenGroupHDF( "EnergyGrid", .false., file_id, group_id )
      CALL ReadGridHDF( OpacityTable % EnergyGrid, group_id )
      CALL CloseGroupHDF( group_id )

      ! --- IF( ReadOpacity_EmAb )THEN

      CALL OpenGroupHDF &
             ( "EmAb_CorrectedAbsorption", .false., file_id, group_id )

      datasize1d(1) = 2
      CALL ReadHDF &
             ( "Offsets", OpacityTable % EmAb % Offsets, group_id, datasize1d )

      CALL ReadHDF &
             ( "Units",   OpacityTable % EmAb % Units,   group_id, datasize1d )

      datasize4d(1) = nPointsE
      datasize4d(2:4) = OpacityTable % TS % nPoints

      OpacityTable % EmAb % Names(1) = "Electron Neutrino";
      CALL ReadHDF &
             ( TRIM( OpacityTable % EmAb % Names(1) ), &
               OpacityTable % EmAb % Absorptivity(1) % Values, &
               group_id, datasize4d )

      OpacityTable % EmAb % Names(2) = "Electron Antineutrino";
      CALL ReadHDF &
         ( "Electron Antineutrino", &
            OpacityTable % EmAb % Absorptivity(2) % Values, &
            group_id, datasize4d )
      CALL CloseGroupHDF( group_id )

      CALL CloseFileHDF( file_id )
!    ELSE 
!      WRITE(*,*) "ERROR!"
!      WRITE(*,*) "EquationOfStateTable is not consistent with OpacityTable!"
!      CALL CloseFileHDF( file_id )
!      STOP
!    END IF

  END SUBROUTINE ReadOpacityTableHDF_New


  SUBROUTINE ReadOpacityTableHDF( OpacityTable, FileName, Verbose_Option )
 
    TYPE(OpacityTableType), INTENT(inout)        :: OpacityTable
    CHARACTER(len=*),       INTENT(in)           :: FileName
    LOGICAL,                INTENT(in), OPTIONAL :: Verbose_Option

    CHARACTER(LEN=32) :: buffer_string(1)
    LOGICAL           :: Verbose
    INTEGER           :: nPointsTS(3)
    INTEGER           :: nPointsE
    INTEGER           :: nPointsEta
    INTEGER           :: nOpacA
    INTEGER           :: nOpacB, nMomB
    INTEGER           :: nOpacB_NES, nMomB_NES
    INTEGER           :: nOpacB_TP, nMomB_TP
    INTEGER           :: nOpacC, nMomC
    INTEGER           :: buffer(1)
    INTEGER           :: hdfreadErr
    INTEGER(HID_T)    :: file_id
    INTEGER(HID_T)    :: group_id
    INTEGER(HID_T)    :: subgroup_id
    INTEGER(HSIZE_T)  :: datasize1d(1)

    IF( PRESENT( Verbose_Option ) )THEN
      Verbose = Verbose_Option
    ELSE
      Verbose = .FALSE.
    END IF

    hdfreadErr = 0

    IF( Verbose )THEN
      WRITE(*,*)
      WRITE(*,'(A4,A9,A)') '', 'Reading: ', TRIM( FileName )
    END IF

    CALL OpenFileHDF( FileName, .false., file_id )

    datasize1d(1) = 1
    CALL ReadHDF( "nOpacitiesA", buffer, file_id, datasize1d )
    nOpacA = buffer(1)   

    CALL ReadHDF( "nOpacitiesB", buffer, file_id, datasize1d )
    nOpacB = buffer(1)

    CALL ReadHDF( "nMomentsB", buffer, file_id, datasize1d )
    nMomB = buffer(1)

    CALL ReadHDF( "nOpacitiesB_NES", buffer, file_id, datasize1d )
    nOpacB_NES = buffer(1)

    CALL ReadHDF( "nMomentsB_NES", buffer, file_id, datasize1d )
    nMomB_NES = buffer(1)

    CALL ReadHDF &
           ( "nOpacitiesB_TP", buffer, file_id, datasize1d, hdfreadErr )
    IF( hdfreadErr == 0 ) THEN 
      nOpacB_TP = buffer(1)
    ELSE
      nOpacB_TP = 0
      PRINT*, "hdfreadErr =", hdfreadErr
      PRINT*, "nOpacB_TP =", nOpacB_TP
    END IF

    CALL ReadHDF &
           ( "nMomentsB_TP", buffer, file_id, datasize1d, hdfreadErr )
    IF( hdfreadErr == 0 ) THEN 
      nMomB_TP = buffer(1)
    ELSE
      nMomB_TP = 0
      PRINT*, "hdfreadErr =", hdfreadErr
      PRINT*, "nMomB_TP =", nMomB_TP
    END IF

    CALL ReadHDF( "nOpacitiesC", buffer, file_id, datasize1d )
    nOpacC = buffer(1)

    CALL ReadHDF( "nMomentsC", buffer, file_id, datasize1d )
    nMomC = buffer(1)

    CALL ReadHDF( "nPointsE", buffer, file_id, datasize1d )
    nPointsE = buffer(1)

    CALL ReadHDF( "nPointsEta", buffer, file_id, datasize1d )
    nPointsEta = buffer(1)

    CALL ReadHDF( "nPointsTS", nPointsTS, file_id, datasize1d )

    CALL AllocateOpacityTable &
           ( OpacityTable, nOpacA, nOpacB, nMomB, nOpacB_NES, nMomB_NES, &
             nOpacB_TP, nMomB_TP, nOpacC, nMomC, nPointsE, nPointsEta, &
             Verbose_Option = Verbose_Option )  

    ASSOCIATE &
      ( nPointsTS_EOS => OpacityTable % EOSTable % TS % nPoints )

    IF( ALL( nPointsTS_EOS .EQ. nPointsTS ) )THEN

      CALL OpenGroupHDF( "EnergyGrid", .false., file_id, group_id )
      CALL ReadGridHDF( OpacityTable % EnergyGrid, group_id )
      CALL CloseGroupHDF( group_id )

      CALL OpenGroupHDF( "EtaGrid", .false., file_id, group_id )
      CALL ReadGridHDF( OpacityTable % EtaGrid, group_id )
      CALL CloseGroupHDF( group_id )
 
      CALL ReadThermoStateHDF( OpacityTable % TS, file_id )

      IF( nOpacA .ne. 0 ) THEN
        CALL OpenGroupHDF( "EmAb", .false., file_id, group_id )
        CALL ReadOpacityTypeEmAbHDF( OpacityTable % EmAb, group_id )
        CALL CloseGroupHDF( group_id )
      END IF

      IF( nOpacB .ne. 0 ) THEN
        CALL OpenGroupHDF( "Scat_Iso", .false., file_id, group_id )
        CALL ReadOpacityTypeScatHDF( OpacityTable % Scat_Iso, group_id )
        CALL CloseGroupHDF( group_id )
      END IF

      IF( nOpacB_NES .ne. 0 ) THEN
        CALL OpenGroupHDF( "Scat_NES", .false., file_id, group_id )
        CALL ReadOpacityTypeScatHDF( OpacityTable % Scat_NES, group_id )
        CALL CloseGroupHDF( group_id )
      END IF

      IF( nOpacB_TP .ne. 0 ) THEN
        CALL OpenGroupHDF( "Scat_Pair", .false., file_id, group_id )
        CALL ReadOpacityTypeScatHDF( OpacityTable % Scat_Pair, group_id )
        CALL CloseGroupHDF( group_id )
      END IF

      CALL OpenGroupHDF( "Scat_nIso", .false., file_id, group_id )
      CALL ReadOpacityTypeCHDF( OpacityTable % Scat_nIso , group_id )
      CALL CloseGroupHDF( group_id )

      CALL CloseFileHDF( file_id )
    
    ELSE 
      WRITE(*,*) "ERROR!"
      WRITE(*,*) "EquationOfStateTable is not consistent with OpacityTable!"
      STOP
    END IF

    END ASSOCIATE ! nPointsTS_EOS

  END SUBROUTINE ReadOpacityTableHDF


  SUBROUTINE ReadOpacityTypeEmAbHDF( EmAb, group_id )

    TYPE(OpacityTypeEmAb),INTENT(inout)                 :: EmAb
    INTEGER(HID_T), INTENT(in)                       :: group_id

    INTEGER(HSIZE_T), DIMENSION(1)                   :: datasize1d
    INTEGER(HSIZE_T), DIMENSION(4)                   :: datasize4d
    INTEGER                                          :: i
    INTEGER, DIMENSION(1)                            :: buffer
    REAL(dp), DIMENSION(1)                           :: bufferReal
    INTEGER(HID_T)                                   :: subgroup_id

    datasize1d(1) = 1
    CALL ReadHDF( "nOpacities", buffer, group_id, datasize1d )
    EmAb % nOpacities = buffer(1)

    datasize1d = buffer(1)
    CALL ReadHDF( "Offsets", EmAb % Offsets, group_id, datasize1d )

    Call ReadHDF( "Names", EmAb % Names, group_id, datasize1d )

    Call ReadHDF( "Units", EmAb % Units, group_id, datasize1d )

    datasize1d(1) = 4
    CALL ReadHDF( "nPoints", EmAb % nPoints, group_id, datasize1d )

    datasize4d = EmAb % nPoints

    CALL OpenGroupHDF( "Absorptivity", .false., group_id, subgroup_id )

    DO i = 1, EmAb % nOpacities

      CALL ReadHDF &
             ( EmAb % Names(i), &
               EmAb % Absorptivity(i) % Values, &
               subgroup_id, datasize4d )

    END DO ! nOpacities

    CALL CloseGroupHDF( subgroup_id )

  END SUBROUTINE ReadOpacityTypeEmAbHDF


  SUBROUTINE ReadOpacityTypeScatHDF( Scat_Iso, group_id )

    TYPE(OpacityTypeScat),INTENT(inout)                 :: Scat_Iso
    INTEGER(HID_T), INTENT(in)                       :: group_id

    INTEGER(HSIZE_T), DIMENSION(1)                   :: datasize1d
    INTEGER(HSIZE_T), DIMENSION(2)                   :: datasize2d
    INTEGER(HSIZE_T), DIMENSION(4)                   :: datasize4d
    INTEGER(HSIZE_T), DIMENSION(5)                   :: datasize5d
    INTEGER                                          :: i
    INTEGER, DIMENSION(1)                            :: buffer
    REAL(dp), DIMENSION(1)                           :: bufferReal
    INTEGER(HID_T)                                   :: subgroup_id

    datasize1d(1) = 1
    CALL ReadHDF( "nOpacities", buffer, group_id, datasize1d )
    Scat_Iso % nOpacities = buffer(1)

    CALL ReadHDF( "nMoments", buffer, group_id, datasize1d )
    Scat_Iso % nMoments   = buffer(1)

    datasize1d = buffer(1)
    Call ReadHDF( "Names", Scat_Iso % Names, group_id, datasize1d )

    Call ReadHDF( "Units", Scat_Iso % Units, group_id, datasize1d )

    CALL ReadHDF( "Species", Scat_Iso % Species, group_id, datasize1d )

    datasize1d(1) = 4
    CALL ReadHDF( "nPoints", Scat_Iso % nPoints, group_id, datasize1d )

    datasize2d = (/Scat_Iso % nOpacities, Scat_Iso % nMoments/)
    CALL ReadHDF( "Offsets", Scat_Iso % Offsets, group_id, datasize2d )

    datasize5d(1:4) = Scat_Iso % nPoints
    datasize5d(5) = Scat_Iso % nMoments

    CALL OpenGroupHDF( "Kernel", .false., group_id, subgroup_id )

    DO i = 1, Scat_Iso % nOpacities

      CALL ReadHDF &
             ( Scat_Iso % Names(i), &
               Scat_Iso % Kernel(i) % Values, &
               subgroup_id, datasize5d )

    END DO ! nOpacities

    CALL CloseGroupHDF( subgroup_id )

  END SUBROUTINE ReadOpacityTypeScatHDF


  SUBROUTINE ReadOpacityTypeCHDF( EmAb, group_id )

    TYPE(OpacityTypeC),INTENT(inout)                 :: EmAb
    INTEGER(HID_T), INTENT(in)                       :: group_id

  END SUBROUTINE ReadOpacityTypeCHDF


  SUBROUTINE ReadGridHDF( Grid, group_id )

    TYPE(GridType), INTENT(inout)               :: Grid
    INTEGER(HID_T), INTENT(in)                  :: group_id

    INTEGER(HSIZE_T), DIMENSION(1)              :: datasize1d
    INTEGER, DIMENSION(1)                       :: buffer
    CHARACTER(LEN=32), DIMENSION(1)             :: buffer_string

    datasize1d(1) = 1
    Call ReadHDF( "Name", buffer_string, group_id, datasize1d )
    Grid % Name = buffer_string(1)

    Call ReadHDF( "Unit", buffer_string, group_id, datasize1d )
    Grid % Unit = buffer_string(1)

    CALL ReadHDF( "nPoints", buffer, group_id, datasize1d )
    Grid % nPoints = buffer(1)

    CALL ReadHDF( "LogInterp", buffer, group_id, datasize1d )
    Grid % LogInterp = buffer(1)
 
    datasize1d = Grid % nPoints
    CALL ReadHDF( "Values", Grid % Values, &
                              group_id, datasize1d )

    Grid % minValue = MINVAL( Grid % Values )
    
    Grid % maxValue = MAXVAL( Grid % Values )

  END SUBROUTINE ReadGridHDF

END MODULE wlOpacityTableIOModuleHDF
