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
!
!    WeakLib ver:  
!                  08/02/21 	Added Scat_Brem		Vassilios Mewes
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
    OpacityTypeScat
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
  USE wlThermoStateModule, ONLY:  &
    ThermoStateType, &
    AllocateThermoState, &
    DeAllocateThermoState
  USE HDF5

  IMPLICIT NONE
  PRIVATE

  INTEGER :: hdferr

  PUBLIC WriteOpacityTableHDF
  PUBLIC ReadOpacityTableHDF

CONTAINS

  SUBROUTINE WriteOpacityTableHDF &
    ( OpacityTable, FileName, WriteOpacity_EmAb_Option, &
      WriteOpacity_Iso_Option, WriteOpacity_NES_Option, &
      WriteOpacity_Pair_Option, WriteOpacity_Brem_Option )
 
    TYPE(OpacityTableType), INTENT(inout)        :: OpacityTable
    CHARACTER(len=*),       INTENT(in)           :: FileName
    LOGICAL,                INTENT(in), OPTIONAL :: WriteOpacity_EmAb_Option
    LOGICAL,                INTENT(in), OPTIONAL :: WriteOpacity_Iso_Option
    LOGICAL,                INTENT(in), OPTIONAL :: WriteOpacity_NES_Option
    LOGICAL,                INTENT(in), OPTIONAL :: WriteOpacity_Pair_Option
    LOGICAL,                INTENT(in), OPTIONAL :: WriteOpacity_Brem_Option

    LOGICAL           :: WriteOpacity_EmAb
    LOGICAL           :: WriteOpacity_Iso 
    LOGICAL           :: WriteOpacity_NES 
    LOGICAL           :: WriteOpacity_Pair
    LOGICAL           :: WriteOpacity_Brem
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
  
    IF( PRESENT( WriteOpacity_Iso_Option ) )THEN
      WriteOpacity_Iso = WriteOpacity_Iso_Option
    ELSE
      WriteOpacity_Iso = .FALSE.
    END IF

    IF( PRESENT( WriteOpacity_NES_Option ) )THEN
      WriteOpacity_NES = WriteOpacity_NES_Option
    ELSE
      WriteOpacity_NES = .FALSE.
    END IF

    IF( PRESENT( WriteOpacity_Pair_Option ) )THEN
      WriteOpacity_Pair = WriteOpacity_Pair_Option
    ELSE
      WriteOpacity_Pair = .FALSE.
    END IF

    IF( PRESENT( WriteOpacity_Brem_Option ) )THEN
      WriteOpacity_Brem = WriteOpacity_Brem_Option
    ELSE
      WriteOpacity_Brem = .FALSE.
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

    WRITE(*,*) "Writting out EmAb"
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

    IF( WriteOpacity_Iso ) THEN

    WRITE(*,*) "Writting out Iso"
       IF( .NOT. ALLOCATED( OpacityTable % Scat_Iso % Names ) )THEN

        ! --- Insert Appropriate Reaction ---
        WRITE(*,'(A4,A)') &
          '', 'OpacityTable % Scat_Iso not allocated.  Write Skipped.'

      ELSE

        CALL OpenGroupHDF &
               ( "Scat_Iso_Kernels", .true., file_id, group_id )
        CALL WriteOpacityTableHDF_Scat( OpacityTable % Scat_Iso, group_id )
        CALL CloseGroupHDF( group_id )

      END IF

    END IF

    IF( WriteOpacity_NES .or. WriteOpacity_Pair) THEN

      CALL OpenGroupHDF( "EtaGrid", .true., file_id, group_id )
      CALL WriteGridHDF( OpacityTable % EtaGrid, group_id )
      CALL CloseGroupHDF( group_id )

      IF( WriteOpacity_NES ) THEN
      WRITE(*,*) "Writting out NES"

        IF( .NOT. ALLOCATED( OpacityTable % Scat_NES % Names ) )THEN
  
          ! --- Insert Appropriate Reaction ---
          WRITE(*,'(A4,A)') &
            '', 'OpacityTable % Scat_NES not allocated.  Write Skipped.'
  
        ELSE
  
          CALL OpenGroupHDF &
                 ( "Scat_NES_Kernels", .true., file_id, group_id )
          CALL WriteOpacityTableHDF_Scat( OpacityTable % Scat_NES, group_id )
          CALL CloseGroupHDF( group_id )
  
        END IF

      END IF

      IF( WriteOpacity_Pair ) THEN

        IF( .NOT. ALLOCATED( OpacityTable % Scat_Pair % Names ) )THEN

          ! --- Insert Appropriate Reaction ---
          WRITE(*,'(A4,A)') &
            '', 'OpacityTable % Scat_Pair not allocated.  Write Skipped.'

        ELSE

          CALL OpenGroupHDF &
                 ( "Scat_Pair_Kernels", .true., file_id, group_id )
          CALL WriteOpacityTableHDF_Scat( OpacityTable % Scat_Pair, group_id )
          CALL CloseGroupHDF( group_id )

        END IF

      END IF

    END IF           

    IF( WriteOpacity_Brem ) THEN

      IF( .NOT. ALLOCATED( OpacityTable % Scat_Brem % Names ) )THEN

        ! --- Insert Appropriate Reaction ---
        WRITE(*,'(A4,A)') &
                '', 'OpacityTable % Scat_Brem not allocated.  Write Skipped.'

      ELSE

        CALL OpenGroupHDF &
                 ( "Scat_Brem_Kernels", .true., file_id, group_id )
        CALL WriteOpacityTableHDF_Scat( OpacityTable % Scat_Brem, group_id )
        CALL CloseGroupHDF( group_id )

      END IF

    END IF

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
    INTEGER :: ii

    INTEGER, DIMENSION(1)             :: tempInteger

    datasize1d = 1
    tempInteger(1) = EmAb % nOpacities
    CALL WriteHDF( "nOpacities", tempInteger, group_id, datasize1d )

    datasize1d = EmAb % nOpacities
    CALL WriteHDF &
           ( "Units", EmAb % Units, group_id, datasize1d ) 

    CALL WriteHDF &
           ( "Offsets", EmAb % Offsets, group_id, datasize1d )

    datasize4d = EmAb % nPoints

    DO ii = 1, EmAb % nOpacities

      CALL WriteHDF &
             ( TRIM( EmAb % Names(ii) ), &
               EmAb % Opacity(ii) % Values(:,:,:,:), group_id, datasize4d )

    END DO
  
  END SUBROUTINE WriteOpacityTableHDF_EmAb


  SUBROUTINE WriteOpacityTableHDF_Scat( Scat , group_id )

    TYPE(OpacityTypeScat), INTENT(in)    :: Scat
    INTEGER(HID_T),        INTENT(in)    :: group_id

    INTEGER(HSIZE_T)                     :: datasize1d(1)
    INTEGER(HSIZE_T)                     :: datasize2d(2)
    INTEGER(HSIZE_T)                     :: datasize4d(4)
    INTEGER(HSIZE_T)                     :: datasize5d(5)
    INTEGER                              :: i
    INTEGER                              :: buffer(1)

    CHARACTER(LEN=32), DIMENSION(1)             :: tempString
    INTEGER, DIMENSION(1)                       :: tempInteger
    REAL(dp), DIMENSION(1)                      :: tempReal
    INTEGER(HSIZE_T), DIMENSION(1)              :: datasize1dtemp

    datasize1d = 1
    tempInteger(1) = Scat % nOpacities
    CALL WriteHDF( "nOpacities", tempInteger, group_id, datasize1d )

    tempInteger(1) = Scat % nMoments
    CALL WriteHDF( "nMoments", tempInteger, group_id, datasize1d )

    datasize1dtemp(1) = Scat % nOpacities
    CALL WriteHDF&
         ( "Units", Scat % Units, group_id, datasize1dtemp )

    datasize2d = (/Scat % nOpacities, Scat % nMoments/)
    CALL WriteHDF&
         ( "Offsets", Scat % Offsets, group_id, datasize2d )

    datasize5d(1:5) = Scat % nPoints

    DO i = 1, Scat % nOpacities
     CALL WriteHDF&
        ( Scat % Names(i), Scat % Kernel(i) % Values(:,:,:,:,:),&
                            group_id, datasize5d )
    END DO

  END SUBROUTINE WriteOpacityTableHDF_Scat

  SUBROUTINE ReadOpacityTableHDF &
    ( OpacityTable, FileName_EmAb_Option, FileName_Iso_Option, &
      FileName_NES_Option, FileName_Pair_Option, FileName_Brem_Option, &
      EquationOfStateTableName_Option, Verbose_Option )
 
    TYPE(OpacityTableType), INTENT(inout)          :: OpacityTable
    CHARACTER(len=*),       INTENT(in),   OPTIONAL :: FileName_EmAb_Option
    CHARACTER(len=*),       INTENT(in),   OPTIONAL :: FileName_Iso_Option
    CHARACTER(len=*),       INTENT(in),   OPTIONAL :: FileName_NES_Option
    CHARACTER(len=*),       INTENT(in),   OPTIONAL :: FileName_Pair_Option
    CHARACTER(len=*),       INTENT(in),   OPTIONAL :: FileName_Brem_Option
    CHARACTER(LEN=*),       INTENT(in),   OPTIONAL :: EquationOfStateTableName_Option
    LOGICAL,                INTENT(in),   OPTIONAL :: Verbose_Option

    INTEGER, PARAMETER :: iEmAb = 1
    INTEGER, PARAMETER :: iIso  = 2
    INTEGER, PARAMETER :: iNES  = 3
    INTEGER, PARAMETER :: iPair = 4
    INTEGER, PARAMETER :: iBrem = 5

    LOGICAL            :: ReadOpacity(5)
    LOGICAL            :: Verbose
    CHARACTER(128)     :: FileName(5)
    CHARACTER(128)     :: EquationOfStateTableName
    INTEGER            :: iOp
    INTEGER            :: nPointsE
    INTEGER            :: nPointsEta
    INTEGER            :: nPointsTS(3)
    INTEGER            :: nOpac_EmAb
    INTEGER            :: nOpac_Iso
    INTEGER            :: nMom_Iso
    INTEGER            :: nOpac_NES
    INTEGER            :: nMom_NES
    INTEGER            :: nOpac_Pair
    INTEGER            :: nMom_Pair
    INTEGER            :: nOpac_Brem
    INTEGER            :: nMom_Brem
    INTEGER            :: buffer(1)
    INTEGER(HID_T)     :: file_id
    INTEGER(HID_T)     :: group_id
    INTEGER(HSIZE_T)   :: datasize1d(1)
    INTEGER(HSIZE_T)   :: datasize2d(2)
    INTEGER(HSIZE_T)   :: datasize3d(3)
    INTEGER(HSIZE_T)   :: datasize4d(4)
    INTEGER(HSIZE_T)   :: datasize5d(5)

    TYPE(ThermoStateType) :: TS

    IF( PRESENT( EquationOfStateTableName_Option ) &
        .AND. ( LEN( EquationOfStateTableName_Option ) > 1 ) )THEN
       EquationOfStateTableName = TRIM( EquationOfStateTableName_Option )
    ELSE
       EquationOfStateTableName = 'EquationOfStateTable.h5'
    END IF

    IF( PRESENT( FileName_EmAb_Option ) &
        .AND. ( LEN( FileName_EmAb_Option ) > 1 ) )THEN
      ReadOpacity(iEmAb) = .TRUE.
      FileName   (iEmAb) = TRIM( FileName_EmAb_Option )
    ELSE
      ReadOpacity(iEmAb) = .FALSE.
      nOpac_EmAb = 0
    END IF

    IF( PRESENT( FileName_Iso_Option ) &
        .AND. ( LEN( FileName_Iso_Option ) > 1 ) )THEN
      ReadOpacity(iIso) = .TRUE.
      FileName   (iIso) = TRIM( FileName_Iso_Option )
    ELSE
      ReadOpacity(iIso) = .FALSE.
      nOpac_Iso = 0
      nMom_Iso  = 0
    END IF

    IF( PRESENT( FileName_NES_Option ) &
        .AND. ( LEN( FileName_NES_Option ) > 1 ) )THEN
      ReadOpacity(iNES) = .TRUE.
      FileName   (iNES) = TRIM( FileName_NES_Option )
    ELSE
      ReadOpacity(iNES) = .FALSE.
      nOpac_NES = 0
      nMom_NES  = 0
    END IF

    IF( PRESENT( FileName_Pair_Option ) &
        .AND. ( LEN( FileName_Pair_Option ) > 1 ) )THEN
      ReadOpacity(iPair) = .TRUE.
      FileName   (iPair) = TRIM( FileName_Pair_Option )
    ELSE
      ReadOpacity(iPair) = .FALSE.
      nOpac_Pair = 0
      nMom_Pair  = 0
    END IF

    IF( PRESENT( FileName_Brem_Option ) &
        .AND. ( LEN( FileName_Brem_Option ) > 1 ) )THEN
      ReadOpacity(iBrem) = .TRUE.
      FileName   (iBrem) = TRIM( FileName_Brem_Option )
    ELSE
      ReadOpacity(iBrem) = .FALSE.
      nOpac_Brem = 0
      nMom_Brem  = 0
    END IF

    IF( PRESENT( Verbose_Option ) )THEN
      Verbose = Verbose_Option
    ELSE
      Verbose = .FALSE.
    END IF

    IF( Verbose )THEN
      WRITE(*,*)
      WRITE(*,'(A4,A)') '', 'ReadOpacityTableHDF'
      WRITE(*,*)
      DO iOp = 1, 4
        IF( ReadOpacity(iOp) ) WRITE(*,'(A6,A)') '', TRIM( FileName(iOp) )
      END DO
      WRITE(*,*)
    END IF

    IF( .NOT. ANY( ReadOpacity ) )THEN
      WRITE(*,'(A4,A)') '', 'ERROR: No Opacity Table Provided. Returning'
      RETURN
    END IF

    ! --- Get Number of Energy Points ---

    nPointsE = 0
    DO iOp = iEmAb, iBrem

      IF( ReadOpacity(iOp) )THEN

        CALL OpenFileHDF( FileName(iOp), .FALSE., file_id )

        CALL OpenGroupHDF( "EnergyGrid", .FALSE., file_id, group_id )

        CALL ReadHDF( "nPoints", buffer, group_id, datasize1d )

        CALL CloseGroupHDF( group_id )

        CALL CloseFileHDF( file_id )

        nPointsE = buffer(1)

        EXIT

      END IF

    END DO

    ! --- Get Number of Eta (Chem_e/kT) Points ---

    nPointsEta = 0
    DO iOp = iNES, iPair

      IF( ReadOpacity(iOp) )THEN

        CALL OpenFileHDF( FileName(iOp), .FALSE., file_id )

        CALL OpenGroupHDF( "EtaGrid",    .FALSE., file_id, group_id )

        CALL ReadHDF( "nPoints", buffer, group_id, datasize1d )

        CALL CloseGroupHDF( group_id )

        CALL CloseFileHDF( file_id )

        nPointsEta = buffer(1)

        EXIT

      END IF

    END DO

    ! --- Get Number of ThermoState Points ---

    DO iOp = iEmAb, iBrem

      IF( ReadOpacity(iOp) )THEN

        CALL OpenFileHDF( FileName(iOp), .FALSE., file_id )

        CALL OpenGroupHDF( "ThermoState",.FALSE., file_id, group_id )

        CALL ReadHDF( "Dimensions", nPointsTS, group_id, datasize3d )

        CALL AllocateThermoState( TS, nPointsTS )

        CALL ReadThermoStateHDF( TS, file_id )

        CALL CloseFileHDF( file_id )
 
        EXIT

      END IF

    END DO

    ! --- Get Number of Opacities and Moments ---
    IF(ReadOpacity(iEmAb)) THEN

      buffer(1)  = 2 ! for old opacity table safe

      CALL OpenFileHDF( FileName(iEmAb), .FALSE., file_id )

      CALL OpenGroupHDF( "EmAb_CorrectedAbsorption", .FALSE., file_id, group_id )

      CALL ReadHDF( "nOpacities", buffer, group_id, datasize1d )

      CALL CloseGroupHDF( group_id )

      CALL CloseFileHDF( file_id )

      nOpac_EmAb = buffer(1)

    ENDIF

    IF (ReadOpacity(iIso)) THEN

      buffer(1)  = 2 ! for old opacity table safe

      CALL OpenFileHDF( FileName(iIso), .FALSE., file_id )

      CALL OpenGroupHDF( "Scat_Iso_Kernels", .FALSE., file_id, group_id )

      CALL ReadHDF( "nOpacities", buffer, group_id, datasize1d )

      nOpac_Iso = buffer(1)

      buffer(1)  = 2 ! for old opacity table safe

      CALL ReadHDF( "nMoments", buffer, group_id, datasize1d )

      nMom_Iso = buffer(1)

      CALL CloseGroupHDF( group_id )

      CALL CloseFileHDF( file_id )

    END IF

    IF (ReadOpacity(iNES)) THEN

      buffer(1)  = 1 ! for old opacity table safe

      CALL OpenFileHDF( FileName(iNES), .FALSE., file_id )

      CALL OpenGroupHDF( "Scat_NES_Kernels", .FALSE., file_id, group_id )

      CALL ReadHDF( "nOpacities", buffer, group_id, datasize1d )

      nOpac_NES = buffer(1)

      buffer(1)  = 4 ! for old opacity table safe

      CALL ReadHDF( "nMoments", buffer, group_id, datasize1d )

      nMom_NES = buffer(1)

      CALL CloseGroupHDF( group_id )

      CALL CloseFileHDF( file_id )

    END IF

    IF (ReadOpacity(iPair)) THEN

      buffer(1)  = 1 ! for old opacity table safe

      CALL OpenFileHDF( FileName(iPair), .FALSE., file_id )

      CALL OpenGroupHDF( "Scat_Pair_Kernels", .FALSE., file_id, group_id )

      CALL ReadHDF( "nOpacities", buffer, group_id, datasize1d )

      nOpac_Pair = buffer(1)

      buffer(1)  = 4 ! for old opacity table safe

      CALL ReadHDF( "nMoments", buffer, group_id, datasize1d )

      nMom_Pair = buffer(1)

      CALL CloseGroupHDF( group_id )

      CALL CloseFileHDF( file_id )

    END IF

    IF (ReadOpacity(iBrem)) THEN

      CALL OpenFileHDF( FileName(iBrem), .FALSE., file_id )

      CALL OpenGroupHDF( "Scat_Brem_Kernels", .FALSE., file_id, group_id )

      CALL ReadHDF( "nOpacities", buffer, group_id, datasize1d )

      nOpac_Brem = buffer(1)

      CALL ReadHDF( "nMoments", buffer, group_id, datasize1d )

      nMom_Brem = buffer(1)

      CALL CloseGroupHDF( group_id )

      CALL CloseFileHDF( file_id )

    END IF

    CALL AllocateOpacityTable &
           ( OpacityTable, nOpac_EmAb, nOpac_Iso, nMom_Iso, nOpac_NES, &
             nMom_NES, nOpac_Pair, nMom_Pair, nOpac_Brem, nMom_Brem, &
             nPointsE, nPointsEta, &
             EquationOfStateTableName_Option = EquationOfStateTableName, &
             OpacityThermoState_Option = TS, &
             Verbose_Option = Verbose )

    CALL DeAllocateThermoState( TS )

    ! --- Read Energy Grid ---

    DO iOp = iEmAb, iBrem

      IF( ReadOpacity(iOp) )THEN

        CALL OpenFileHDF( FileName(iOp), .FALSE., file_id )

        CALL OpenGroupHDF( "EnergyGrid", .FALSE., file_id, group_id )

        CALL ReadGridHDF( OpacityTable % EnergyGrid, group_id )

        CALL CloseGroupHDF( group_id )

        CALL CloseFileHDF( file_id )

        EXIT

      END IF

    END DO

    ! --- Read Eta (Chem_e/kT) Grid ---

    DO iOp = iNES, iPair

      IF( ReadOpacity(iOp) )THEN

        CALL OpenFileHDF( FileName(iOp), .FALSE., file_id )

        CALL OpenGroupHDF( "EtaGrid",    .FALSE., file_id, group_id )

        CALL ReadGridHDF( OpacityTable % EtaGrid, group_id )

        CALL CloseGroupHDF( group_id )

        CALL CloseFileHDF( file_id )

        EXIT

      END IF

    END DO

    ! --- Read Opacities ---

    IF( ReadOpacity(iEmAb) )THEN

      IF( Verbose )THEN

        WRITE(*,'(A6,A9,A)') '', 'Reading: ', TRIM( FileName(iEmAb) )

      END IF

      CALL OpenFileHDF( FileName(iEmAb), .FALSE., file_id )

      CALL OpenGroupHDF &
             ( "EmAb_CorrectedAbsorption", .FALSE., file_id, group_id )

      datasize1d(1) = nOpac_EmAb
      CALL ReadHDF &
             ( "Offsets", OpacityTable % EmAb % Offsets, group_id, datasize1d )

      CALL ReadHDF &
             ( "Units",   OpacityTable % EmAb % Units,   group_id, datasize1d )

      datasize4d(1)   = OpacityTable % EnergyGrid % nPoints
      datasize4d(2:4) = OpacityTable % TS % nPoints

      OpacityTable % EmAb % Names(1) = "Electron Neutrino"

      CALL ReadHDF &
             ( TRIM( OpacityTable % EmAb % Names(1) ), &
               OpacityTable % EmAb % Opacity(1) % Values, &
               group_id, datasize4d )

      OpacityTable % EmAb % Names(2) = "Electron Antineutrino"

      CALL ReadHDF &
         ( TRIM( OpacityTable % EmAb % Names(2) ), &
            OpacityTable % EmAb % Opacity(2) % Values, &
            group_id, datasize4d )

      CALL CloseGroupHDF( group_id )

      CALL CloseFileHDF( file_id )

    END IF


    IF( ReadOpacity(iIso) ) THEN

      IF( Verbose )THEN

        WRITE(*,'(A6,A9,A)') '', 'Reading: ', TRIM( FileName(iIso) )

      END IF

      CALL OpenFileHDF( FileName(iIso), .FALSE., file_id )

      CALL OpenGroupHDF &
             ( "Scat_Iso_Kernels", .FALSE., file_id, group_id )

      datasize2d(1) = OpacityTable % Scat_Iso % nOpacities
      datasize2d(2) = OpacityTable % Scat_Iso % nMoments

      CALL ReadHDF &
             ( "Offsets", OpacityTable % Scat_Iso % Offsets, &
               group_id, datasize2d )

      datasize1d = OpacityTable % Scat_Iso % nOpacities
      CALL ReadHDF &
             ( "Units",   OpacityTable % Scat_Iso % Units,   &
               group_id, datasize1d )

      datasize5d(1)   = OpacityTable % EnergyGrid % nPoints
      datasize5d(2)   = OpacityTable % Scat_Iso % nMoments
      datasize5d(3:5) = OpacityTable % TS % nPoints

      OpacityTable % Scat_Iso % Names(1) = "Electron Neutrino"

      CALL ReadHDF &
             ( TRIM( OpacityTable % Scat_Iso % Names(1) ), &
               OpacityTable % Scat_Iso % Kernel(1) % Values, &
               group_id, datasize5d )

      OpacityTable % Scat_Iso % Names(2) = "Electron Antineutrino"

      CALL ReadHDF &
             ( TRIM( OpacityTable % Scat_Iso % Names(2) ), &
               OpacityTable % Scat_Iso % Kernel(2) % Values, &
               group_id, datasize5d )

      CALL CloseGroupHDF( group_id )

      CALL CloseFileHDF( file_id )

    END IF


    IF( ReadOpacity(iNES) ) THEN

      IF( Verbose )THEN

        WRITE(*,'(A6,A9,A)') '', 'Reading: ', TRIM( FileName(iNES) )

      END IF

      CALL OpenFileHDF( FileName(iNES), .FALSE., file_id )

      CALL OpenGroupHDF &
             ( "Scat_NES_Kernels", .FALSE., file_id, group_id )

      datasize2d(1) = OpacityTable % Scat_NES % nOpacities
      datasize2d(2) = OpacityTable % Scat_NES % nMoments

      CALL ReadHDF &
             ( "Offsets", OpacityTable % Scat_NES % Offsets, &
               group_id, datasize2d )

      CALL ReadHDF &
             ( "Units",   OpacityTable % Scat_NES % Units,   &
               group_id, datasize2d )

      datasize5d(1:2) = OpacityTable % EnergyGrid % nPoints
      datasize5d(3)   = OpacityTable % Scat_NES % nMoments
      datasize5d(4)   = &
        OpacityTable % TS % nPoints(OpacityTable % TS % Indices % iT)
      datasize5d(5)   = OpacityTable % EtaGrid % nPoints

      OpacityTable % Scat_NES % Names(1) = "Kernels";

      CALL ReadHDF &
             ( TRIM( OpacityTable % Scat_NES % Names(1) ), &
               OpacityTable % Scat_NES % Kernel(1) % Values, &
               group_id, datasize5d )

      CALL CloseGroupHDF( group_id )

      CALL CloseFileHDF( file_id )

    END IF


    IF( ReadOpacity(iPair) ) THEN

      IF( Verbose )THEN

        WRITE(*,'(A6,A9,A)') '', 'Reading: ', TRIM( FileName(iPair) )

      END IF

      CALL OpenFileHDF( FileName(iPair), .FALSE., file_id )

      CALL OpenGroupHDF &
             ( "Scat_Pair_Kernels", .FALSE., file_id, group_id )

      datasize2d(1) = OpacityTable % Scat_Pair % nOpacities
      datasize2d(2) = OpacityTable % Scat_Pair % nMoments

      CALL ReadHDF &
             ( "Offsets", OpacityTable % Scat_Pair % Offsets, &
               group_id, datasize2d )

      CALL ReadHDF &
             ( "Units",   OpacityTable % Scat_Pair % Units,   &
               group_id, datasize2d )

      datasize5d(1:2) = OpacityTable % EnergyGrid % nPoints
      datasize5d(3)   = OpacityTable % Scat_Pair % nMoments
      datasize5d(4)   = &
        OpacityTable % TS % nPoints(OpacityTable % TS % Indices % iT)
      datasize5d(5)   = OpacityTable % EtaGrid % nPoints

      OpacityTable % Scat_Pair % Names(1) = "Kernels";

      CALL ReadHDF &
             ( TRIM( OpacityTable % Scat_Pair % Names(1) ), &
               OpacityTable % Scat_Pair % Kernel(1) % Values, &
               group_id, datasize5d )

      CALL CloseGroupHDF( group_id )

      CALL CloseFileHDF( file_id )

    END IF

    IF( ReadOpacity(iBrem) ) THEN

      IF( Verbose )THEN

        WRITE(*,'(A6,A9,A)') '', 'Reading: ', TRIM( FileName(iBrem) )

      END IF

      CALL OpenFileHDF( FileName(iBrem), .FALSE., file_id )

      CALL OpenGroupHDF &
             ( "Scat_Brem_Kernels", .FALSE., file_id, group_id )

      datasize2d(1) = OpacityTable % Scat_Brem % nOpacities
      datasize2d(2) = OpacityTable % Scat_Brem % nMoments

      CALL ReadHDF &
             ( "Offsets", OpacityTable % Scat_Brem % Offsets, &
               group_id, datasize2d )

      CALL ReadHDF &
             ( "Units",   OpacityTable % Scat_Brem % Units,   &
               group_id, datasize2d )

      datasize5d(1:2) = OpacityTable % EnergyGrid % nPoints
      datasize5d(3)   = OpacityTable % Scat_Brem % nMoments
      datasize5d(4)   = & 
        OpacityTable % TS % nPoints(OpacityTable % TS % Indices % iRho)
      datasize5d(5)   = &
        OpacityTable % TS % nPoints(OpacityTable % TS % Indices % iT)

      OpacityTable % Scat_Brem % Names(1) = "S_sigma";

      CALL ReadHDF &
             ( TRIM( OpacityTable % Scat_Brem % Names(1) ), &
               OpacityTable % Scat_Brem % Kernel(1) % Values, &
               group_id, datasize5d )

      CALL CloseGroupHDF( group_id )

      CALL CloseFileHDF( file_id )

    END IF

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

    CALL OpenGroupHDF( "Opacity", .false., group_id, subgroup_id )

    DO i = 1, EmAb % nOpacities

      CALL ReadHDF &
             ( EmAb % Names(i), &
               EmAb % Opacity(i) % Values, &
               subgroup_id, datasize4d )

    END DO ! nOpacities

    CALL CloseGroupHDF( subgroup_id )

  END SUBROUTINE ReadOpacityTypeEmAbHDF


  SUBROUTINE ReadOpacityTypeScatHDF( Scat, group_id )

    TYPE(OpacityTypeScat),INTENT(inout)              :: Scat
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
    Scat % nOpacities = buffer(1)

    CALL ReadHDF( "nMoments", buffer, group_id, datasize1d )
    Scat % nMoments   = buffer(1)

    datasize1d = buffer(1)
    Call ReadHDF( "Names", Scat % Names, group_id, datasize1d )

    Call ReadHDF( "Units", Scat % Units, group_id, datasize1d )

    datasize1d(1) = 4
    CALL ReadHDF( "nPoints", Scat % nPoints, group_id, datasize1d )

    datasize2d = (/Scat % nOpacities, Scat % nMoments/)
    CALL ReadHDF( "Offsets", Scat % Offsets, group_id, datasize2d )

    datasize5d(1:5) = Scat % nPoints

    CALL OpenGroupHDF( "Kernel", .false., group_id, subgroup_id )

    DO i = 1, Scat % nOpacities

      CALL ReadHDF &
             ( Scat % Names(i), &
               Scat % Kernel(i) % Values, &
               subgroup_id, datasize5d )

    END DO ! nOpacities

    CALL CloseGroupHDF( subgroup_id )

  END SUBROUTINE ReadOpacityTypeScatHDF

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
