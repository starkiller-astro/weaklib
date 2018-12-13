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

  PUBLIC WriteOpacityTableHDF_New
  PUBLIC ReadOpacityTableHDF_New

CONTAINS

  SUBROUTINE WriteOpacityTableHDF_New &
    ( OpacityTable, FileName, WriteOpacity_EmAb_Option, &
      WriteOpacity_Iso_Option, WriteOpacity_NES_Option, &
      WriteOpacity_Pair_Option )
 
    TYPE(OpacityTableType), INTENT(inout)        :: OpacityTable
    CHARACTER(len=*),       INTENT(in)           :: FileName
    LOGICAL,                INTENT(in), OPTIONAL :: WriteOpacity_EmAb_Option
    LOGICAL,                INTENT(in), OPTIONAL :: WriteOpacity_Iso_Option
    LOGICAL,                INTENT(in), OPTIONAL :: WriteOpacity_NES_Option
    LOGICAL,                INTENT(in), OPTIONAL :: WriteOpacity_Pair_Option

    LOGICAL           :: WriteOpacity_EmAb
    LOGICAL           :: WriteOpacity_Iso 
    LOGICAL           :: WriteOpacity_NES 
    LOGICAL           :: WriteOpacity_Pair
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

    IF( WriteOpacity_Iso ) THEN

       IF( .NOT. ALLOCATED( OpacityTable % Scat_Iso % Names ) )THEN

        ! --- Insert Appropriate Reaction ---
        WRITE(*,'(A4,A)') &
          '', 'OpacityTable % Scat_Iso not allocated.  Write Skipped.'

      ELSE

        CALL OpenGroupHDF &
               ( "Scat_Iso_Kernel", .true., file_id, group_id )
        CALL WriteOpacityTableHDF_Scat( OpacityTable % Scat_Iso, group_id )
        CALL CloseGroupHDF( group_id )

      END IF

    END IF

    IF( WriteOpacity_NES .or. WriteOpacity_Pair) THEN

      CALL OpenGroupHDF( "EtaGrid", .true., file_id, group_id )
      CALL WriteGridHDF( OpacityTable % EtaGrid, group_id )
      CALL CloseGroupHDF( group_id )

      IF( WriteOpacity_NES ) THEN

        IF( .NOT. ALLOCATED( OpacityTable % Scat_NES % Names ) )THEN
  
          ! --- Insert Appropriate Reaction ---
          WRITE(*,'(A4,A)') &
            '', 'OpacityTable % Scat_NES not allocated.  Write Skipped.'
  
        ELSE
  
          CALL OpenGroupHDF &
                 ( "NES_MomentsComponents", .true., file_id, group_id )
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
                 ( "Pair_MomentsComponents", .true., file_id, group_id )
          CALL WriteOpacityTableHDF_Scat( OpacityTable % Scat_Pair, group_id )
          CALL CloseGroupHDF( group_id )

        END IF

      END IF

    END IF           

    CALL CloseFileHDF( file_id )

  END SUBROUTINE WriteOpacityTableHDF_New


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

    datasize1dtemp(1) = Scat % nOpacities
    CALL WriteHDF&
         ( "Units", Scat % Units, group_id, datasize1dtemp )

    datasize2d = (/Scat % nOpacities, Scat % nMoments/)
    CALL WriteHDF&
         ( "Offsets", Scat % Offsets, group_id, datasize2d )

    datasize5d(1:4) = Scat % nPoints
    datasize5d(5) = Scat % nMoments

    DO i = 1, Scat % nOpacities
     CALL WriteHDF&
        ( Scat % Names(i), Scat % Kernel(i) % Values(:,:,:,:,:),&
                            group_id, datasize5d )
    END DO

  END SUBROUTINE WriteOpacityTableHDF_Scat


  SUBROUTINE WriteOpacityTableTypeCHDF( scattn, group_id )

    TYPE(OpacityTypeC), INTENT(in)              :: scattn
    INTEGER(HID_T), INTENT(in)                  :: group_id

  END SUBROUTINE WriteOpacityTableTypeCHDF


  SUBROUTINE ReadOpacityTableHDF_New &
    ( OpacityTable, FileName, ReadOpacity_EmAb_Option, &
      ReadOpacity_Iso_Option, ReadOpacity_NES_Option,  &
      ReadOpacity_Pair_Option )
 
    LOGICAL,                INTENT(in),   OPTIONAL :: ReadOpacity_EmAb_Option
    LOGICAL,                INTENT(in),   OPTIONAL :: ReadOpacity_Iso_Option
    LOGICAL,                INTENT(in),   OPTIONAL :: ReadOpacity_NES_Option
    LOGICAL,                INTENT(in),   OPTIONAL :: ReadOpacity_Pair_Option
    CHARACTER(len=*),       INTENT(in)             :: FileName
    TYPE(OpacityTableType), INTENT(inout)          :: OpacityTable

    LOGICAL            :: ReadOpacity_EmAb
    LOGICAL            :: ReadOpacity_Iso
    LOGICAL            :: ReadOpacity_NES
    LOGICAL            :: ReadOpacity_Pair
    INTEGER            :: nPointsE
    INTEGER            :: nPointsEta
    INTEGER            :: nOpacA
    INTEGER            :: nOpacB, nMomB
    INTEGER            :: nOpacB_NES, nMomB_NES
    INTEGER            :: nOpacB_TP, nMomB_TP
    INTEGER            :: hdfreadErr
    INTEGER            :: buffer(1)
    INTEGER(HID_T)     :: file_id
    INTEGER(HID_T)     :: group_id
    INTEGER(HID_T)     :: subgroup_id
    INTEGER(HSIZE_T)   :: datasize1d(1)
    INTEGER(HSIZE_T)   :: datasize2d(2)
    INTEGER(HSIZE_T)   :: datasize4d(4)
    INTEGER(HSIZE_T)   :: datasize5d(5)
    CHARACTER(LEN=32)  :: buffer_string(1)

    hdfreadErr = 0

    IF( PRESENT( ReadOpacity_EmAb_Option ) )THEN
      ReadOpacity_EmAb = ReadOpacity_EmAb_Option
      nOpacA = 2
    ELSE
      ReadOpacity_EmAb = .FALSE.
      nOpacA = 0
    END IF

    IF( PRESENT( ReadOpacity_Iso_Option ) )THEN
      ReadOpacity_Iso = ReadOpacity_Iso_Option
      nOpacB = 2
      nMomB  = 2
    ELSE
      ReadOpacity_Iso = .FALSE.
      nOpacB = 0
      nMomB  = 0
    END IF

    IF( PRESENT( ReadOpacity_NES_Option ) )THEN
      ReadOpacity_NES = ReadOpacity_NES_Option
      nOpacB_NES = 2
      nMomB_NES  = 2
    ELSE
      ReadOpacity_NES = .FALSE.
      nOpacB_NES = 0
      nMomB_NES  = 0
    END IF

       IF( PRESENT( ReadOpacity_Pair_Option ) )THEN
      ReadOpacity_Pair = ReadOpacity_Pair_Option
      nOpacB_TP = 2
      nMomB_TP  = 2
    ELSE
      ReadOpacity_Pair = .FALSE.
      nOpacB_TP = 0
      nMomB_TP  = 0
    END IF

    WRITE(*,*) "           File in"
    WRITE(*,*) " Reading ", FileName, " hdf5 file ... "

    CALL OpenFileHDF( FileName, .false., file_id )

    CALL OpenGroupHDF( "EnergyGrid", .false., file_id, group_id )
    CALL ReadHDF( "nPoints", buffer, group_id, datasize1d )
    nPointsE = buffer(1)
    CALL CloseGroupHDF( group_id )

    IF( ReadOpacity_NES .OR. ReadOpacity_Pair ) THEN
      CALL OpenGroupHDF( "EtaGrid", .false., file_id, group_id )
      CALL ReadHDF( "nPoints", buffer, group_id, datasize1d )
      nPointsEta = buffer(1)
      CALL CloseGroupHDF( group_id )
    ELSE
      nPointsEta = 0
    END IF

    CALL AllocateOpacityTable &
           ( OpacityTable, nOpacA, nOpacB, nMomB, nOpacB_NES, nMomB_NES, &
             nOpacB_TP, nMomB_TP, 0, 0, nPointsE, nPointsEta )
   
    IF( hdfreadErr .NE. 0 ) THEN
      WRITE(*,*) "ERROR!"
      WRITE(*,*) "EquationOfStateTable is not consistent with OpacityTable!"
      CALL CloseFileHDF( file_id )
      STOP
    END IF

    CALL OpenGroupHDF( "EnergyGrid", .false., file_id, group_id )
    CALL ReadGridHDF( OpacityTable % EnergyGrid, group_id )
    CALL CloseGroupHDF( group_id )

    IF( ReadOpacity_EmAb )THEN
      WRITE(*,*) "Now read-in EmAb table"

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
         ( TRIM( OpacityTable % EmAb % Names(2) ), &
            OpacityTable % EmAb % Absorptivity(2) % Values, &
            group_id, datasize4d )
      CALL CloseGroupHDF( group_id )

    END IF ! ReadOpacity_EmAb

    IF( ReadOpacity_NES ) THEN
      WRITE(*,*) "Now read-in NES table"

      CALL OpenGroupHDF( "EtaGrid", .false., file_id, group_id )
      CALL ReadGridHDF( OpacityTable % EtaGrid, group_id )
      CALL CloseGroupHDF( group_id )    

      CALL OpenGroupHDF &
             ( "NES_MomentsComponents", .false., file_id, group_id )

      datasize2d(1:2) = 2
      CALL ReadHDF &
             ( "Offsets", OpacityTable % Scat_NES % Offsets, group_id, datasize2d )

      CALL ReadHDF &
             ( "Units",   OpacityTable % Scat_NES % Units,   group_id, datasize2d )

      datasize5d(1:2) = nPointsE
      datasize4d(3)   = OpacityTable % TS % nPoints(2) !! fix me
      datasize5d(4)   = OpacityTable % EtaGrid % nPoints
      datasize5d(5)   = nMomB_NES

      OpacityTable % Scat_NES % Names(1) = "NES Kernel Moment H0i H0ii";
      CALL ReadHDF &
             ( TRIM( OpacityTable % Scat_NES % Names(1) ), &
               OpacityTable % Scat_NES % Kernel(1) % Values, &
               group_id, datasize5d )

      OpacityTable % Scat_NES % Names(2) = "NES Kernel Moment H1i H1ii";
      CALL ReadHDF &
         ( TRIM( OpacityTable % Scat_NES % Names(2) ), &
            OpacityTable % Scat_NES % Kernel(2) % Values, &
            group_id, datasize5d )
      CALL CloseGroupHDF( group_id )

    END IF ! ReadOpacity_NES

    CALL CloseFileHDF( file_id )

  END SUBROUTINE ReadOpacityTableHDF_New


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
