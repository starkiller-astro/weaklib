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

  USE wlKindModule, ONLY:            &
    dp
  USE wlGridModule, ONLY:            &
    GridType
  USE wlOpacityTableModule, ONLY:    &
    OpacityTableType,                &
    AllocateOpacityTable
  USE wlOpacityFieldsModule, ONLY:   &
    OpacityTypeEmAb,                 &
    OpacityTypeScat,                 &
    OpacityTypeScatIso,              &
    OpacityTypeScatNES,              &
    OpacityTypeScatNNS
  USE wlIOModuleHDF, ONLY:           &
    ReadHDF,                         &
    WriteHDF,                        &
    WriteGroupAttributeHDF_string,   &
    WriteVersionAttribute,           &
    OpenFileHDF,                     &
    CloseFileHDF,                    &
    OpenGroupHDF,                    &
    CloseGroupHDF,                   &
    WriteThermoStateHDF,             &
    ReadThermoStateHDF
  USE wlEquationOfStateTableModule
  USE wlThermoStateModule, ONLY:     &
    ThermoStateType,                 &
    AllocateThermoState,             &
    DeAllocateThermoState
  USE HDF5
  USE wlParallelModule, ONLY:        &
    myid, ierr

  IMPLICIT NONE
  PRIVATE

  INTEGER :: hdferr

  PUBLIC WriteOpacityTableHDF
  PUBLIC ReadOpacityTableHDF

CONTAINS

  SUBROUTINE WriteOpacityTableHDF &
    ( OpacityTable, FileName, WriteOpacity_EmAb_Option, &
      WriteOpacity_Iso_Option, WriteOpacity_NES_Option, &
      WriteOpacity_NNS_Option, WriteOpacity_Pair_Option, &
      WriteOpacity_Brem_Option )
 
    TYPE(OpacityTableType), INTENT(inout)        :: OpacityTable
    CHARACTER(len=*),       INTENT(in)           :: FileName
    LOGICAL,                INTENT(in), OPTIONAL :: WriteOpacity_EmAb_Option
    LOGICAL,                INTENT(in), OPTIONAL :: WriteOpacity_Iso_Option
    LOGICAL,                INTENT(in), OPTIONAL :: WriteOpacity_NES_Option
    LOGICAL,                INTENT(in), OPTIONAL :: WriteOpacity_NNS_Option
    LOGICAL,                INTENT(in), OPTIONAL :: WriteOpacity_Pair_Option
    LOGICAL,                INTENT(in), OPTIONAL :: WriteOpacity_Brem_Option

    LOGICAL           :: WriteOpacity_EmAb
    LOGICAL           :: WriteOpacity_Iso 
    LOGICAL           :: WriteOpacity_NES 
    LOGICAL           :: WriteOpacity_NNS 
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

    IF( PRESENT( WriteOpacity_NNS_Option ) )THEN
      WriteOpacity_NNS = WriteOpacity_NNS_Option
    ELSE
      WriteOpacity_NNS = .FALSE.
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
  
    CALL OpenGroupHDF( "ThermoState", .true., file_id, group_id )
    CALL WriteThermoStateHDF( OpacityTable % TS, group_id )
    CALL CloseGroupHDF( group_id )

    IF( WriteOpacity_EmAb )THEN

    WRITE(*,*) "Writing out EmAb"
      IF( .NOT. ALLOCATED( OpacityTable % EmAb % Names ) )THEN

        ! --- Insert Appropriate Reaction ---
        WRITE(*,'(A4,A)') &
          '', 'OpacityTable % EmAb not allocated.  Write Skipped.'

      ELSE

        CALL OpenGroupHDF &
               ( "EmAb Parameters", .true., file_id, group_id )

        CALL WriteOpacityTableHDF_EmAb_parameters( OpacityTable % EmAb, group_id )
        CALL CloseGroupHDF( group_id )

        CALL OpenGroupHDF &
               ( "EmAb", .true., file_id, group_id )

        CALL WriteVersionAttribute(group_id)

        CALL WriteOpacityTableHDF_EmAb( OpacityTable % EmAb, group_id )
        CALL CloseGroupHDF( group_id )

      END IF

      IF(OpacityTable % EmAb % nuclei_EC_table .gt. 0) THEN

        WRITE(*,*) "Writing EC table spectrum and rate."
        IF( .NOT. ALLOCATED( OpacityTable % EmAb % EC_table_Names ) )THEN

          ! --- Insert Appropriate Reaction ---
          WRITE(*,'(A4,A)') &
            '', 'OpacityTable % EmAb on nuclei using EC table not allocated.  Write Skipped.'

        ELSE

          CALL OpenGroupHDF &
                 ( "EC_table", .true., file_id, group_id )

          CALL WriteOpacityTableHDF_EC_table( OpacityTable % EmAb, group_id )
          CALL CloseGroupHDF( group_id )

        END IF

      ENDIF

    END IF

    IF( WriteOpacity_Iso ) THEN

    WRITE(*,*) "Writing out Iso"
       IF( .NOT. ALLOCATED( OpacityTable % Scat_Iso % Names ) )THEN

        ! --- Insert Appropriate Reaction ---
        WRITE(*,'(A4,A)') &
          '', 'OpacityTable % Scat_Iso not allocated.  Write Skipped.'

      ELSE

        CALL OpenGroupHDF &
               ( "Scat_Iso_Kernels", .true., file_id, group_id )
        CALL WriteOpacityTableHDF_Scat( OpacityTable % Scat_Iso, group_id )

          BLOCK

            CHARACTER(LEN=100), DIMENSION(3) :: tempString

            tempString(1) = &
            "Opacity from isoenergetic scattering, Bruenn and Mezzacappa (1997), Horowitz (1997)"
            tempString(2) = &
            "https://ui.adsabs.harvard.edu/link_gateway/1997PhRvD..56.7529B/doi:10.1103/PhysRevD.56.7529"
            tempString(3) = &
            "https://ui.adsabs.harvard.edu/link_gateway/1997PhRvD..55.4577H/doi:10.1103/PhysRevD.55.4577"

            CALL WriteGroupAttributeHDF_string("Opacity description", tempString, group_id) 

          END BLOCK
    
          BLOCK

            CHARACTER(LEN=100), DIMENSION(3) :: tempString

            tempString(1) = "Weak magnetism corrections for isoenergetic scattering, Horowitz (2002)"
            tempString(2) = "https://ui.adsabs.harvard.edu/link_gateway/2002PhRvD..65d3001H/doi:10.1103/PhysRevD.65.043001"
            IF(OpacityTable % Scat_Iso % weak_magnetism_corrections .gt. 0) THEN
              tempString(3) = "Included."
            ELSE
              tempString(3) = "Not included."
            ENDIF

            CALL WriteGroupAttributeHDF_string("weak_magnetism_corrections", tempString, group_id) 

          END BLOCK

          BLOCK

            CHARACTER(LEN=100), DIMENSION(6) :: tempString

            tempString(1) = "Ion-ion correlations corrections for isoenergetic scattering"
            tempString(2) = "Itoh (1997), Horowitz (1997), Bruenn and Mezzacappa (1997)"
            tempString(3) = &
            "https://ui.adsabs.harvard.edu/link_gateway/1975PThPh..54.1580I/doi:10.1143/PTP.54.1580"
            tempString(4) = &
            "https://ui.adsabs.harvard.edu/link_gateway/1997PhRvD..56.7529B/doi:10.1103/PhysRevD.56.7529"
            tempString(5) = &
            "https://ui.adsabs.harvard.edu/link_gateway/1997PhRvD..55.4577H/doi:10.1103/PhysRevD.55.4577"
            IF(OpacityTable % Scat_Iso % ion_ion_corrections .gt. 0) THEN
              tempString(6) = "Included."
            ELSE
              tempString(6) = "Not included."
            ENDIF

            CALL WriteGroupAttributeHDF_string("ion_ion_corrections", tempString, group_id) 

          END BLOCK

          BLOCK

            CHARACTER(LEN=100), DIMENSION(3) :: tempString

            tempString(1) = "Many-body effects corrections for isoenergetic scattering, Horowitz et al (2017)"
            tempString(2) = "https://ui.adsabs.harvard.edu/link_gateway/2017PhRvC..95b5801H/doi:10.1103/PhysRevC.95.025801"
            IF(OpacityTable % Scat_Iso % many_body_corrections .gt. 0) THEN
              tempString(3) = "Included."
            ELSE
              tempString(3) = "Not included."
            ENDIF

            CALL WriteGroupAttributeHDF_string("many_body_corrections", tempString, group_id) 

          END BLOCK

          BLOCK

            CHARACTER(LEN=100), DIMENSION(5) :: tempString

            WRITE(tempString(5),'(A13,ES11.3E3)') "ga_strange = ", OpacityTable % Scat_Iso % ga_strange

            tempString(1) = "Strange quark controbutions to neutral-current neutrino nucleon interactions"
            tempString(2) = "Suggested value is ga_strange = -0.1 or -0.04 from Hobbs et al 2016"
            tempString(3) = "https://ui.adsabs.harvard.edu/link_gateway/2016PhRvC..93e2801H/doi:10.1103/PhysRevC.93.052801"
            tempString(4) = "Not active if ga_strange = 0"
            !IF(OpacityTable % Scat_Iso % ga_strange .ne. 0.0d0) THEN
            !  tempString(4) = "Included."
            !ELSE
            !  tempString(4) = "Not included."
            !ENDIF

            CALL WriteGroupAttributeHDF_string("strange_quark_contributions", tempString, group_id) 

          END BLOCK
        
          CALL WriteVersionAttribute(group_id)

        CALL CloseGroupHDF( group_id )

      END IF

    END IF

    IF( WriteOpacity_NES .or. WriteOpacity_NNS .or. WriteOpacity_Pair) THEN

      CALL OpenGroupHDF( "EtaGrid", .true., file_id, group_id )
      CALL WriteGridHDF( OpacityTable % EtaGrid, group_id )
      CALL CloseGroupHDF( group_id )

      IF( WriteOpacity_NES ) THEN
      WRITE(*,*) "Writing out NES"

        IF( .NOT. ALLOCATED( OpacityTable % Scat_NES % Names ) )THEN
  
          ! --- Insert Appropriate Reaction ---
          WRITE(*,'(A4,A)') &
            '', 'OpacityTable % Scat_NES not allocated.  Write Skipped.'
  
        ELSE
  
          CALL OpenGroupHDF &
                 ( "Scat_NES_Kernels", .true., file_id, group_id )
          CALL WriteOpacityTableHDF_Scat( OpacityTable % Scat_NES, group_id )

          BLOCK

            CHARACTER(LEN=100), DIMENSION(3) :: tempString

            tempString(1) = "Opacity from NES, Bruenn (1985), Mezzacappa and Bruenn (1993)"
            tempString(2) = "https://ui.adsabs.harvard.edu/link_gateway/1985ApJS...58..771B/doi:10.1086/191056"
            tempString(3) = "https://ui.adsabs.harvard.edu/link_gateway/1993ApJ...410..740M/doi:10.1086/172791"

            CALL WriteGroupAttributeHDF_string("Opacity description", tempString, group_id) 

          END BLOCK

          BLOCK

            CHARACTER(LEN=100), DIMENSION(4) :: tempString

            tempString(1) = "Opacity from NPS, Bruenn (1985), Mezzacappa and Bruenn (1993)"
            tempString(2) = "https://ui.adsabs.harvard.edu/link_gateway/1985ApJS...58..771B/doi:10.1086/191056"
            tempString(3) = "https://ui.adsabs.harvard.edu/link_gateway/1993ApJ...410..740M/doi:10.1086/172791"
            IF(OpacityTable % Scat_NES % NPS .gt. 0) THEN
              tempString(4) = "Included."
            ELSE
              tempString(4) = "Not included."
            ENDIF

            CALL WriteGroupAttributeHDF_string("Neutrino positron scattering", tempString, group_id) 

          END BLOCK

          CALL WriteVersionAttribute(group_id)

          CALL CloseGroupHDF( group_id )
  
        END IF

      END IF

      IF( WriteOpacity_NNS ) THEN
      WRITE(*,*) "Writing out NNS"

        IF( .NOT. ALLOCATED( OpacityTable % Scat_NNS % Names ) )THEN
  
          ! --- Insert Appropriate Reaction ---
          WRITE(*,'(A4,A)') &
            '', 'OpacityTable % Scat_NNS not allocated.  Write Skipped.'
  
        ELSE
  
          CALL OpenGroupHDF &
                 ( "Scat_NNS_Kernels", .true., file_id, group_id )
          CALL WriteOpacityTableHDF_Scat( OpacityTable % Scat_NNS, group_id )

          BLOCK

            CHARACTER(LEN=100), DIMENSION(2) :: tempString

            tempString(1) = "Opacity from NNS, Bruenn et al. (2020)"
            tempString(2) = "https://ui.adsabs.harvard.edu/link_gateway/2020ApJS..248...11B/doi:10.3847/1538-4365/ab7aff"

            CALL WriteGroupAttributeHDF_string("Opacity description", tempString, group_id) 

          END BLOCK

          BLOCK

            CHARACTER(LEN=100), DIMENSION(3) :: tempString

            tempString(1) = "Weak magnetism corrections for NNS, Bruenn et al. (2020)"
            tempString(2) = "https://ui.adsabs.harvard.edu/link_gateway/2020ApJS..248...11B/doi:10.3847/1538-4365/ab7aff"
            IF(OpacityTable % Scat_NNS % weak_magnetism_corrections .gt. 0) THEN
              tempString(3) = "Included."
            ELSE
              tempString(3) = "Not included."
            ENDIF

            CALL WriteGroupAttributeHDF_string("weak_magnetism_corrections", tempString, group_id) 

          END BLOCK

          BLOCK

            CHARACTER(LEN=100), DIMENSION(3) :: tempString

            tempString(1) = "Many-body effects corrections for isoenergetic scattering, Horowitz et al (2017)"
            tempString(2) = "https://ui.adsabs.harvard.edu/link_gateway/2017PhRvC..95b5801H/doi:10.1103/PhysRevC.95.025801"
            IF(OpacityTable % Scat_NNS % many_body_corrections .gt. 0) THEN
              tempString(3) = "Included."
            ELSE
              tempString(3) = "Not included."
            ENDIF

            CALL WriteGroupAttributeHDF_string("many_body_corrections", tempString, group_id) 

          END BLOCK
        
          CALL WriteVersionAttribute(group_id)

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

          BLOCK

            CHARACTER(LEN=100), DIMENSION(3) :: tempString

            tempString(1) = "Opacity from pair production, Bruenn (1985), Mezzacappa and Bruenn (1993)"
            tempString(2) = "https://ui.adsabs.harvard.edu/link_gateway/1985ApJS...58..771B/doi:10.1086/191056"
            tempString(3) = "https://ui.adsabs.harvard.edu/link_gateway/1993ApJ...410..740M/doi:10.1086/172791"

            CALL WriteGroupAttributeHDF_string("Opacity description", tempString, group_id) 

          END BLOCK
 
          CALL WriteVersionAttribute(group_id)

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

        BLOCK

          CHARACTER(LEN=100), DIMENSION(3) :: tempString

          tempString(1) = "Nucleon-nucleon Bremsstrahlung, Hannestad and Raffelt (1998)"
          tempString(2) = "The full spin-density autocorrelation function S_sigma(eps+eps') in units [1/MeV]"
          tempString(3) = "https://ui.adsabs.harvard.edu/link_gateway/1998ApJ...507..339H/doi:10.1086/306303"

          CALL WriteGroupAttributeHDF_string("Opacity description", tempString, group_id) 

        END BLOCK

        CALL WriteVersionAttribute(group_id)

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

  SUBROUTINE WriteOpacityTableHDF_EmAb_parameters( EmAb, group_id )

    TYPE(OpacityTypeEmAb), INTENT(in) :: EmAb
    INTEGER(HID_T),        INTENT(in) :: group_id

    INTEGER(HSIZE_T) :: datasize1d(1)
    INTEGER(HSIZE_T) :: datasize4d(4)
    INTEGER :: ii

    INTEGER, DIMENSION(1)             :: tempInteger

    datasize1d = 1

    tempInteger(1) = EmAb % np_FK
    CALL WriteHDF( "np_FK", tempInteger, group_id, datasize1d )
    tempInteger(1) = EmAb % np_FK_inv_n_decay
    CALL WriteHDF( "np_FK_inv_n_decay", tempInteger, group_id, datasize1d )
    tempInteger(1) = EmAb % np_isoenergetic
    CALL WriteHDF( "np_isoenergetic", tempInteger, group_id, datasize1d )
    tempInteger(1) = EmAb % np_non_isoenergetic
    CALL WriteHDF( "np_non_isoenergetic", tempInteger, group_id, datasize1d )
    tempInteger(1) = EmAb % np_weak_magnetism
    CALL WriteHDF( "np_weak_magnetism", tempInteger, group_id, datasize1d )
    tempInteger(1) = EmAb % nuclei_EC_FFN
    CALL WriteHDF( "nuclei_EC_FFN", tempInteger, group_id, datasize1d )
    tempInteger(1) = EmAb % nuclei_EC_table
    CALL WriteHDF( "nuclei_EC_table", tempInteger, group_id, datasize1d )

    BLOCK
      character(len=100), dimension(3) :: tempString

      tempString(1) = "A list of all tabulated opcities with references."
      tempString(2) = "The attributes as well as the parameters in the hdf5"
      tempString(3) = "file indicate which options were used to build the table."

      CALL WriteGroupAttributeHDF_string("EmAb table description", tempString, group_id) 
    END BLOCK

    BLOCK

      character(len=100), dimension(4) :: tempString

      tempString(1) = "Isoenergetic approximation, Bruenn (1985), Mezzacappa and Bruenn (1993)"
      tempString(2) = "https://ui.adsabs.harvard.edu/link_gateway/1985ApJS...58..771B/doi:10.1086/191056"
      tempString(3) = "https://ui.adsabs.harvard.edu/link_gateway/1993ApJ...410..740M/doi:10.1086/172791"
      IF(EmAb % np_isoenergetic .gt. 0) THEN
        tempString(4) = "Included."
      ELSE
        tempString(4) = "Not included."
      ENDIF

      CALL WriteGroupAttributeHDF_string("np_isoenergetic", tempString, group_id) 

    END BLOCK

    BLOCK

      character(len=100), dimension(5) :: tempString

      tempString(1) = "Recoil, nucleon final-state blocking, and special relativity corrections to EmAb on nucleons"
      tempString(2) = "Reddy et al 1998"
      tempString(3) = "Only used for rho > 1e9, otherwise the isoenergetic approximation is used"
      tempString(4) = "https://ui.adsabs.harvard.edu/link_gateway/1998PhRvD..58a3009R/doi:10.1103/PhysRevD.58.013009"
      IF(EmAb % np_non_isoenergetic .gt. 0) THEN
        tempString(5) = "Included."
      ELSE
        tempString(5) = "Not included."
      ENDIF

      CALL WriteGroupAttributeHDF_string("np_non_isoenergetic", tempString, group_id) 

    END BLOCK

    BLOCK

      character(len=100), dimension(3) :: tempString

      tempString(1) = "Weak magnetism corrections for EmAb on free nucleons, Horowitz (1997)"
      tempString(2) = "https://ui.adsabs.harvard.edu/link_gateway/1997PhRvD..55.4577H/doi:10.1103/PhysRevD.55.4577"
      IF(EmAb % np_weak_magnetism .gt. 0) THEN
        tempString(3) = "Included."
      ELSE
        tempString(3) = "Not included."
      ENDIF

      CALL WriteGroupAttributeHDF_string("np_weak_magnetism", tempString, group_id) 

    END BLOCK

    BLOCK

      character(len=100), dimension(3) :: tempString

      tempString(1) = "Full kinematics rates for EmAb on free nucleons, Fischer et al (2020)"
      tempString(2) = "https://ui.adsabs.harvard.edu/link_gateway/2020PhRvC.101b5804F/doi:10.1103/PhysRevC.101.025804"
      IF(EmAb % np_FK .gt. 0) THEN
        tempString(3) = "Included."
      ELSE
        tempString(3) = "Not included."
      ENDIF

      CALL WriteGroupAttributeHDF_string("np_FK", tempString, group_id) 

    END BLOCK

    BLOCK

      character(len=100), dimension(3) :: tempString

      tempString(1) = "Full kinematics rates for inverse neutron decay, Fischer et al (2020)"
      tempString(2) = "https://ui.adsabs.harvard.edu/link_gateway/2020PhRvC.101b5804F/doi:10.1103/PhysRevC.101.025804"
      IF(EmAb % np_FK_inv_n_decay .gt. 0) THEN
        tempString(3) = "Included."
      ELSE
        tempString(3) = "Not included."
      ENDIF

      CALL WriteGroupAttributeHDF_string("np_FK_inv_n_decay", &
            tempString, group_id) 

    END BLOCK

    BLOCK

      character(len=100), dimension(4) :: tempString

      tempString(1) = "Electron capture on nuclei from Bruenn 1985 using the FFN formalism, Fuller et al (1982)."
      tempString(2) = "https://ui.adsabs.harvard.edu/link_gateway/1982ApJ...252..715F/doi:10.1086/159597"
      tempString(3) = "https://ui.adsabs.harvard.edu/link_gateway/1985ApJS...58..771B/doi:10.1086/191056"
      IF(EmAb % nuclei_EC_FFN .gt. 0) THEN
        tempString(4) = "Included."
      ELSE
        tempString(4) = "Not included."
      ENDIF

      CALL WriteGroupAttributeHDF_string("nuclei_EC_FFN", tempString, group_id) 

    END BLOCK

    BLOCK

      character(len=100), dimension(11) :: tempString

      tempString(1)  = "NSE-folded tabular data, Langanke et al. (2003), Hix et al. (2003)"
      tempString(2)  = "https://ui.adsabs.harvard.edu/link_gateway/2003PhRvL..90x1102L/doi:10.1103/PhysRevLett.90.241102"
      tempString(3)  = "https://ui.adsabs.harvard.edu/link_gateway/2003PhRvL..91t1102H/doi:10.1103/PhysRevLett.91.201102"
      tempString(4)  = "The table stores the rate of electron capture on nuclei"
      tempString(5)  = "and the spectrum of emissivity of electron-neutrinos."
      tempString(6)  = "The rate and spectrum are tabulated independently so that users can interpolated the spectrum"
      tempString(7)  = "and afterwards normalise it to 1."
      tempString(8)  = "Tabulated is jec * e_nu^2 on an equidistant grid from 0..100 MeV with a resolution of dE=0.5 MeV."
      tempString(9)  = "The tabulated spectrum is normalised so that it integrates to 1."
      tempstring(10) = "The spectrum has not been multiplied by the heavy nucleus abundance."
      IF(EmAb % nuclei_EC_table .gt. 0) THEN
        tempString(11) = "Included."
      ELSE
        tempString(11) = "Not included."
      ENDIF
        

      CALL WriteGroupAttributeHDF_string("nuclei_EC_table", tempString, group_id) 

    END BLOCK
  
  END SUBROUTINE WriteOpacityTableHDF_EmAb_parameters

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

  SUBROUTINE WriteOpacityTableHDF_EC_table( EmAb, group_id )

    TYPE(OpacityTypeEmAb), INTENT(in) :: EmAb
    INTEGER(HID_T),        INTENT(in) :: group_id

    INTEGER(HSIZE_T) :: datasize1d(1)
    INTEGER(HSIZE_T) :: datasize3d(3)
    INTEGER(HSIZE_T) :: datasize4d(4)
    INTEGER :: ii

    INTEGER, DIMENSION(1)             :: tempInteger
    REAL(dp)                          :: tmp_real(1)

    CHARACTER(LEN=100), DIMENSION(5) :: tempString

    datasize1d = 1

    datasize1d = EmAb % EC_table_nOpacities !1 !EmAb % nOpacities
    CALL WriteHDF &
           ( "Units", EmAb % EC_table_Units, group_id, datasize1d ) 

    CALL WriteHDF &
           ( "spec_Offsets", EmAb % EC_table_spec_Offsets, group_id, datasize1d )
    CALL WriteHDF &
           ( "rate_Offsets", EmAb % EC_table_rate_Offsets, group_id, datasize1d )

    datasize1d = 1

    tempInteger(1) = EmAb % EC_Table_nE
    CALL WriteHDF &
           ( "nPointsE",   tempInteger, group_id, datasize1d )

    tempInteger(1) = EmAb % EC_Table_nRho
    CALL WriteHDF &
           ( "nPointsRho", tempInteger, group_id, datasize1d )

    tempInteger(1) = EmAb % EC_Table_nT
    CALL WriteHDF &
           ( "nPointsT",   tempInteger, group_id, datasize1d )

    tempInteger(1) = EmAb % EC_Table_nYe
    CALL WriteHDF &
           ( "nPointsYe",  tempInteger, group_id, datasize1d )

    datasize1d = EmAb % EC_Table_nE
    CALL WriteHDF &
           ( "nu_E", &
             EmAb % EC_table_E(:), group_id, datasize1d )

    datasize1d = EmAb % EC_Table_nRho
    CALL WriteHDF &
           ( "rho", &
             EmAb % EC_table_rho(:), group_id, datasize1d )

    datasize1d = EmAb % EC_Table_nT
    CALL WriteHDF &
           ( "T", &
             EmAb % EC_table_T(:), group_id, datasize1d )

    datasize1d = EmAb % EC_Table_nYe
    CALL WriteHDF &
           ( "Ye", &
             EmAb % EC_table_Ye(:), group_id, datasize1d )

    datasize1d = 1

    tmp_real(1) = EmAb % EC_table_rho_min
    CALL WriteHDF &
           ( "minRho",  tmp_real, group_id, datasize1d )
    tmp_real(1) = EmAb % EC_table_rho_max
    CALL WriteHDF &
           ( "maxRho",  tmp_real, group_id, datasize1d )

    tmp_real(1) = EmAb % EC_table_T_min
    CALL WriteHDF &
           ( "minT",  tmp_real, group_id, datasize1d )
    tmp_real(1) = EmAb % EC_table_T_max
    CALL WriteHDF &
           ( "maxT",  tmp_real, group_id, datasize1d )

    tmp_real(1) = EmAb % EC_table_Ye_min
    CALL WriteHDF &
           ( "minYe",  tmp_real, group_id, datasize1d )
    tmp_real(1) = EmAb % EC_table_Ye_max
    CALL WriteHDF &
           ( "maxYe",  tmp_real, group_id, datasize1d )

    datasize4d = [EmAb % EC_Table_nRho, EmAb % EC_Table_nT, EmAb % EC_Table_nYe, EmAb % EC_Table_nE]

    CALL WriteHDF &
           ( "Spectrum", &
             EmAb % EC_table_spec(1) % Values(:,:,:,:), group_id, datasize4d )

    datasize3d = [EmAb % EC_Table_nRho, EmAb % EC_Table_nT, EmAb % EC_Table_nYe]

    CALL WriteHDF &
           ( "Rate", &
             EmAb % EC_table_rate(1) % Values(:,:,:), group_id, datasize3d )

  END SUBROUTINE WriteOpacityTableHDF_EC_table


  SUBROUTINE WriteOpacityTableHDF_Scat( Scat , group_id )

    CLASS(OpacityTypeScat), INTENT(in)    :: Scat
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

    SELECT TYPE ( Scat )

      TYPE IS ( OpacityTypeScatIso )

        datasize1d = 1
        tempInteger(1) = Scat % weak_magnetism_corrections
        CALL WriteHDF( "weak_magnetism_corr", tempInteger, group_id, datasize1d )
        tempInteger(1) = Scat % ion_ion_corrections
        CALL WriteHDF( "ion_ion_corr", tempInteger, group_id, datasize1d )

        tempInteger(1) = Scat % many_body_corrections
        CALL WriteHDF( "many_body_corr", tempInteger, group_id, datasize1d )

        tempReal(1)    = Scat % ga_strange
        CALL WriteHDF( "ga_strange", tempReal, group_id, datasize1d )

        tempInteger(1) = Scat % includes_nucleons
        CALL WriteHDF( "includes_nucleons", tempInteger, group_id, datasize1d )

      TYPE IS ( OpacityTypeScatNES )

        datasize1d = 1
        tempInteger(1) = Scat % NPS
        CALL WriteHDF( "NPS", tempInteger, group_id, datasize1d )

      TYPE IS ( OpacityTypeScatNNS )

        datasize1d = 1
        tempInteger(1) = Scat % weak_magnetism_corrections
        CALL WriteHDF( "weak_magnetism_corr", tempInteger, group_id, datasize1d )

        tempInteger(1) = Scat % many_body_corrections
        CALL WriteHDF( "many_body_corr", tempInteger, group_id, datasize1d )

     END SELECT

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
      FileName_NES_Option, FileName_NNS_Option, FileName_Pair_Option, &
      FileName_Brem_Option, EquationOfStateTableName_Option, Verbose_Option )

    USE MPI
 
    TYPE(OpacityTableType), INTENT(inout)          :: OpacityTable
    CHARACTER(len=*),       INTENT(in),   OPTIONAL :: FileName_EmAb_Option
    CHARACTER(len=*),       INTENT(in),   OPTIONAL :: FileName_Iso_Option
    CHARACTER(len=*),       INTENT(in),   OPTIONAL :: FileName_NES_Option
    CHARACTER(len=*),       INTENT(in),   OPTIONAL :: FileName_NNS_Option
    CHARACTER(len=*),       INTENT(in),   OPTIONAL :: FileName_Pair_Option
    CHARACTER(len=*),       INTENT(in),   OPTIONAL :: FileName_Brem_Option
    CHARACTER(LEN=*),       INTENT(in),   OPTIONAL :: EquationOfStateTableName_Option
    LOGICAL,                INTENT(in),   OPTIONAL :: Verbose_Option

    INTEGER, PARAMETER :: iEmAb = 1
    INTEGER, PARAMETER :: iIso  = 2
    INTEGER, PARAMETER :: iNES  = 3
    INTEGER, PARAMETER :: iNNS  = 4
    INTEGER, PARAMETER :: iPair = 5
    INTEGER, PARAMETER :: iBrem = 6

    LOGICAL            :: ReadOpacity(6)
    LOGICAL            :: Verbose
    CHARACTER(128)     :: FileName(6)
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
    INTEGER            :: nOpac_NNS
    INTEGER            :: nMom_NNS
    INTEGER            :: nOpac_Pair
    INTEGER            :: nMom_Pair
    INTEGER            :: nOpac_Brem
    INTEGER            :: nMom_Brem
    INTEGER            :: buffer(1)
    REAL(dp)           :: tmp_real(1)
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

    IF( PRESENT( FileName_NNS_Option ) &
        .AND. ( LEN( FileName_NNS_Option ) > 1 ) )THEN
      ReadOpacity(iNNS) = .TRUE.
      FileName   (iNNS) = TRIM( FileName_NNS_Option )
    ELSE
      ReadOpacity(iNNS) = .FALSE.
      nOpac_NNS = 0
      nMom_NNS  = 0
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

    !Get parameters for EmAb in order to check if the spectrum
    !and rate for electron capture on nuclei using the LMSH table
    !is present in the hdf5 table
    IF( ReadOpacity(iEmAb) )THEN

      CALL OpenFileHDF( FileName(iEmAb), .FALSE., file_id )

      !We need to check check before if this group exists to be 
      !compatible with legacy tables, as they do not contain
      !the group "EmAb Parameters".

      BLOCK
        CHARACTER(len=150) :: FileName
        INTEGER(SIZE_T)    :: flength

        CALL h5fget_name_f( file_id, FileName, flength, hdferr )

        CALL h5eset_auto_f( 0, hdferr )
  
        CALL h5gopen_f( file_id, "EmAb Parameters", group_id, hdferr )

        IF ( hdferr .ne. 0 ) THEN

          CALL MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )

          IF(myid == 0) THEN

            WRITE(*,*) 'Group EmAb Parameters not found in ', TRIM( FileName )
            WRITE(*,*) 'This most likely means you are using legacy weaklib tables.'
            WRITE(*,*) 'Neither electron capture on nulcei using a LMSH table'
            WRITE(*,*) 'nor corrections to Bruenn85 opacities for EmAb are included.'

          ENDIF

          OpacityTable % EmAb % np_FK               = -1
          OpacityTable % EmAb % np_FK_inv_n_decay   = -1
          OpacityTable % EmAb % np_isoenergetic     = -1
          OpacityTable % EmAb % np_non_isoenergetic = -1
          OpacityTable % EmAb % np_weak_magnetism   = -1
          OpacityTable % EmAb % nuclei_EC_FFN       = -1
          OpacityTable % EmAb % nuclei_EC_table     = -1

          CALL h5eclear_f( hdferr )
          CALL h5eset_auto_f( 1, hdferr )

        ELSE


          CALL OpenGroupHDF &
                 ( "EmAb Parameters", .FALSE., file_id, group_id )

          datasize1d(1) = 1
          Call ReadHDF ( "np_FK", buffer, group_id, datasize1d)
          OpacityTable % EmAb % np_FK = buffer(1)
          Call ReadHDF ( "np_FK_inv_n_decay", buffer, group_id, datasize1d)
          OpacityTable % EmAb % np_FK_inv_n_decay = buffer(1)
          Call ReadHDF ( "np_isoenergetic", buffer, group_id, datasize1d)
          OpacityTable % EmAb % np_isoenergetic = buffer(1)
          Call ReadHDF ( "np_non_isoenergetic", buffer, group_id, datasize1d)
          OpacityTable % EmAb % np_non_isoenergetic = buffer(1)
          Call ReadHDF ( "np_weak_magnetism", buffer, group_id, datasize1d)
          OpacityTable % EmAb % np_weak_magnetism = buffer(1)
          Call ReadHDF ( "nuclei_EC_FFN", buffer, group_id, datasize1d)
          OpacityTable % EmAb % nuclei_EC_FFN = buffer(1)
          Call ReadHDF ( "nuclei_EC_table", buffer, group_id, datasize1d)
          OpacityTable % EmAb % nuclei_EC_table = buffer(1)

          CALL CloseGroupHDF( group_id )

        ENDIF

      END BLOCK

      IF(OpacityTable % EmAb % nuclei_EC_table .gt. 0) THEN

        CALL OpenGroupHDF &
               ( "EC_table", .FALSE., file_id, group_id )


        CALL ReadHDF( "nPointsE", buffer, group_id, datasize1d )
        OpacityTable % EmAb % EC_table_nE = buffer(1)

        CALL ReadHDF( "nPointsRho", buffer, group_id, datasize1d )
        OpacityTable % EmAb % EC_table_nRho = buffer(1)

        CALL ReadHDF( "nPointsT", buffer, group_id, datasize1d )
        OpacityTable % EmAb % EC_table_nT = buffer(1)

        CALL ReadHDF( "nPointsYe", buffer, group_id, datasize1d )
        OpacityTable % EmAb % EC_table_nYe = buffer(1)

        CALL CloseGroupHDF( group_id )

      ENDIF

      CALL CloseFileHDF( file_id )

    ENDIF

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

      !The EmAb group name changed, so we need to check for legacy 
      !tables and use the respective group name here!
      IF (OpacityTable % EmAb % nuclei_EC_table == -1) THEN

        CALL OpenGroupHDF( "EmAb_CorrectedAbsorption", .FALSE., file_id, group_id )
      ELSE

        CALL OpenGroupHDF( "EmAb", .FALSE., file_id, group_id )

      ENDIF

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

    IF (ReadOpacity(iNNS)) THEN

      buffer(1)  = 1 ! for old opacity table safe

      CALL OpenFileHDF( FileName(iNNS), .FALSE., file_id )

      CALL OpenGroupHDF( "Scat_NNS_Kernels", .FALSE., file_id, group_id )

      CALL ReadHDF( "nOpacities", buffer, group_id, datasize1d )

      nOpac_NNS = buffer(1)

      buffer(1)  = 4 ! for old opacity table safe

      CALL ReadHDF( "nMoments", buffer, group_id, datasize1d )

      nMom_NNS = buffer(1)

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
             nMom_NES, nOpac_NNS, nMom_NNS, nOpac_Pair, nMom_Pair, &
             nOpac_Brem, nMom_Brem, nPointsE, nPointsEta, &
             EquationOfStateTableName_Option = EquationOfStateTableName, &
             OpacityThermoState_Option = TS, &
             Verbose_Option = Verbose )

    IF(OpacityTable % EmAb % nuclei_EC_table .gt. 0) THEN

      ALLOCATE(OpacityTable % EmAb % EC_table_rho(OpacityTable % EmAb % EC_table_nRho) )
      ALLOCATE(OpacityTable % EmAb % EC_table_T  (OpacityTable % EmAb % EC_table_nT  ) )
      ALLOCATE(OpacityTable % EmAb % EC_table_Ye (OpacityTable % EmAb % EC_table_nYe ) )
      ALLOCATE(OpacityTable % EmAb % EC_table_E  (OpacityTable % EmAb % EC_table_nE  ) )

      ALLOCATE( OpacityTable % EmAb % EC_table_spec(1) % Values &
              ( OpacityTable % EmAb % EC_table_nRho,   &
                OpacityTable % EmAb % EC_table_nT, &
                OpacityTable % EmAb % EC_table_nYe,   &
                OpacityTable % EmAb % EC_table_nE) )
      ALLOCATE( OpacityTable % EmAb % EC_table_rate(1) % Values &
              ( OpacityTable % EmAb % EC_table_nRho, &
                OpacityTable % EmAb % EC_table_nT,   &
                OpacityTable % EmAb % EC_table_nYe) )

    ENDIF

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

      !The EmAb group name changed, so we need to check for legacy 
      !tables and use the respective group name here!
      IF (OpacityTable % EmAb % nuclei_EC_table == -1) THEN

        CALL OpenGroupHDF( "EmAb_CorrectedAbsorption", .FALSE., file_id, group_id )
      ELSE

        CALL OpenGroupHDF( "EmAb", .FALSE., file_id, group_id )

      ENDIF

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

      IF(OpacityTable % EmAb % nuclei_EC_table .gt. 0) THEN

        CALL OpenGroupHDF &
               ( "EC_table", .FALSE., file_id, group_id )

        datasize1d(1) = 1 !nOpac_EmAb
        CALL ReadHDF &
               ( "spec_Offsets", OpacityTable % EmAb % EC_table_spec_Offsets, group_id, datasize1d )
        CALL ReadHDF &
               ( "rate_Offsets", OpacityTable % EmAb % EC_table_rate_Offsets, group_id, datasize1d )

        CALL ReadHDF &
               ( "Units",   OpacityTable % EmAb % EC_table_Units,   group_id, datasize1d )

        datasize1d = OpacityTable % EmAb % EC_table_nE
        CALL ReadHDF &
               ( "nu_E", &
                 OpacityTable % EmAb % EC_table_E, &
                 group_id, datasize1d )

        datasize1d = OpacityTable % EmAb % EC_table_nRho
        CALL ReadHDF &
               ( "rho", &
                 OpacityTable % EmAb % EC_table_rho, &
                 group_id, datasize1d )

        datasize1d = OpacityTable % EmAb % EC_table_nT
        CALL ReadHDF &
               ( "T", &
                 OpacityTable % EmAb % EC_table_T, &
                 group_id, datasize1d )

        datasize1d = OpacityTable % EmAb % EC_table_nYe
        CALL ReadHDF &
               ( "Ye", &
                 OpacityTable % EmAb % EC_table_Ye, &
                 group_id, datasize1d )

        datasize1d = 1
        CALL ReadHDF( "minRho", tmp_real, group_id, datasize1d )
        OpacityTable % EmAb % EC_table_rho_min = tmp_real(1)

        CALL ReadHDF( "maxRho", tmp_real, group_id, datasize1d )
        OpacityTable % EmAb % EC_table_rho_max = tmp_real(1)

        CALL ReadHDF( "maxT", tmp_real, group_id, datasize1d )
        OpacityTable % EmAb % EC_table_T_max = tmp_real(1)

        CALL ReadHDF( "minT", tmp_real, group_id, datasize1d )
        OpacityTable % EmAb % EC_table_T_min = tmp_real(1)

        CALL ReadHDF( "maxYe", tmp_real, group_id, datasize1d )
        OpacityTable % EmAb % EC_table_Ye_max = tmp_real(1)

        CALL ReadHDF( "minYe", tmp_real, group_id, datasize1d )
        OpacityTable % EmAb % EC_table_Ye_min = tmp_real(1)

        datasize4d = [OpacityTable % EmAb % EC_table_nRho, OpacityTable % EmAb % EC_table_nT, &
                      OpacityTable % EmAb % EC_table_nYe,  OpacityTable % EmAb % EC_table_nE]

        OpacityTable % EmAb % EC_table_Names = "Electron Neutrino"

        CALL ReadHDF &
               ( "Spectrum", &
                 OpacityTable % EmAb % EC_table_spec(1) % Values, &
                 group_id, datasize4d )

        datasize3d = [OpacityTable % EmAb % EC_table_nRho, &
                      OpacityTable % EmAb % EC_table_nT,   &
                      OpacityTable % EmAb % EC_table_nYe]

        CALL ReadHDF &
               ( "Rate", &
                 OpacityTable % EmAb % EC_table_rate(1) % Values, &
                 group_id, datasize3d )

        CALL CloseGroupHDF( group_id )

      ENDIF

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


    IF( ReadOpacity(iNNS) ) THEN

      IF( Verbose )THEN

        WRITE(*,'(A6,A9,A)') '', 'Reading: ', TRIM( FileName(iNNS) )

      END IF

      CALL OpenFileHDF( FileName(iNNS), .FALSE., file_id )

      CALL OpenGroupHDF &
             ( "Scat_NNS_Kernels", .FALSE., file_id, group_id )

      datasize2d(1) = OpacityTable % Scat_NNS % nOpacities
      datasize2d(2) = OpacityTable % Scat_NNS % nMoments

      CALL ReadHDF &
             ( "Offsets", OpacityTable % Scat_NNS % Offsets, &
               group_id, datasize2d )

      CALL ReadHDF &
             ( "Units",   OpacityTable % Scat_NNS % Units,   &
               group_id, datasize2d )

      datasize5d(1:2) = OpacityTable % EnergyGrid % nPoints
      datasize5d(3)   = OpacityTable % Scat_NNS % nMoments
      datasize5d(4)   = &
        OpacityTable % TS % nPoints(OpacityTable % TS % Indices % iT)
      datasize5d(5)   = OpacityTable % EtaGrid % nPoints

      OpacityTable % Scat_NNS % Names(1) = "Kernels";

      CALL ReadHDF &
             ( TRIM( OpacityTable % Scat_NNS % Names(1) ), &
               OpacityTable % Scat_NNS % Kernel(1) % Values, &
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
    Call ReadHDF( "Names",   EmAb % Names,   group_id, datasize1d )
    Call ReadHDF( "Units",   EmAb % Units,   group_id, datasize1d )

    datasize1d(1) = 1  
    CALL ReadHDF( "Offsets", EmAb % EC_table_spec_Offsets, group_id, datasize1d )
    Call ReadHDF( "Names",   EmAb % EC_table_Names,   group_id, datasize1d )
    Call ReadHDF( "Units",   EmAb % EC_table_Units,   group_id, datasize1d )

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

    USE MPI

    CLASS(OpacityTypeScat),INTENT(inout) :: Scat
    INTEGER(HID_T), INTENT(in)           :: group_id

    INTEGER(HSIZE_T), DIMENSION(1)       :: datasize1d
    INTEGER(HSIZE_T), DIMENSION(2)       :: datasize2d
    INTEGER(HSIZE_T), DIMENSION(4)       :: datasize4d
    INTEGER(HSIZE_T), DIMENSION(5)       :: datasize5d
    INTEGER                              :: i
    INTEGER, DIMENSION(1)                :: buffer
    REAL(dp), DIMENSION(1)               :: bufferReal
    INTEGER(HID_T)                       :: subgroup_id

    CHARACTER(len=150)                   :: FileName
    INTEGER(SIZE_T)                      :: flength
    INTEGER(HID_T)                       :: dataset_id

    SELECT TYPE ( Scat )

      TYPE IS ( OpacityTypeScatIso )

        CALL h5fget_name_f( group_id, FileName, flength, hdferr )
          
        CALL h5eset_auto_f( 0, hdferr )
 
        CALL h5dopen_f( group_id, "weak_magnetism_corr", dataset_id, hdferr )

        IF( hdferr .ne. 0 ) THEN
          
          CALL MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )

          IF(myid == 0) THEN
            WRITE(*,*) 'Dataset weak_magnetism_corr not found in ', TRIM( FileName )
            WRITE(*,*) 'This most likely means you are using legacy weaklib tables.'
          ENDIF

          CALL h5eclear_f( hdferr )
          CALL h5eset_auto_f( 1, hdferr )
            
        ELSE
          datasize1d(1) = 1
          CALL ReadHDF( "weak_magnetism_corr", buffer, group_id, datasize1d )
          Scat % weak_magnetism_corrections = buffer(1)

        ENDIF

        CALL h5eset_auto_f( 0, hdferr )
 
        CALL h5dopen_f( group_id, "ion_ion_corr", dataset_id, hdferr )

        IF( hdferr .ne. 0 ) THEN
          CALL MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )

          IF(myid == 0) THEN
            WRITE(*,*) 'Dataset ion_ion_corr not found in ', TRIM( FileName )
            WRITE(*,*) 'This most likely means you are using legacy weaklib tables.'
          ENDIF

          CALL h5eclear_f( hdferr )
          CALL h5eset_auto_f( 1, hdferr )
        ELSE
          datasize1d(1) = 1
          CALL ReadHDF( "ion_ion_corr", buffer, group_id, datasize1d )
          Scat % ion_ion_corrections = buffer(1)
        ENDIF

        CALL h5eset_auto_f( 0, hdferr )
 
        CALL h5dopen_f( group_id, "many_body_corr", dataset_id, hdferr )

        IF( hdferr .ne. 0 ) THEN
          CALL MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )

          IF(myid == 0) THEN
            WRITE(*,*) 'Dataset many_body_corr not found in ', TRIM( FileName )
            WRITE(*,*) 'This most likely means you are using legacy weaklib tables.'
          ENDIF

          CALL h5eclear_f( hdferr )
          CALL h5eset_auto_f( 1, hdferr )
        ELSE
          datasize1d(1) = 1
          CALL ReadHDF( "many_body_corr", buffer, group_id, datasize1d )
          Scat % many_body_corrections = buffer(1)
        ENDIF

        CALL h5eset_auto_f( 0, hdferr )
 
        CALL h5dopen_f( group_id, "ga_strange", dataset_id, hdferr )

        IF( hdferr .ne. 0 ) THEN
          CALL MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )

          IF(myid == 0) THEN
            WRITE(*,*) 'Dataset ga_strange not found in ', TRIM( FileName )
            WRITE(*,*) 'This most likely means you are using legacy weaklib tables.'
          ENDIF

          CALL h5eclear_f( hdferr )
          CALL h5eset_auto_f( 1, hdferr )
        ELSE
          datasize1d(1) = 1
          CALL ReadHDF( "ga_strange", bufferReal, group_id, datasize1d )
          Scat % ga_strange = bufferReal(1)
        ENDIF

        CALL h5eset_auto_f( 0, hdferr )
 
        CALL h5dopen_f( group_id, "includes_nucleons", dataset_id, hdferr )

        IF( hdferr .ne. 0 ) THEN
          CALL MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )

          IF(myid == 0) THEN
            WRITE(*,*) 'Dataset includes_nucleons not found in ', TRIM( FileName )
            WRITE(*,*) 'This most likely means you are using legacy weaklib tables.'
          ENDIF

          CALL h5eclear_f( hdferr )
          CALL h5eset_auto_f( 1, hdferr )
        ELSE
          datasize1d(1) = 1
          CALL ReadHDF( "includes_nucleons", buffer, group_id, datasize1d )
          Scat % includes_nucleons = buffer(1)
        ENDIF

      TYPE IS ( OpacityTypeScatNES )

        CALL h5fget_name_f( group_id, FileName, flength, hdferr )
          
        CALL h5eset_auto_f( 0, hdferr )
 
        CALL h5dopen_f( group_id, "NPS", dataset_id, hdferr )

        IF( hdferr .ne. 0 ) THEN
          CALL MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )

          IF(myid == 0) THEN
            WRITE(*,*) 'Dataset NPS not found in ', TRIM( FileName )
            WRITE(*,*) 'This most likely means you are using legacy weaklib tables.'
          ENDIF

          CALL h5eclear_f( hdferr )
          CALL h5eset_auto_f( 1, hdferr )
            
        ELSE
          datasize1d(1) = 1
          CALL ReadHDF( "NPS", buffer, group_id, datasize1d )
          Scat % NPS = buffer(1)

        ENDIF

      TYPE IS ( OpacityTypeScatNNS )

        CALL h5fget_name_f( group_id, FileName, flength, hdferr )
          
        CALL h5eset_auto_f( 0, hdferr )
 
        CALL h5dopen_f( group_id, "weak_magnetism_corr", dataset_id, hdferr )

        IF( hdferr .ne. 0 ) THEN
          
          CALL MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )

          IF(myid == 0) THEN
            WRITE(*,*) 'Dataset weak_magnetism_corr not found in ', TRIM( FileName )
            WRITE(*,*) 'This most likely means you are using legacy weaklib tables.'
          ENDIF

          CALL h5eclear_f( hdferr )
          CALL h5eset_auto_f( 1, hdferr )
            
        ELSE
          datasize1d(1) = 1
          CALL ReadHDF( "weak_magnetism_corr", buffer, group_id, datasize1d )
          Scat % weak_magnetism_corrections = buffer(1)

        ENDIF

        CALL h5eset_auto_f( 0, hdferr )
 
        CALL h5dopen_f( group_id, "many_body_corr", dataset_id, hdferr )

        IF( hdferr .ne. 0 ) THEN
          CALL MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )

          IF(myid == 0) THEN
            WRITE(*,*) 'Dataset many_body_corr not found in ', TRIM( FileName )
            WRITE(*,*) 'This most likely means you are using legacy weaklib tables.'
          ENDIF

          CALL h5eclear_f( hdferr )
          CALL h5eset_auto_f( 1, hdferr )
        ELSE
          datasize1d(1) = 1
          CALL ReadHDF( "many_body_corr", buffer, group_id, datasize1d )
          Scat % many_body_corrections = buffer(1)
        ENDIF

    END SELECT

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
