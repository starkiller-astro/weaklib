MODULE wlIOModuleCHIMERA

  USE wlKindModule, ONLY: dp
  USE HDF5
  USE wlThermoStateModule
  USE wlDependentVariablesModule
  USE wlEquationOfStateTableModule
  USE wlInterpolationModule
  USE wlExtNumericalModule, ONLY: epsilon, zero
  USE wlIOModuleHDF
  USE wlEOSIOModuleHDF

  implicit none

  PUBLIC ReadChimeraProfile1D 
  PUBLIC ReadCHIMERAHDF

CONTAINS

  SUBROUTINE ReadComposeTableHDF( EOSTable, FileName )
    TYPE(EquationOfStateTableType), INTENT(inout) :: EOSTable
    CHARACTER(len=*), INTENT(in)                  :: FileName

    INTEGER, DIMENSION(3)                         :: nPoints
    INTEGER                                       :: nVariables
    INTEGER(HID_T)                                :: file_id
    INTEGER(HID_T)                                :: group_id

    CALL OpenFileHDF( FileName, .false., file_id )

    CALL OpenGroupHDF( "pointsnb", .false., file_id, group_id )

    CALL ReadDimensionsHDF( nPoints, group_id )
    CALL ReadNumberVariablesHDF( nVariables, group_id )
    CALL CloseGroupHDF( group_id )

    CALL AllocateEquationOfStateTable( EOSTable, nPoints , nVariables )

    CALL ReadThermoStateHDF( EOSTable % TS, file_id )

    CALL ReadDependentVariablesHDF( EOSTable % DV, file_id )

    CALL DescribeEquationOfStateTable( EOSTable )

    CALL CloseFileHDF( file_id )

  END SUBROUTINE ReadComposeTableHDF

  SUBROUTINE ReadChimeraProfile1D( FileName, MaxZone, r, rho, T, Ye, p, s, &
                                     e_internal, xn, xp, xhe, xa, chem_n,  &
                                     chem_p, chem_e, a_heavy, z_heavy,     &
                                     be_heavy, u, rstmss, vsound,          &
                                     SkipLinesOption )

    CHARACTER(len=*), INTENT(in)  :: FileName
    INTEGER, INTENT(in)           :: MaxZone
    INTEGER, INTENT(in), OPTIONAL :: SkipLinesOption
    REAL(dp), DIMENSION(:), INTENT(out) :: r, rho, T, Ye, p, s, e_internal, xn, &
                                             xp, xhe, xa, chem_n, chem_p, chem_e, &
                                             a_heavy, z_heavy, be_heavy, u, &
                                             rstmss, vsound  

    INTEGER :: i
    INTEGER :: j, FileUnit, SkipLines, istate
    LOGICAL, DIMENSION(MaxZone) :: L, shock
    CHARACTER(len=1), DIMENSION(MaxZone) :: nnse, EOS
    CHARACTER(LEN=255) :: stanza_name   ! Name of stanza headers
    REAL(dp), DIMENSION(MaxZone) :: v, w, dr, dummy14, flat, lum, &
                                      Tmev, ah, dum1, dum2, dum3, &
                                      dum4, dum5, dum6, dum7,     &
                                      dum8, dum9, dum10, dum11,   &
                                      dum12, dum13

    257 FORMAT (1x,i4,6es11.3,a1,es11.3,L3,es11.3,es14.6,4(es11.3),es10.2,a1,es11.3,1x,a1)
    258 FORMAT (1x,i4,10es12.4)
    259 FORMAT (1x,i4,8es11.4)
    260 FORMAT (1x,i4,11es11.4)
    261 FORMAT (1x,i4,14es11.4)
    262 FORMAT (1x,i4,10es11.4)
    263 FORMAT (1x,i4,18es11.4)
    264 FORMAT (1x,i4,7es11.4)

    SkipLines = 0
    IF ( Present ( SkipLinesOption ) ) &
      SkipLines = SkipLinesOption

    OPEN( newunit = FileUnit, FILE = FileName, status ="old", iostat=istate )

    DO i = 1, SkipLines
      READ(FileUnit,*)
    END DO

    DO
      IF ( istate /= 0 ) EXIT

      READ(FileUnit,'(a)',iostat=istate) stanza_name

      IF ( TRIM(ADJUSTL(stanza_name)) == 'Model configuration' ) THEN
        PRINT*, TRIM(ADJUSTL(stanza_name))
        DO i = 1, 4 
          READ(FileUnit,*)
        END DO

        DO i = MaxZone, 1, -1
          READ(FileUnit, 257) j, u(i), v(i), w(i), r(i), dr(i), dummy14(i), shock(i), &
                              flat(i), L(i), lum(i), rstmss(i), rho(i), Tmev(i),   &
                              T(i), s(i), ah(i), nnse(i), Ye(i), EOS(i)
        END DO
      END IF

      IF ( TRIM(ADJUSTL(stanza_name)) == 'Energy and shock ram pressure data' ) THEN
        PRINT*, TRIM(ADJUSTL(stanza_name))
        DO i = 1, 4 
          READ(FileUnit,*)
        END DO

        DO i = MaxZone, 1, -1
          READ(FileUnit, 258) j, dum1(i), dum2(i), dum3(i), dum4(i), p(i), &
                              dum5(i), dum6(i), dum7(i), dum8(i), dum9(i) 
        END DO
      END IF

      IF ( TRIM(ADJUSTL(stanza_name)) == 'EVH1 Energy Data' ) THEN
        PRINT*, TRIM(ADJUSTL(stanza_name))
        DO i = 1, 4 
          READ(FileUnit,*)
        END DO

        DO i = MaxZone, 1, -1
          READ(FileUnit, 259) j, e_internal(i), dum1(i), dum2(i), dum3(i), dum4(i), &
                              dum5(i), dum6(i), dum7(i)
        END DO
      END IF

      IF ( TRIM(ADJUSTL(stanza_name)) == 'Hydrodynamic data' ) THEN
        PRINT*, TRIM(ADJUSTL(stanza_name))
        DO i = 1, 4 
          READ(FileUnit,*)
        END DO

        DO i = MaxZone, 1, -1
          READ(FileUnit, 261) j, dum1(i), dum2(i), dum3(i), dum4(i), &
                              dum5(i), vsound(i), dum6(i), dum7(i), dum8(i), &
                              dum9(i), dum10(i), dum11(i), dum12(i), dum13(i)
        END DO
      END IF

      IF ( TRIM(ADJUSTL(stanza_name)) == 'Chemical potential data' ) THEN
        PRINT*, TRIM(ADJUSTL(stanza_name))
        DO i = 1, 4 
          READ(FileUnit,*)
        END DO

        DO i = MaxZone, 1, -1
          READ(FileUnit, 262) j, chem_n(i), chem_p(i), dum1(i), chem_e(i), &
                              dum2(i), dum3(i), dum4(i), dum5(i), dum6(i), &
                              dum7(i)
        END DO
      END IF

      IF ( TRIM(ADJUSTL(stanza_name)) == 'Comp. & Mass Data' ) THEN
        PRINT*, TRIM(ADJUSTL(stanza_name))
        DO i = 1, 8 
          READ(FileUnit,*)
        END DO

        DO i = MaxZone, 1, -1
          READ(FileUnit, 263) j, xn(i), xp(i), xa(i), xhe(i), a_heavy(i), &
                              dum1(i), z_heavy(i), dum2(i), dum3(i), dum4(i), &
                              dum5(i), dum6(i), dum7(i), dum12(i), dum8(i),&
                              dum9(i), dum10(i), dum11(i)
        END DO
      END IF

      IF ( TRIM(ADJUSTL(stanza_name)) == 'Auxilary Nucleus Comp. Data and Burn Data' ) THEN
        PRINT*, TRIM(ADJUSTL(stanza_name))
        DO i = 1, 8 
          READ(FileUnit,*)
        END DO

        DO i = MaxZone, 1, -1
          READ(FileUnit, 262) j, dum1(i), dum2(i), be_heavy(i), dum3(i), dum4(i), &
                              dum5(i), dum6(i), dum7(i), dum8(i), dum9(i)
        END DO
      END IF

    END DO

    CLOSE(FileUnit)

  END SUBROUTINE ReadChimeraProfile1D

  SUBROUTINE ReadCHIMERAHDF( Rho, T, Ye, E_Int, Entropy, NSE, imax, nx, ny, &
                             nz, FileName)

    CHARACTER(len=*), INTENT(in)                :: FileName
    INTEGER, INTENT(out) :: imax, nx, ny, nz
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(out) :: Rho
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(out) :: T
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(out) :: Ye
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(out) :: E_Int
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(out) :: Entropy
    INTEGER, DIMENSION(:,:,:), ALLOCATABLE, INTENT(out) :: NSE
    INTEGER, DIMENSION(2) :: indices

    INTEGER(HSIZE_T), DIMENSION(1)              :: datasize1d
    INTEGER(HSIZE_T), DIMENSION(3)              :: datasize3d
    INTEGER(HID_T)                              :: file_id
    INTEGER(HID_T)                              :: group_id

    CALL OpenFileHDF( FileName, .false., file_id )

    CALL OpenGroupHDF( "mesh", .false., file_id, group_id )

    datasize1d(1) = 2
    CALL read_1d_slab_int('radial_index_bound', indices, group_id, &
           datasize1d)
    imax = indices(2)
    nx = imax + 2

    CALL read_1d_slab_int('theta_index_bound', indices, group_id, &
           datasize1d)
    ny = indices(2)

    CALL read_1d_slab_int('phi_index_bound', indices, group_id, &
           datasize1d)
    nz = indices(2)

    CALL CloseGroupHDF( group_id )

    ALLOCATE( Rho( nx, ny, nz ), T( nx, ny, nz ), Ye( nx, ny, nz ),           &
              E_Int( nx, ny, nz ), Entropy( nx, ny, nz ), NSE( nx + 1, ny, nz ) )

    CALL OpenGroupHDF( "fluid", .false., file_id, group_id )

    datasize3d = (/nx,ny,nz/)
    CALL ReadHDF( "rho_c", Rho(:,:,:), &
                              group_id, datasize3d )

    CALL ReadHDF( "t_c", T(:,:,:), &
                              group_id, datasize3d )

    CALL ReadHDF( "ye_c", Ye(:,:,:), &
                              group_id, datasize3d )

    CALL ReadHDF( "e_int", E_Int(:,:,:), &
                              group_id, datasize3d )

    CALL ReadHDF( "entropy", Entropy(:,:,:), &
                              group_id, datasize3d )

    CALL CloseGroupHDF( group_id )

    CALL OpenGroupHDF( "abundance", .false., file_id, group_id )

    datasize3d = SHAPE(NSE)
    CALL ReadHDF( "nse_c", NSE(:,:,:), &
                              group_id, datasize3d )

    CALL CloseGroupHDF( group_id )

  END SUBROUTINE ReadCHIMERAHDF

END MODULE wlIOModuleCHIMERA
