MODULE wlIOModuleCHIMERA

  USE wlKindModule, ONLY: dp
  !USE wlThermoStateModule
  !USE wlDependentVariablesModule
  !USE wlEquationOfStateTableModule

  implicit none

  PUBLIC ReadChimeraProfile1D 

CONTAINS

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

    257 FORMAT (1x,i4,6es11.3,a1,es11.3,L3,es11.3,es14.6,4(es11.3),es10.2,a1,es11.3,2x,a1)
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

    OPEN(unit = FileUnit, file = FileName, status ="old", iostat=istate )

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

END MODULE wlIOModuleCHIMERA
