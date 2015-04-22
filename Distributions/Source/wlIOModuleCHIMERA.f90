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

    CHARACTER(len=*), INTENT(in) :: FileName
    INTEGER, INTENT(in)          :: MaxZone
    INTEGER, INTENT(in), OPTIONAL :: SkipLinesOption
!    REAL(dp), DIMENSION(:), INTENT(out) :: r, rho, T, Ye, s
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
    !REAL(dp), DIMENSION(MaxZone) :: u, v, w, dr, dummy14, flat, lum, &
    !                                       rstmss, Tmev, ah

    257 FORMAT (1x,i4,6es11.3,a1,es11.3,L3,es11.3,es14.6,4(es11.3),es10.2,a1,es11.3,2x,a1)
    258 FORMAT (1x,i4,9es11.3)
    259 FORMAT (1x,i4,8es11.4)
    260 FORMAT (1x,i4,11es11.3)
    261 FORMAT (1x,i4,14es11.3)

    SkipLines = 0
    IF ( Present ( SkipLinesOption ) ) &
      SkipLines = SkipLinesOption

    OPEN(unit = FileUnit, file = FileName, status ="old", iostat=istate )

    DO i = 1, SkipLines
      READ(FileUnit,*)
    END DO

    DO
      IF ( istate /= 0) EXIT

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

      IF ( TRIM(ADJUSTL(stanza_name))&
           == 'Hydrodynamic data -  Pressure, energy, and derivatives' ) THEN
        PRINT*, TRIM(ADJUSTL(stanza_name))
        DO i = 1, 4 
          READ(FileUnit,*)
        END DO

        DO i = MaxZone, 1, -1
          READ(FileUnit, 258) j, dum1(i), p(i), dum2(i), dum3(i), dum4(i), &
                              dum5(i), dum6(i), dum7(i), dum8(i)
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

      IF ( TRIM(ADJUSTL(stanza_name)) == 'Masses and Mass averaged data' ) THEN
        PRINT*, TRIM(ADJUSTL(stanza_name))
        DO i = 1, 2 
          READ(FileUnit,*)
        END DO

        READ(FileUnit,'(a)',iostat=istate) stanza_name
        PRINT*, TRIM(ADJUSTL(stanza_name(11:13)))
        IF ( TRIM(ADJUSTL(stanza_name(11:13))) == 'rho' ) THEN 
          DO i = 1, 2 
            READ(FileUnit,*)
          END DO
          DO i = MaxZone, 1, -1
            READ(FileUnit, 260) j, dum1(i), dum2(i), xn(i), xp(i), xhe(i), xa(i), &
                                   dum3(i), dum4(i), dum5(i), dum6(i), dum7(i)
          END DO
        END IF
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

    a_heavy = 0.d0
    z_heavy = 0.d0
    be_heavy = 0.d0
 
    END DO

    ! Repeat this loop, for each variable if necessary, as reads with lots of dummy arguments
    ! For each required quantity.  It may be necessary to do multiple skiplines, or rather, to 
    ! just "GO TO LINE "  
    ! (1) REUSE 257 w/ dummy arguments: rho T Ye and S are on line 47
    ! (2) p is on line 12024 DONE
    ! (3) e_internal is on line 3093 DONE
    ! (4) xn          xp         xhe       xa are on line 7636 DONE
    ! (5) vsound is on line 9829 DONE
    ! (6) Chemical potentials n_cm_pot p_cm_pot e_cm_pot are on line 28958
    ! (7) a, n, z, ahoy, zhvy, ahehvy, zhehvy, a_nphehvy, z_nphehvy are on line 26032
    ! Create formats for each read

    ! Collect independent, dependent, and new variables and put in format needed for interpolationtest
    CLOSE(FileUnit)

  END SUBROUTINE ReadChimeraProfile1D

END MODULE wlIOModuleCHIMERA
