MODULE wlIOModuleCHIMERA

  USE wlKindModule, ONLY: dp
  !USE wlThermoStateModule
  !USE wlDependentVariablesModule
  !USE wlEquationOfStateTableModule

  implicit none

  PUBLIC ReadChimeraProfile1D 

CONTAINS

  SUBROUTINE ReadChimeraProfile1D( FileName, MaxZone, rho, T, Ye, SkipLinesOption )

    CHARACTER(len=*), INTENT(in) :: FileName
    INTEGER, INTENT(in)          :: MaxZone
    INTEGER, INTENT(in), OPTIONAL :: SkipLinesOption
    REAL(dp), DIMENSION(:), INTENT(out) :: rho, T, Ye

    INTEGER :: i
    INTEGER :: j, FileUnit, SkipLines
    LOGICAL, DIMENSION(MaxZone) :: L, shock
    CHARACTER(len=1), DIMENSION(MaxZone) :: nnse, EOS
    REAL(dp), DIMENSION(MaxZone) :: u, v, w, r, dr, dummy14, flat, lum, &
                                           rstmss, Tmev, s, ah

    257 FORMAT (1x,i4,6es11.3,a1,es11.3,L3,es11.3,es14.6,4(es11.3),es10.2,a1,es11.3,2x,a1)

    SkipLines = 0
    IF ( Present ( SkipLinesOption ) ) &
      SkipLines = SkipLinesOption

    OPEN(unit = FileUnit, file = FileName, status ="old" )

    DO i = 1, SkipLines
      READ(FileUnit,*)
    END DO

    DO i = MaxZone, 1, -1
      READ(FileUnit, 257) j, u(i), v(i), w(i), r(i), dr(i), dummy14(i), shock(i), &
                          flat(i), L(i), lum(i), rstmss(i), rho(i), Tmev(i),   &
                          T(i), s(i), ah(i), nnse(i), Ye(i), EOS(i)

    END DO

    CLOSE(FileUnit)

  END SUBROUTINE ReadChimeraProfile1D

END MODULE wlIOModuleCHIMERA
