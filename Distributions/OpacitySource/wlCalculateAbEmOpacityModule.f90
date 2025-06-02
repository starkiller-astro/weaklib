MODULE wlCalculateAbEmOpacityModule

    USE wlKindModule, ONLY: dp
    USE physcnst_module, ONLY: Gw, mp, hbar, cvel, pi

  IMPLICIT NONE
  PRIVATE

  REAL(double), PARAMETER     :: g2    = ( Gw/mp**2 )**2 * ( hbar * cvel )**2   ! square of the Fermi constant in units of cm^{2} MeV^{-2}

CONTAINS

  SUBROUTINE ElasticAbsorptionOpacity(D, T, Ye, Ym)



  END SUBROUTINE ElasticAbsorptionOpacity


END MODULE wlCalculateAbEmOpacityModule
