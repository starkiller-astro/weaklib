MODULE wlExtNumericalModule

!  replicates "numerical_module" functions for routines/modules extracted from Chimera

  USE wlKindModule, ONLY: dp

  REAL(dp), PARAMETER, PUBLIC :: zero    = 0.d0
  REAL(dp), PARAMETER, PUBLIC :: one     = 1.d0
  REAL(dp), PARAMETER, PUBLIC :: half    = 1.d0/2.d0
  REAL(dp), PARAMETER, PUBLIC :: third   = 1.d0/3.d0

  REAL(dp), PARAMETER, PUBLIC :: pi      = 3.1415926535897932385d0 ! pi

  REAL(dp), PARAMETER, PUBLIC :: epsilon = 1.d-100 ! small number

  PRIVATE

END MODULE wlExtNumericalModule
