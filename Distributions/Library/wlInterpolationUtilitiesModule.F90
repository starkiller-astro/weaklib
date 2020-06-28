MODULE wlInterpolationUtilitiesModule

  USE wlKindModule, ONLY: dp
  USE wlInterpolationModule
#if defined(WEAKLIB_OACC)
  USE openacc, ONLY: acc_async_sync
#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: TEST_wlInterpolationUtilitiesModule

CONTAINS


  SUBROUTINE TEST_wlInterpolationUtilitiesModule

    PRINT*, 'TEST_wlInterpolationUtilitiesModule'

  END SUBROUTINE TEST_wlInterpolationUtilitiesModule

END MODULE wlInterpolationUtilitiesModule
