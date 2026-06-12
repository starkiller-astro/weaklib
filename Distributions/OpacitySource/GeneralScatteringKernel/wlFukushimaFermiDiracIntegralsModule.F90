MODULE wlFukushimaFermiDiracIntegralsModule

  USE wlKindModule, ONLY: DP

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: F0z_Fukushima, F1z_Fukushima, F2z_Fukushima

CONTAINS

PURE ELEMENTAL FUNCTION F0z_Fukushima(y) RESULT(fd)
    !=============================================================================
    ! Computes the Fermi-Dirac integral of order k=0 (F_0) using a minimax 
    ! rational approximation.
    !
    ! Reference: Fukushima, T. (2014, App. Math. Comp.) 
    !=============================================================================

    ! Assuming you are using your project's kind module. 
    ! If not, use: USE iso_fortran_env, ONLY: DP => real64
    USE wlKindModule, ONLY: DP 
    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: y
    REAL(DP)             :: fd

    REAL(DP) :: x, ex, t, s
    REAL(DP) :: num, den

    ! 7.389056... is exactly exp(2.0)
    REAL(DP), PARAMETER :: EXP_TWO = 7.38905609893065023_DP

    x = -ABS(y)

    IF (x < -2.0_DP) THEN
      ! =======================================================================
      ! Regime 1: Strongly to moderately non-degenerate (x < -2)
      ! =======================================================================
      ex = EXP(x)
      t  = ex * EXP_TWO
      
      ! Evaluate numerator polynomial using Horner's method
      num = 22696.2126132366633_DP          &
          + t * (5222.0667923565138_DP      &
          + t * (357.623326425354522_DP     &
          + t * (6.9167792879948140_DP      &
          + t * 0.00200096064827815813_DP)))
          
      ! Evaluate denominator polynomial using Horner's method
      den = 45392.4252264733267_DP          &
          + t * (14539.5980679273792_DP     &
          + t * (1611.36476693109675_DP     &
          + t * (71.072178562726798_DP      &
          + t)))
          
      fd = ex * (1.0_DP - ex * (num / den))

    ELSE
      ! =======================================================================
      ! Regime 2: Transition / moderately degenerate (-2 <= x <= 0)
      ! =======================================================================
      s = -0.5_DP * x
      t = 1.0_DP - s
      
      ! Evaluate numerator polynomial
      num = 159.601717762460980_DP          &
          + t * (23.7193942338278703_DP     &
          + t * (0.377783268730614356_DP    &
          + t * (10.5181677709577503_DP     &
          + t * (3.78181326142271599_DP     &
          + t * (-0.441998676614933572_DP   &
          + t * (-0.450072959113928254_DP   &
          + t * (-0.0734798777056723512_DP  &
          + t * 0.000915454570009894267_DP)))))))
          
      ! Evaluate denominator polynomial
      den = 284.26032127745967_DP           &
          + s * (315.2592651624449_DP       &
          + s * (310.2713981221035_DP       &
          + s * (206.21640678892182_DP      &
          + s * (96.77898293084927_DP       &
          + s * (35.456591489081173_DP      &
          + s * (8.1762315442738975_DP      &
          + s))))))
          
      fd = num / den
    END IF

    ! =========================================================================
    ! Regime 3: Strongly degenerate mapping (y > 0)
    ! Invokes the analytical inversion: F_0(y) = y + F_0(-y)
    ! =========================================================================
    IF (y > 0.0_DP) THEN
      fd = fd + y
    END IF

  END FUNCTION F0z_Fukushima

  PURE ELEMENTAL FUNCTION F1z_Fukushima(y) RESULT(fd)
    !=============================================================================
    ! Computes the Fermi-Dirac integral of order k=1 using a minimax rational 
    ! approximation.
    !
    ! Reference: Fukushima, T. (2014, App. Math. Comp.) 
    ! Original Author: Toshio Fukushima <Toshio.Fukushima@nao.ac.jp>
    !=============================================================================
    
    REAL(DP), INTENT(IN) :: y
    REAL(DP)             :: fd

    REAL(DP) :: x, ex, t, s
    REAL(DP) :: num, den

    ! --- Mathematical Constants ---
    ! 1.6449340... is exactly pi^2 / 6 (the asymptotic limit of F_1 for large y)
    REAL(DP), PARAMETER :: PI_SQ_OVER_6 = 1.64493406684822644_DP
    ! 7.389056... is exactly exp(2.0)
    REAL(DP), PARAMETER :: EXP_TWO      = 7.38905609893065023_DP

    x = -ABS(y)

    IF (x < -2.0_DP) THEN
      ! =======================================================================
      ! Regime 1: Strongly to moderately non-degenerate (x < -2)
      ! =======================================================================
      ex = EXP(x)
      t  = ex * EXP_TWO
      
      ! Evaluate numerator polynomial using Horner's method
      num = 22189.1070807945062_DP          &
          + t * (4915.92700908746777_DP     &
          + t * (322.901386168881348_DP     &
          + t * (5.9897442965804548_DP      &
          + t * 0.00397641173774375092_DP)))
          
      ! Evaluate denominator polynomial using Horner's method
      den = 88756.428323178025_DP           &
          + t * (25002.3197546553836_DP     &
          + t * (2389.06277237306633_DP     &
          + t * (88.376214553692756_DP      &
          + t)))
          
      fd = ex * (1.0_DP - ex * (num / den))

    ELSE
      ! =======================================================================
      ! Regime 2: Transition / moderately degenerate (-2 <= x <= 0)
      ! =======================================================================
      s = -0.5_DP * x
      t = 1.0_DP - s
      
      ! Evaluate numerator polynomial
      num = 145.488167182330098_DP          &
          + t * (251.392824471576922_DP     &
          + t * (56.6537141912783024_DP     &
          + t * (17.9918985363509694_DP     &
          + t * (20.1369115558099802_DP     &
          + t * (7.09659390228556164_DP     &
          + t * (0.199701180197912643_DP    &
          + t * (-0.403173132925886253_DP   &
          - t * 0.0792966701498222697_DP)))))))
          
      ! Evaluate denominator polynomial
      den = 606.0757707716040_DP            &
          + s * (374.1806357435014_DP       &
          + s * (252.1367051536344_DP       &
          + s * (27.2746245830016_DP        &
          + s * (-61.57766112137513_DP      &
          + s * (-53.72117554363975_DP      &
          + s * (-25.678454878692950_DP     &
          + s * (-7.1995819520154718_DP     &
          - s)))))))
          
      fd = num / den
    END IF

    ! =========================================================================
    ! Regime 3: Strongly degenerate mapping (y > 0)
    ! Invokes the analytical inversion: F_1(y) = y^2/2 + pi^2/6 - F_1(-y)
    ! =========================================================================
    IF (y > 0.0_DP) THEN
      fd = -fd + PI_SQ_OVER_6 + 0.5_DP * y * y
    END IF

  END FUNCTION F1z_Fukushima

  PURE ELEMENTAL FUNCTION F2z_Fukushima(y) RESULT(fd)
    !=============================================================================
    ! Computes the Fermi-Dirac integral of order k=2 (F_2) using a minimax 
    ! rational approximation.
    ! Note: Fukushima denotes this as 'F2z_Fukushima' due to half-integer indexing schemes.
    !
    ! Reference: Fukushima, T. (2014, App. Math. Comp.) 
    !=============================================================================

    REAL(DP), INTENT(IN) :: y
    REAL(DP)             :: fd

    REAL(DP) :: x, ex, t, s
    REAL(DP) :: num, den

    ! --- Mathematical Constants ---
    REAL(DP), PARAMETER :: EXP_TWO      = 7.38905609893065023_DP
    ! 3.289868... is exactly pi^2 / 3
    REAL(DP), PARAMETER :: PI_SQ_OVER_3 = 3.28986813369645287_DP
    REAL(DP), PARAMETER :: ONE_THIRD    = 1.0_DP / 3.0_DP

    x = -ABS(y)

    IF (x < -2.0_DP) THEN
      ! =======================================================================
      ! Regime 1: Strongly to moderately non-degenerate (x < -2)
      ! =======================================================================
      ex = EXP(x)
      t  = ex * EXP_TWO
      
      ! Evaluate numerator polynomial using Horner's method
      num = 1914.06748184935743_DP          &
          + t * (273.085756700981399_DP     &
          + t * (8.5861610217850095_DP      &
          + t * 0.0161890243763741414_DP))
          
      ! Evaluate denominator polynomial using Horner's method
      den = 7656.2699273974454_DP           &
          + t * (1399.35442210906621_DP     &
          + t * (72.929152915475392_DP      &
          + t))
          
      fd = ex * (2.0_DP - ex * (num / den))

    ELSE
      ! =======================================================================
      ! Regime 2: Transition / moderately degenerate (-2 <= x <= 0)
      ! =======================================================================
      s = -0.5_DP * x
      t = 1.0_DP - s
      
      ! Evaluate numerator polynomial
      num = 2711.49678259128843_DP          &
          + t * (1299.85460914884154_DP     &
          + t * (222.606134197895041_DP     &
          + t * (172.881855215582924_DP     &
          + t * (112.951038040682055_DP     &
          + t * (24.0376482128898634_DP     &
          + t * (-2.68393549333878715_DP    &
          + t * (-2.14077421411719935_DP    &
          - t * 0.326188299771397236_DP)))))))
          
      ! Evaluate denominator polynomial
      den = 2517.1726659917047_DP           &
          + s * (3038.7689794575778_DP       &
          + s * (2541.7823512372631_DP       &
          + s * (1428.0589853413436_DP       &
          + s * (531.62378035996132_DP       &
          + s * (122.54595216479181_DP       &
          + s * (8.395768655115050_DP        &
          + s * (-3.9142702096919080_DP      &
          - s)))))))
          
      fd = num / den
    END IF

    ! =========================================================================
    ! Regime 3: Strongly degenerate mapping (y > 0)
    ! Invokes the analytical inversion: F_2(y) = y^3/3 + y*(pi^2/3) + F_2(-y)
    ! =========================================================================
    IF (y > 0.0_DP) THEN
      fd = fd + y * (PI_SQ_OVER_3 + ONE_THIRD * y * y)
    END IF

  END FUNCTION F2z_Fukushima

END MODULE wlFukushimaFermiDiracIntegralsModule