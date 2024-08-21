MODULE wlEOSComponentsInversionModule

  USE wlKindModule, ONLY: &
    dp
  USE wlInterpolationModule, ONLY: &
    LogInterpolateSingleVariable_2D_Custom_Point
  USE wlInterpolationUtilitiesModule, ONLY: &
    Index1D_Lin, &
    Index1D_Log
	
  IMPLICIT NONE
  PRIVATE
	
  PUBLIC :: InitializeEOSComponentsInversion
  PUBLIC :: ComputeTemperatureWith_DEYpYl
  PUBLIC :: ComputeTemperatureWith_DEYpYl_Many
  PUBLIC :: ComputeTemperatureWith_DEYpYl_Single
  PUBLIC :: ComputeTemperatureWith_DEYpYl_Single_Guess
  PUBLIC :: ComputeTemperatureWith_DEYpYl_Single_Guess_Error
  PUBLIC :: ComputeTemperatureWith_DEYpYl_Single_Guess_NoError
  PUBLIC :: ComputeTemperatureWith_DEYpYl_Single_NoGuess
  PUBLIC :: ComputeTemperatureWith_DEYpYl_Single_NoGuess_Error
  PUBLIC :: ComputeTemperatureWith_DEYpYl_Single_NoGuess_NoError
  PUBLIC :: ComputeTemperatureWith_DPYpYl
  PUBLIC :: ComputeTemperatureWith_DPYpYl_Many
  PUBLIC :: ComputeTemperatureWith_DPYpYl_Single_Guess
  PUBLIC :: ComputeTemperatureWith_DPYpYl_Single_Guess_Error
  PUBLIC :: ComputeTemperatureWith_DPYpYl_Single_Guess_NoError
  PUBLIC :: ComputeTemperatureWith_DPYpYl_Single_NoGuess
  PUBLIC :: ComputeTemperatureWith_DPYpYl_Single_NoGuess_Error
  PUBLIC :: ComputeTemperatureWith_DPYpYl_Single_NoGuess_NoError
  PUBLIC :: ComputeTemperatureWith_DSYpYl
  PUBLIC :: ComputeTemperatureWith_DSYpYl_Many
  PUBLIC :: ComputeTemperatureWith_DSYpYl_Single_Guess
  PUBLIC :: ComputeTemperatureWith_DSYpYl_Single_Guess_Error
  PUBLIC :: ComputeTemperatureWith_DSYpYl_Single_Guess_NoError
  PUBLIC :: ComputeTemperatureWith_DSYpYl_Single_NoGuess
  PUBLIC :: ComputeTemperatureWith_DSYpYl_Single_NoGuess_Error
  PUBLIC :: ComputeTemperatureWith_DSYpYl_Single_NoGuess_NoError
  PUBLIC :: DescribeEOSComponentsInversionError
	
  LOGICAL  :: InversionComponentsInitialized
  REAL(dp), PUBLIC :: MinD, MaxD
  REAL(dp), PUBLIC :: MinT, MaxT
  REAL(dp), PUBLIC :: MinYp, MaxYp
  REAL(dp), PUBLIC :: MinE, MaxE
  REAL(dp), PUBLIC :: MinP, MaxP
  REAL(dp), PUBLIC :: MinS, MaxS
	
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET( &
    !$OMP   InversionInitialized, MinD, MaxD, MinT, MaxT, &
    !$OMP   MinY, MaxY, MinE, MaxE, MinP, MaxP, MinS, MaxS )
#elif defined(WEAKLIB_OACC)
    !$ACC DECLARE CREATE( &
    !$ACC   InversionInitialized, MinD, MaxD, MinT, MaxT, &
    !$ACC   MinY, MaxY, MinE, MaxE, MinP, MaxP, MinS, MaxS )
#endif

  INTERFACE ComputeTemperatureWith_DEY
    MODULE PROCEDURE ComputeTemperatureWith_DEYpYl_Many
    MODULE PROCEDURE ComputeTemperatureWith_DEYpYl_Single_Guess_Error
    MODULE PROCEDURE ComputeTemperatureWith_DEYpYl_Single_Guess_NoError
    MODULE PROCEDURE ComputeTemperatureWith_DEYpYl_Single_NoGuess_Error
    MODULE PROCEDURE ComputeTemperatureWith_DEYpYl_Single_NoGuess_NoError
  END INTERFACE ComputeTemperatureWith_DEY

  INTERFACE ComputeTemperatureWith_DEYpYl_Single
    MODULE PROCEDURE ComputeTemperatureWith_DEYpYl_Single_Guess_Error
    MODULE PROCEDURE ComputeTemperatureWith_DEYpYl_Single_Guess_NoError
    MODULE PROCEDURE ComputeTemperatureWith_DEYpYl_Single_NoGuess_Error
    MODULE PROCEDURE ComputeTemperatureWith_DEYpYl_Single_NoGuess_NoError
  END INTERFACE ComputeTemperatureWith_DEYpYl_Single

  INTERFACE ComputeTemperatureWith_DEYpYl_Single_Guess
    MODULE PROCEDURE ComputeTemperatureWith_DEYpYl_Single_Guess_Error
    MODULE PROCEDURE ComputeTemperatureWith_DEYpYl_Single_Guess_NoError
  END INTERFACE ComputeTemperatureWith_DEYpYl_Single_Guess

  INTERFACE ComputeTemperatureWith_DEYpYl_Single_NoGuess
    MODULE PROCEDURE ComputeTemperatureWith_DEYpYl_Single_NoGuess_Error
    MODULE PROCEDURE ComputeTemperatureWith_DEYpYl_Single_NoGuess_NoError
  END INTERFACE ComputeTemperatureWith_DEYpYl_Single_NoGuess

  INTERFACE ComputeTemperatureWith_DPY
    MODULE PROCEDURE ComputeTemperatureWith_DPYpYl_Many
    MODULE PROCEDURE ComputeTemperatureWith_DPYpYl_Single_Guess_Error
    MODULE PROCEDURE ComputeTemperatureWith_DPYpYl_Single_Guess_NoError
    MODULE PROCEDURE ComputeTemperatureWith_DPYpYl_Single_NoGuess_Error
    MODULE PROCEDURE ComputeTemperatureWith_DPYpYl_Single_NoGuess_NoError
  END INTERFACE ComputeTemperatureWith_DPY

  INTERFACE ComputeTemperatureWith_DPYpYl_Single
    MODULE PROCEDURE ComputeTemperatureWith_DPYpYl_Single_Guess_Error
    MODULE PROCEDURE ComputeTemperatureWith_DPYpYl_Single_Guess_NoError
    MODULE PROCEDURE ComputeTemperatureWith_DPYpYl_Single_NoGuess_Error
    MODULE PROCEDURE ComputeTemperatureWith_DPYpYl_Single_NoGuess_NoError
  END INTERFACE ComputeTemperatureWith_DPYpYl_Single

  INTERFACE ComputeTemperatureWith_DPYpYl_Single_Guess
    MODULE PROCEDURE ComputeTemperatureWith_DPYpYl_Single_Guess_Error
    MODULE PROCEDURE ComputeTemperatureWith_DPYpYl_Single_Guess_NoError
  END INTERFACE ComputeTemperatureWith_DPYpYl_Single_Guess

  INTERFACE ComputeTemperatureWith_DPYpYl_Single_NoGuess
    MODULE PROCEDURE ComputeTemperatureWith_DPYpYl_Single_NoGuess_Error
    MODULE PROCEDURE ComputeTemperatureWith_DPYpYl_Single_NoGuess_NoError
  END INTERFACE ComputeTemperatureWith_DPYpYl_Single_NoGuess

  INTERFACE ComputeTemperatureWith_DSY
    MODULE PROCEDURE ComputeTemperatureWith_DSYpYl_Many
    MODULE PROCEDURE ComputeTemperatureWith_DSYpYl_Single_Guess_Error
    MODULE PROCEDURE ComputeTemperatureWith_DSYpYl_Single_Guess_NoError
    MODULE PROCEDURE ComputeTemperatureWith_DSYpYl_Single_NoGuess_Error
    MODULE PROCEDURE ComputeTemperatureWith_DSYpYl_Single_NoGuess_NoError
  END INTERFACE ComputeTemperatureWith_DSY

  INTERFACE ComputeTemperatureWith_DSYpYl_Single
    MODULE PROCEDURE ComputeTemperatureWith_DSYpYl_Single_Guess_Error
    MODULE PROCEDURE ComputeTemperatureWith_DSYpYl_Single_Guess_NoError
    MODULE PROCEDURE ComputeTemperatureWith_DSYpYl_Single_NoGuess_Error
    MODULE PROCEDURE ComputeTemperatureWith_DSYpYl_Single_NoGuess_NoError
  END INTERFACE ComputeTemperatureWith_DSYpYl_Single

  INTERFACE ComputeTemperatureWith_DSYpYl_Single_Guess
    MODULE PROCEDURE ComputeTemperatureWith_DSYpYl_Single_Guess_Error
    MODULE PROCEDURE ComputeTemperatureWith_DSYpYl_Single_Guess_NoError
  END INTERFACE ComputeTemperatureWith_DSYpYl_Single_Guess

  INTERFACE ComputeTemperatureWith_DSYpYl_Single_NoGuess
    MODULE PROCEDURE ComputeTemperatureWith_DSYpYl_Single_NoGuess_Error
    MODULE PROCEDURE ComputeTemperatureWith_DSYpYl_Single_NoGuess_NoError
  END INTERFACE ComputeTemperatureWith_DSYpYl_Single_NoGuess


CONTAINS

  ! COME BACK TO THIS. GOTTA DECIDE HOW TO DO THIS. PROBABLY TO MAKE IT EASIER FOR THE USER 
  ! AND CONSISTENT WITH THE OTHER ROUTINES YOU CALL MUON AND ELECTRON EOS HERE, WITH
  ! MAX(YE)=MAX(YP) AND MIN(YE)=MIN(YP). FOR MUONS IT'S A BIT UNCLEAR TO ME, MAYBE DO
  ! MIN(YMU)=0 AND MAX(YMU)=MAX(YP), BUT THE PROBLEM IS THAT LARGE YMU ARE NOT REALLY 
  ! REALIZED, BUT SUCH A LARGE MUON FRACTION WILL CONTRIBUTE SIGNIFICANTLY TO PRESSURE
  ! AND ENERGY. BUT MAYBE TO DETERIMINE MAXIMUM AND MINIMUM OF P,S,E THIS IS NOT SO IMPORTANT
  SUBROUTINE InitializeEOSComponentsInversion( Ds, Ts, Yps, Es, Ps, Ss, Verbose_Option )

    REAL(dp), INTENT(in) :: Ds(1:)      , Ts(1:)      , Yps(1:)
    REAL(dp), INTENT(in) :: Es(1:,1:,1:), Ps(1:,1:,1:), Ss(1:,1:,1:)
    LOGICAL,  INTENT(in), OPTIONAL :: Verbose_Option

    LOGICAL :: Verbose

    IF( PRESENT( Verbose_Option ) )THEN
      Verbose = Verbose_Option
    ELSE
      Verbose = .FALSE.
    END IF

	! Calculate muons and electrons everywhere
    MinD = MINVAL( Ds ); MaxD = MAXVAL( Ds )
    MinT = MINVAL( Ts ); MaxT = MAXVAL( Ts )
    MinYp = MINVAL( Yps ); MaxYp = MAXVAL( Yps )
    MinE = MINVAL( Es ); MaxE = MAXVAL( Es )
    MinP = MINVAL( Ps ); MaxP = MAXVAL( Ps )
    MinS = MINVAL( Ss ); MaxS = MAXVAL( Ss )

    IF( Verbose )THEN
      WRITE(*,*)
      WRITE(*,'(A4,A)') '', 'InitializeEOSInversion'
      WRITE(*,*)
      WRITE(*,'(A6,A24,2ES11.3E3)') '', 'Min/Max D   [g cm^-3] = ', MinD, MaxD
      WRITE(*,'(A6,A24,2ES11.3E3)') '', 'Min/Max T         [K] = ', MinT, MaxT
      WRITE(*,'(A6,A24,2ES11.3E3)') '', 'Min/Max Y             = ', MinY, MaxY
      WRITE(*,'(A6,A24,2ES11.3E3)') '', 'Min/Max E  [erg g^-1] = ', MinE, MaxE
      WRITE(*,'(A6,A24,2ES11.3E3)') '', 'Min/Max P [dyn cm^-2] = ', MinP, MaxP
      WRITE(*,'(A6,A24,2ES11.3E3)') '', 'Min/Max S       [k_B] = ', MinS, MaxS
      WRITE(*,*)
    END IF

    InversionComponentsInitialized = .TRUE.

#if defined(WEAKLIB_OMP_OL)
    !$OMP TARGET UPDATE TO( &
    !$OMP   InversionInitialized, MinD, MaxD, MinT, MaxT, &
    !$OMP   MinY, MaxY, MinE, MaxE, MinP, MaxP, MinS, MaxS )
#elif defined(WEAKLIB_OACC)
    !$ACC UPDATE DEVICE( &
    !$ACC   InversionInitialized, MinD, MaxD, MinT, MaxT, &
    !$ACC   MinY, MaxY, MinE, MaxE, MinP, MaxP, MinS, MaxS )
#endif

  END SUBROUTINE InitializeEOSComponentsInversion
  
END MODULE wlEOSComponentsInversionModule
