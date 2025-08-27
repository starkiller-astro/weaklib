MODULE wlEOSComponentsSeparateInversionModule

  USE wlKindModule, ONLY: &
    dp
  USE wlEosTemperatureSeparateInversionModule, ONLY: &
    InitializeLeptonTables, &
    InvertTemperatureWith_DEYpYl_Guess, &
    InvertTemperatureWith_DEYpYl_NoGuess, &
    InvertTemperatureWith_DPYpYl_Guess, &
    InvertTemperatureWith_DPYpYl_NoGuess, &
    InvertTemperatureWith_DSYpYl_Guess, &
    InvertTemperatureWith_DSYpYl_NoGuess
  USE wlMuonEOS, ONLY: &
    MuonStateType, FullMuonEOS
  USE wlElectronPhotonEOS, ONLY: &
    ElectronPhotonStateType, ElectronPhotonEOS
  USE wlLeptonEOSModule, ONLY: &
    HelmTableType, MuonTableType
  USE wlHelmMuonIOModuleHDF, ONLY: &
    ReadHelmholtzTableHDF, ReadMuonTableHDF
    
  IMPLICIT NONE
  PRIVATE
  
  PUBLIC :: InitializeEOSComponentsInversion
  PUBLIC :: ComputeTemperatureWith_DEYpYl
  PUBLIC :: ComputeTemperatureWith_DEYpYl_Single
  PUBLIC :: ComputeTemperatureWith_DEYpYl_Single_Guess
  PUBLIC :: ComputeTemperatureWith_DEYpYl_Single_Guess_Error
  PUBLIC :: ComputeTemperatureWith_DEYpYl_Single_Guess_NoError
  PUBLIC :: ComputeTemperatureWith_DEYpYl_Single_NoGuess
  PUBLIC :: ComputeTemperatureWith_DEYpYl_Single_NoGuess_Error
  PUBLIC :: ComputeTemperatureWith_DEYpYl_Single_NoGuess_NoError
  PUBLIC :: ComputeTemperatureWith_DPYpYl
  PUBLIC :: ComputeTemperatureWith_DPYpYl_Single
  PUBLIC :: ComputeTemperatureWith_DPYpYl_Single_Guess
  PUBLIC :: ComputeTemperatureWith_DPYpYl_Single_Guess_Error
  PUBLIC :: ComputeTemperatureWith_DPYpYl_Single_Guess_NoError
  PUBLIC :: ComputeTemperatureWith_DPYpYl_Single_NoGuess
  PUBLIC :: ComputeTemperatureWith_DPYpYl_Single_NoGuess_Error
  PUBLIC :: ComputeTemperatureWith_DPYpYl_Single_NoGuess_NoError
  PUBLIC :: ComputeTemperatureWith_DSYpYl
  PUBLIC :: ComputeTemperatureWith_DSYpYl_Single
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
    !$OMP   InversionComponentsInitialized, MinD, MaxD, MinT, MaxT, &
    !$OMP   MinYp, MaxYp, MinE, MaxE, MinP, MaxP, MinS, MaxS )
#elif defined(WEAKLIB_OACC)
    !$ACC DECLARE CREATE( &
    !$ACC   InversionComponentsInitialized, MinD, MaxD, MinT, MaxT, &
    !$ACC   MinYp, MaxYp, MinE, MaxE, MinP, MaxP, MinS, MaxS )
#endif

  INTERFACE ComputeTemperatureWith_DEYpYl
    MODULE PROCEDURE ComputeTemperatureWith_DEYpYl_Single_Guess_Error
    MODULE PROCEDURE ComputeTemperatureWith_DEYpYl_Single_Guess_NoError
    MODULE PROCEDURE ComputeTemperatureWith_DEYpYl_Single_NoGuess_Error
    MODULE PROCEDURE ComputeTemperatureWith_DEYpYl_Single_NoGuess_NoError
  END INTERFACE ComputeTemperatureWith_DEYpYl

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

  INTERFACE ComputeTemperatureWith_DPYpYl
    MODULE PROCEDURE ComputeTemperatureWith_DPYpYl_Single_Guess_Error
    MODULE PROCEDURE ComputeTemperatureWith_DPYpYl_Single_Guess_NoError
    MODULE PROCEDURE ComputeTemperatureWith_DPYpYl_Single_NoGuess_Error
    MODULE PROCEDURE ComputeTemperatureWith_DPYpYl_Single_NoGuess_NoError
  END INTERFACE ComputeTemperatureWith_DPYpYl

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

  INTERFACE ComputeTemperatureWith_DSYpYl
    MODULE PROCEDURE ComputeTemperatureWith_DSYpYl_Single_Guess_Error
    MODULE PROCEDURE ComputeTemperatureWith_DSYpYl_Single_Guess_NoError
    MODULE PROCEDURE ComputeTemperatureWith_DSYpYl_Single_NoGuess_Error
    MODULE PROCEDURE ComputeTemperatureWith_DSYpYl_Single_NoGuess_NoError
  END INTERFACE ComputeTemperatureWith_DSYpYl

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

  SUBROUTINE InitializeEOSComponentsInversion( Ds, Ts, Yps, Es, Ps, Ss, &
      HelmTable, MuonTable, Verbose_Option )

    REAL(dp), INTENT(in) :: Ds(1:)      , Ts(1:)      , Yps(1:)
    REAL(dp), INTENT(in) :: Es(1:,1:,1:), Ps(1:,1:,1:), Ss(1:,1:,1:)
    TYPE(HelmTableType), POINTER, INTENT(IN) :: HelmTable
    TYPE(MuonTableType), POINTER, INTENT(IN) :: MuonTable
    LOGICAL,  INTENT(in), OPTIONAL :: Verbose_Option
    LOGICAL :: Verbose

    IF( PRESENT( Verbose_Option ) )THEN
      Verbose = Verbose_Option
    ELSE
      Verbose = .FALSE.
    END IF

    ! Initialize Helmholtz and Muon Local Tables
    CALL InitializeLeptonTables( HelmTable, MuonTable )
    
    MinD  = MINVAL( Ds  ); MaxD  = MAXVAL( Ds  )
    MinT  = MINVAL( Ts  ); MaxT  = MAXVAL( Ts  )
    MinYp = MINVAL( Yps ); MaxYp = MAXVAL( Yps )
    MinE  = MINVAL( Es  ); MaxE  = MAXVAL( Es  )
    MinP  = MINVAL( Ps  ); MaxP  = MAXVAL( Ps  )
    MinS  = MINVAL( Ss  ); MaxS  = MAXVAL( Ss  )

    IF( Verbose )THEN
      WRITE(*,*)
      WRITE(*,'(A4,A)') '', 'InitializeEOSInversion Separate Interpolation'
      WRITE(*,*)
      WRITE(*,'(A6,A24,2ES11.3E3)') '', 'Min/Max D   [g cm^-3] = ', MinD , MaxD
      WRITE(*,'(A6,A24,2ES11.3E3)') '', 'Min/Max T         [K] = ', MinT , MaxT
      WRITE(*,'(A6,A24,2ES11.3E3)') '', 'Min/Max Yp            = ', MinYp, MaxYp
      WRITE(*,'(A6,A24,2ES11.3E3)') '', 'Min/Max E  [erg g^-1] = ', MinE , MaxE
      WRITE(*,'(A6,A24,2ES11.3E3)') '', 'Min/Max P [dyn cm^-2] = ', MinP , MaxP
      WRITE(*,'(A6,A24,2ES11.3E3)') '', 'Min/Max S       [k_B] = ', MinS , MaxS
      WRITE(*,*)
    END IF

    InversionComponentsInitialized = .TRUE.

#if defined(WEAKLIB_OMP_OL)
    !$OMP TARGET UPDATE TO( &
    !$OMP   InversionComponentsInitialized, MinD, MaxD, MinT, MaxT, &
    !$OMP   MinYp, MaxYp, MinE, MaxE, MinP, MaxP, MinS, MaxS )
#elif defined(WEAKLIB_OACC)
    !$ACC UPDATE DEVICE( &
    !$ACC   InversionComponentsInitialized, MinD, MaxD, MinT, MaxT, &
    !$ACC   MinYp, MaxYp, MinE, MaxE, MinP, MaxP, MinS, MaxS )
#endif

  END SUBROUTINE InitializeEOSComponentsInversion
  
  INTEGER FUNCTION CheckInputError &
        ( D, X, Ye, Ym, MinX, MaxX )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in) :: D, X, Ye, Ym, MinX, MaxX
    REAL(dp) :: Yp
    
    Yp = Ye + Ym
    
    CheckInputError = 0

    IF ( .NOT. InversionComponentsInitialized ) THEN
      CheckInputError = 10
    ELSE IF ( D /= D .OR. X /= X .OR.  &
        Yp /= Yp  .OR. Ye /= Ye  .OR. Ym /= Ym ) THEN
      CheckInputError = 11
    ELSE IF ( D < MinD .OR. D > MaxD ) THEN
      CheckInputError = 01
#if EOS_DEBUG
      WRITE(*,*) 'ComputeTemperature ERROR 01'
      WRITE(*,*) 'input D, E/P/S, Ye, Ym :', D, X, Ye, Ym
      WRITE(*,*) 'Table Min/Max D:', MinD, MaxD
#endif
    ELSE IF ( X < MinX .OR. X > MaxX ) THEN
      CheckInputError = 02
#if EOS_DEBUG
      WRITE(*,*) 'ComputeTemperature ERROR 02'
      WRITE(*,*) 'input D, E/P/S, Ye, Ym :', D, X, Ye, Ym
      WRITE(*,*) 'Table Min/Max E/P/S:', MinX, MaxX
#endif
    ELSE IF ( Yp < MinYp .OR. Yp > MaxYp ) THEN
      CheckInputError = 03
#if EOS_DEBUG
      WRITE(*,*) 'ComputeTemperature ERROR 03'
      WRITE(*,*) 'input D, E/P/S, Ye, Ym :', D, X, Ye, Ym
      WRITE(*,*) 'Table Min/Max Yp:', MinYp, MaxYp
#endif
    END IF

  END FUNCTION CheckInputError
  
 SUBROUTINE DescribeEOSComponentsInversionError( Error )

    INTEGER, INTENT(in) :: Error

    CHARACTER(64) :: ErrorString(00:13)

    IF( Error > 13 ) STOP 'ERROR in EOSInversionError flag'

    ErrorString(00) = 'Returned Successfully'
    ErrorString(01) = 'First Argument (D) Outside Table Bounds'
    ErrorString(02) = 'Second Argument (E, P, or S) Outside Table Bounds'
    ErrorString(03) = 'Third Argument (Yp) Outside Table Bounds'
    ErrorString(10) = 'EOS Inversion Not Initialized'
    ErrorString(11) = 'NAN in Argument(s)'
    ErrorString(13) = 'Unable to Find Any Root'

    WRITE(*,*)
    WRITE(*,*) '  wlEOSSplitInversionModule ERROR: ' // TRIM( ErrorString(Error) )
    WRITE(*,*)

  END SUBROUTINE DescribeEOSComponentsInversionError
  
  SUBROUTINE ComputeTemperatureWith_DEYpYl_Single_Guess_Error &
    ( D, E, Ye, Ym, Ds, Ts, Yps, Es, OS, T, T_Guess, Error )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in)  :: D     , E     , Ye, Ym
    REAL(dp), INTENT(in)  :: Ds(1:), Ts(1:), Yps(1:)
    REAL(dp), INTENT(in)  :: Es(1:,1:,1:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(out) :: T
    REAL(dp), INTENT(in)  :: T_Guess
    INTEGER,  INTENT(out) :: Error
    
    T = 0.0_dp
    Error = CheckInputError( D, E, Ye, Ym, MinE, MaxE )
    IF ( Error == 0 ) THEN
      CALL InvertTemperatureWith_DEYpYl_Guess &
             ( D, E, Ye, Ym, Ds, Ts, Yps, Es, &
               OS, T, T_Guess, Error )
    END IF
    
  END SUBROUTINE ComputeTemperatureWith_DEYpYl_Single_Guess_Error

  SUBROUTINE ComputeTemperatureWith_DEYpYl_Single_Guess_NoError &
    ( D, E, Ye, Ym, Ds, Ts, Yps, Es, OS, T, T_Guess )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in)  :: D     , E     , Ye, Ym
    REAL(dp), INTENT(in)  :: Ds(1:), Ts(1:), Yps(1:)
    REAL(dp), INTENT(in)  :: Es(1:,1:,1:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(out) :: T
    REAL(dp), INTENT(in)  :: T_Guess

    INTEGER  :: Error
    
    T = 0.0_dp
    Error = CheckInputError( D, E, Ye, Ym, MinE, MaxE )
    IF ( Error == 0 ) THEN
      CALL InvertTemperatureWith_DEYpYl_Guess &
             ( D, E, Ye, Ym, Ds, Ts, Yps, Es, OS, T, T_Guess, Error )
    END IF

  END SUBROUTINE ComputeTemperatureWith_DEYpYl_Single_Guess_NoError


  SUBROUTINE ComputeTemperatureWith_DEYpYl_Single_NoGuess_Error &
    ( D, E, Ye, Ym, Ds, Ts, Yps, Es, OS, T, Error )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in)    :: D     , E     , Ye, Ym
    REAL(dp), INTENT(in)    :: Ds(1:), Ts(1:), Yps(1:)
    REAL(dp), INTENT(in)    :: Es(1:,1:,1:)
    REAL(dp), INTENT(in)    :: OS
    REAL(dp), INTENT(out)   :: T
    INTEGER,  INTENT(out)   :: Error
    
    T = 0.0_dp
    Error = CheckInputError( D, E, Ye, Ym, MinE, MaxE )
    IF ( Error == 0 ) THEN
      CALL InvertTemperatureWith_DEYpYl_NoGuess &
             ( D, E, Ye, Ym, Ds, Ts, Yps, Es, OS, T, Error )
    END IF

  END SUBROUTINE ComputeTemperatureWith_DEYpYl_Single_NoGuess_Error


  SUBROUTINE ComputeTemperatureWith_DEYpYl_Single_NoGuess_NoError &
    ( D, E, Ye, Ym, Ds, Ts, Yps, Es, OS, T )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in)    :: D     , E     , Ye, Ym
    REAL(dp), INTENT(in)    :: Ds(1:), Ts(1:), Yps(1:)
    REAL(dp), INTENT(in)    :: Es(1:,1:,1:)
    REAL(dp), INTENT(in)    :: OS
    REAL(dp), INTENT(out)   :: T

    INTEGER  :: Error
    
    T = 0.0_dp
    Error = CheckInputError( D, E, Ye, Ym, MinE, MaxE )
    IF ( Error == 0 ) THEN
      CALL InvertTemperatureWith_DEYpYl_NoGuess &
             ( D, E, Ye, Ym, Ds, Ts, Yps, Es, OS, T, Error )
    END IF

  END SUBROUTINE ComputeTemperatureWith_DEYpYl_Single_NoGuess_NoError

  SUBROUTINE ComputeTemperatureWith_DPYpYl_Single_Guess_Error &
    ( D, P, Ye, Ym, Ds, Ts, Yps, Ps, OS, T, T_Guess, Error )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in)  :: D     , P     , Ye, Ym
    REAL(dp), INTENT(in)  :: Ds(1:), Ts(1:), Yps(1:)
    REAL(dp), INTENT(in)  :: Ps(1:,1:,1:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(out) :: T
    REAL(dp), INTENT(in)  :: T_Guess
    INTEGER,  INTENT(out) :: Error
    
    T = 0.0_dp
    Error = CheckInputError( D, P, Ye, Ym, MinP, MaxP )
    IF ( Error == 0 ) THEN
      CALL InvertTemperatureWith_DPYpYl_Guess &
             ( D, P, Ye, Ym, Ds, Ts, Yps, Ps, OS, T, T_Guess, Error )
    END IF

  END SUBROUTINE ComputeTemperatureWith_DPYpYl_Single_Guess_Error


  SUBROUTINE ComputeTemperatureWith_DPYpYl_Single_Guess_NoError &
    ( D, P, Ye, Ym, Ds, Ts, Yps, Ps, OS, T, T_Guess )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in)  :: D     , P     , Ye, Ym
    REAL(dp), INTENT(in)  :: Ds(1:), Ts(1:), Yps(1:)
    REAL(dp), INTENT(in)  :: Ps(1:,1:,1:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(out) :: T
    REAL(dp), INTENT(in)  :: T_Guess

    INTEGER  :: Error
    
    T = 0.0_dp
    Error = CheckInputError( D, P, Ye, Ym, MinP, MaxP )
    IF ( Error == 0 ) THEN
      CALL InvertTemperatureWith_DPYpYl_Guess &
             ( D, P, Ye, Ym, Ds, Ts, Yps, Ps, OS, T, T_Guess, Error )
    END IF

  END SUBROUTINE ComputeTemperatureWith_DPYpYl_Single_Guess_NoError


  SUBROUTINE ComputeTemperatureWith_DPYpYl_Single_NoGuess_Error &
    ( D, P, Ye, Ym, Ds, Ts, Yps, Ps, OS, T, Error )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in)    :: D     , P     , Ye, Ym
    REAL(dp), INTENT(in)    :: Ds(1:), Ts(1:), Yps(1:)
    REAL(dp), INTENT(in)    :: Ps(1:,1:,1:)
    REAL(dp), INTENT(in)    :: OS
    REAL(dp), INTENT(out)   :: T
    INTEGER,  INTENT(out)   :: Error

    T = 0.0_dp
    Error = CheckInputError( D, P, Ye, Ym, MinP, MaxP )
    IF ( Error == 0 ) THEN
      CALL InvertTemperatureWith_DPYpYl_NoGuess &
             ( D, P, Ye, Ym, Ds, Ts, Yps, Ps, OS, T, Error )
    END IF

  END SUBROUTINE ComputeTemperatureWith_DPYpYl_Single_NoGuess_Error


  SUBROUTINE ComputeTemperatureWith_DPYpYl_Single_NoGuess_NoError &
    ( D, P, Ye, Ym, Ds, Ts, Yps, Ps, OS, T )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in)    :: D     , P     , Ye, Ym
    REAL(dp), INTENT(in)    :: Ds(1:), Ts(1:), Yps(1:)
    REAL(dp), INTENT(in)    :: Ps(1:,1:,1:)
    REAL(dp), INTENT(in)    :: OS
    REAL(dp), INTENT(out)   :: T

    INTEGER  :: Error

    T = 0.0_dp
    Error = CheckInputError( D, P, Ye, Ym, MinP, MaxP )
    IF ( Error == 0 ) THEN
      CALL InvertTemperatureWith_DPYpYl_NoGuess &
             ( D, P, Ye, Ym, Ds, Ts, Yps, Ps, OS, T, Error )
    END IF

  END SUBROUTINE ComputeTemperatureWith_DPYpYl_Single_NoGuess_NoError

  SUBROUTINE ComputeTemperatureWith_DSYpYl_Single_Guess_Error &
    ( D, S, Ye, Ym, Ds, Ts, Yps, Ss, OS, T, T_Guess, Error )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in)  :: D     , S     , Ye, Ym
    REAL(dp), INTENT(in)  :: Ds(1:), Ts(1:), Yps(1:)
    REAL(dp), INTENT(in)  :: Ss(1:,1:,1:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(out) :: T
    REAL(dp), INTENT(in)  :: T_Guess
    INTEGER,  INTENT(out) :: Error
    
    T = 0.0_dp
    Error = CheckInputError( D, S, Ye, Ym, MinS, MaxS )
    IF ( Error == 0 ) THEN
      CALL InvertTemperatureWith_DSYpYl_Guess &
             ( D, S, Ye, Ym, Ds, Ts, Yps, Ss, OS, T, T_Guess, Error )
    END IF

  END SUBROUTINE ComputeTemperatureWith_DSYpYl_Single_Guess_Error


  SUBROUTINE ComputeTemperatureWith_DSYpYl_Single_Guess_NoError &
    ( D, S, Ye, Ym, Ds, Ts, Yps, Ss, OS, T, T_Guess )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in)  :: D     , S     , Ye, Ym
    REAL(dp), INTENT(in)  :: Ds(1:), Ts(1:), Yps(1:)
    REAL(dp), INTENT(in)  :: Ss(1:,1:,1:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(out) :: T
    REAL(dp), INTENT(in)  :: T_Guess

    INTEGER  :: Error
    
    T = 0.0_dp
    Error = CheckInputError( D, S, Ye, Ym, MinS, MaxS )
    IF ( Error == 0 ) THEN
      CALL InvertTemperatureWith_DSYpYl_Guess &
             ( D, S, Ye, Ym, Ds, Ts, Yps, Ss, OS, T, T_Guess, Error )
    END IF

  END SUBROUTINE ComputeTemperatureWith_DSYpYl_Single_Guess_NoError


  SUBROUTINE ComputeTemperatureWith_DSYpYl_Single_NoGuess_Error &
    ( D, S, Ye, Ym, Ds, Ts, Yps, Ss, OS, T, Error )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in)    :: D     , S     , Ye, Ym
    REAL(dp), INTENT(in)    :: Ds(1:), Ts(1:), Yps(1:)
    REAL(dp), INTENT(in)    :: Ss(1:,1:,1:)
    REAL(dp), INTENT(in)    :: OS
    REAL(dp), INTENT(out)   :: T
    INTEGER,  INTENT(out)   :: Error
    
    T = 0.0_dp
    Error = CheckInputError( D, S, Ye, Ym, MinS, MaxS )
    IF ( Error == 0 ) THEN
      CALL InvertTemperatureWith_DSYpYl_NoGuess &
             ( D, S, Ye, Ym, Ds, Ts, Yps, Ss, OS, T, Error )
    END IF

  END SUBROUTINE ComputeTemperatureWith_DSYpYl_Single_NoGuess_Error


  SUBROUTINE ComputeTemperatureWith_DSYpYl_Single_NoGuess_NoError &
    ( D, S, Ye, Ym, Ds, Ts, Yps, Ss, OS, T )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in)    :: D     , S     , Ye, Ym
    REAL(dp), INTENT(in)    :: Ds(1:), Ts(1:), Yps(1:)
    REAL(dp), INTENT(in)    :: Ss(1:,1:,1:)
    REAL(dp), INTENT(in)    :: OS
    REAL(dp), INTENT(out)   :: T

    INTEGER  :: Error
    
    T = 0.0_dp
    Error = CheckInputError( D, S, Ye, Ym, MinS, MaxS )
    IF ( Error == 0 ) THEN
      CALL InvertTemperatureWith_DSYpYl_NoGuess &
             ( D, S, Ye, Ym, Ds, Ts, Yps, Ss, OS, T, Error )
    END IF

  END SUBROUTINE ComputeTemperatureWith_DSYpYl_Single_NoGuess_NoError
  
END MODULE wlEOSComponentsSeparateInversionModule
