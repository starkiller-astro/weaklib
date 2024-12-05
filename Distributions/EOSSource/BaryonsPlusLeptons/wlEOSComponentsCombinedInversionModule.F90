MODULE wlEOSComponentsCombinedInversionModule

  USE wlKindModule, ONLY: &
    dp
  USE wlInterpolationModule, ONLY: &
    LogInterpolateSingleVariable_2D_Custom_Point
  USE wlInterpolationUtilitiesModule, ONLY: &
    Index1D_Lin, &
    Index1D_Log
  USE wlMuonEOS, ONLY: &
    MuonStateType, FullMuonEOS
  USE wlElectronPhotonEOS, ONLY: &
    ElectronPhotonStateType, ElectronPhotonEOS
  USE wlLeptonEOSModule, ONLY: &
    HelmholtzTableType, MuonEOSType
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
    !$OMP   InversionInitialized, MinD, MaxD, MinT, MaxT, &
    !$OMP   MinYp, MaxYp, MinE, MaxE, MinP, MaxP, MinS, MaxS )
#elif defined(WEAKLIB_OACC)
    !$ACC DECLARE CREATE( &
    !$ACC   InversionInitialized, MinD, MaxD, MinT, MaxT, &
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

  TYPE(HelmholtzTableType) :: HelmholtzTable 
  TYPE(MuonEOSType) :: MuonTable

  REAL(dp) :: dummy
  
CONTAINS

  ! COME BACK TO THIS. GOTTA DECIDE HOW TO DO THIS. PROBABLY TO MAKE IT EASIER FOR THE USER 
  ! AND CONSISTENT WITH THE OTHER ROUTINES YOU CALL MUON AND ELECTRON EOS HERE, WITH
  ! MAX(YE)=MAX(YP) AND MIN(YE)=MIN(YP). FOR MUONS IT'S A BIT UNCLEAR TO ME, MAYBE DO
  ! MIN(YMU)=0 AND MAX(YMU)=MAX(YP), BUT THE PROBLEM IS THAT LARGE YMU ARE NOT REALLY 
  ! REALIZED, BUT SUCH A LARGE MUON FRACTION WILL CONTRIBUTE SIGNIFICANTLY TO PRESSURE
  ! AND ENERGY. BUT MAYBE TO DETERIMINE MAXIMUM AND MINIMUM OF P,S,E THIS IS NOT SO IMPORTANT
  SUBROUTINE InitializeEOSComponentsInversion( Ds, Ts, Yps, Es, Ps, Ss, &
      HelmholtzTableName, MuonTableName, Verbose_Option )

    REAL(dp), INTENT(in) :: Ds(1:)      , Ts(1:)      , Yps(1:)
    REAL(dp), INTENT(in) :: Es(1:,1:,1:), Ps(1:,1:,1:), Ss(1:,1:,1:)
    CHARACTER(len=*), INTENT(IN) :: HelmholtzTableName, MuonTableName
    LOGICAL,  INTENT(in), OPTIONAL :: Verbose_Option
    
    LOGICAL :: Verbose

    IF( PRESENT( Verbose_Option ) )THEN
      Verbose = Verbose_Option
    ELSE
      Verbose = .FALSE.
    END IF

    ! Initialize Helmholtz and Muon Local Tables
    CALL ReadHelmholtzTableHDF( HelmholtzTable, HelmholtzTableName )
    CALL ReadMuonTableHDF( MuonTable, MuonTableName )
    
  ! For maximum and minimum you only need to calculate 
    ! maximum and minimum Ye and Ym
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
      WRITE(*,'(A6,A24,2ES11.3E3)') '', 'Min/Max Yp            = ', MinYp, MaxYp
      WRITE(*,'(A6,A24,2ES11.3E3)') '', 'Min/Max E  [erg g^-1] = ', MinE, MaxE
      WRITE(*,'(A6,A24,2ES11.3E3)') '', 'Min/Max P [dyn cm^-2] = ', MinP, MaxP
      WRITE(*,'(A6,A24,2ES11.3E3)') '', 'Min/Max S       [k_B] = ', MinS, MaxS
      WRITE(*,*)
    END IF

    InversionComponentsInitialized = .TRUE.

#if defined(WEAKLIB_OMP_OL)
    !$OMP TARGET UPDATE TO( &
    !$OMP   InversionInitialized, MinD, MaxD, MinT, MaxT, &
    !$OMP   MinYp, MaxYp, MinE, MaxE, MinP, MaxP, MinS, MaxS )
#elif defined(WEAKLIB_OACC)
    !$ACC UPDATE DEVICE( &
    !$ACC   InversionInitialized, MinD, MaxD, MinT, MaxT, &
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
  
  REAL(dp) FUNCTION InverseLogInterp( x_a, x_b, y_a, y_b, Yp, OS )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in) :: x_a, x_b, y_a, y_b, Yp, OS

    InverseLogInterp &
      = 10.0_dp**( LOG10( x_a ) + LOG10( x_b/x_a ) &
                 * LOG10( (Yp+OS)/(y_a+OS) ) / LOG10( (y_b+OS)/(y_a+OS) ) )

    RETURN
  END FUNCTION InverseLogInterp
  

  SUBROUTINE ComputeTemperatureWith_DXYpYl_Guess &
    ( D, X, Ye, Ym, Ds, Ts, Yps, Xs, InputE, InputP, InputS, OS, T, T_Guess, &
    Error )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in)  :: D     , X     , Ye, Ym
    REAL(dp), INTENT(in)  :: Ds(1:), Ts(1:), Yps(1:)
    REAL(dp), INTENT(in)  :: Xs(1:,1:,1:)
    REAL(dp), INTENT(in)  :: OS
    REAL(dp), INTENT(in)  :: InputE, InputP, InputS
    ! CHARACTER(len=*), INTENT(in) :: X_name
    REAL(dp), INTENT(out) :: T
    REAL(dp), INTENT(in)  :: T_Guess
    INTEGER,  INTENT(out) :: Error

    REAL(dp) :: Yp
    INTEGER  :: iD, iT, iYp, iL_D, iL_Y
    INTEGER  :: SizeDs, SizeTs, SizeYps
    INTEGER  :: d_c, d_i
    INTEGER  :: i_a, i_b, i_c, i
    REAL(dp) :: T_a, T_b, T_c, T_i
    REAL(dp) :: X_a, X_b, X_c, X_i
    REAL(dp) :: f_a, f_b, f_c
    REAL(dp) :: LogD
    REAL(dp) :: LogDs_i(2), Yps_i(2)
    REAL(dp) :: Xs_a(2,2), Xs_b(2,2), Xs_c(2,2), Xs_i(2,2)
    REAL(dp) :: tBegin, tEnd
    REAL(dp) :: E_LeptPhot, P_LeptPhot, S_LeptPhot
    
    LOGICAL :: InvertP, InvertE, InvertS
    REAL(dp) :: Is_this_P, Is_this_E, Is_this_S
    
    ! Electron and Muon quantities
    TYPE(ElectronPhotonStateType) :: ElectronPhotonState
    TYPE(MuonStateType) :: MuonState
    
    ! Make sure that Yp = Ye + Ym also at the table level
    REAL(dp) :: Ye_over_Yp, Ym_over_Yp

    Yp = Ye + Ym
    
    Ye_over_Yp = Ye/Yp
    Ym_over_Yp = Ym/Yp

    T = 0.0_dp
    Error = 0

    LogD = LOG10( D )

    SizeDs = SIZE( Ds )
    SizeTs = SIZE( Ts )
    SizeYps = SIZE( Yps )

    iD = Index1D_Log( D,       Ds )
    iT = Index1D_Log( T_Guess, Ts )
    iYp = Index1D_Lin( Yp,       Yps )

    iD = MIN( MAX( 1, iD ), SizeDs - 1 )
    iT = MIN( MAX( 1, iT ), SizeTs - 1 )
    iYp = MIN( MAX( 1, iYp ), SizeYps - 1 )

    LogDs_i = LOG10( Ds(iD:iD+1) )
    Yps_i = Yps(iYp:iYp+1)
    ! -------------------------------------------------------------------

    ! --- First Check if Initial Guess Gives a Solution ---
    i_a = iT
    T_a = Ts(i_a)
    Xs_a = 10.0d0**Xs(iD:iD+1,i_a,iYp:iYp+1) - OS
    ! Calculate electron and muon contribution
    DO iL_D=1,2
      DO iL_Y=1,2
        ElectronPhotonState % t = T_a
        ElectronPhotonState % rho = Ds(iD+iL_D-1)
        ElectronPhotonState % ye = Yps(iYp+iL_Y-1) * Ye_over_Yp
        
        ! CALL FullHelmEOS(1, HelmholtzTable, ElectronPhotonState, .false., .false.)
        CALL ElectronPhotonEOS(HelmholtzTable, ElectronPhotonState)

        MuonState % t = T_a
        MuonState % rhoym = Ds(iD+iL_D-1) * Yps(iYp+iL_Y-1) * Ym_over_Yp
        
        CALL FullMuonEOS(MuonTable, MuonState)
        
        E_LeptPhot = ElectronPhotonState % e + MuonState % e 
        P_LeptPhot = ElectronPhotonState % p + MuonState % p
        S_LeptPhot = ElectronPhotonState % s + MuonState % s
        
        Xs_a(iL_D,iL_Y) = Xs_a(iL_D,iL_Y) + InputE*E_LeptPhot + InputP*P_LeptPhot + InputS*S_LeptPhot

      END DO
    END DO
    Xs_a = LOG10(Xs_a+OS)
    
    CALL LogInterpolateSingleVariable_2D_Custom_Point &
           ( LogD, Yp, LogDs_i, Yps_i, OS, Xs_a, X_a )

    f_a = X - X_a
    
    i_b = i_a + 1
    T_b = Ts(i_b)
    Xs_b = 10.0d0**Xs(iD:iD+1,i_b,iYp:iYp+1) - OS
    ! Calculate electron and muon contribution
    DO iL_D=1,2
      DO iL_Y=1,2
        ElectronPhotonState % t = T_b
        ElectronPhotonState % rho = Ds(iD+iL_D-1)
        ElectronPhotonState % ye = Yps(iYp+iL_Y-1) * Ye_over_Yp
        
        ! CALL FullHelmEOS(1, HelmholtzTable, ElectronPhotonState, .false., .false.)
        CALL ElectronPhotonEOS(HelmholtzTable, ElectronPhotonState)
          
        MuonState % t = T_b
        MuonState % rhoym = Ds(iD+iL_D-1) * Yps(iYp+iL_Y-1) * Ym_over_Yp
        
        CALL FullMuonEOS(MuonTable, MuonState)

        E_LeptPhot = ElectronPhotonState % e + MuonState % e 
        P_LeptPhot = ElectronPhotonState % p + MuonState % p
        S_LeptPhot = ElectronPhotonState % s + MuonState % s
        
        Xs_b(iL_D,iL_Y) = Xs_b(iL_D,iL_Y) + InputE*E_LeptPhot + InputP*P_LeptPhot + InputS*S_LeptPhot

      END DO
    END DO
    Xs_b = LOG10(Xs_b+OS)

    CALL LogInterpolateSingleVariable_2D_Custom_Point &
           ( LogD, Yp, LogDs_i, Yps_i, OS, Xs_b, X_b )
    
    f_b = X - X_b
    
    IF ( f_a * f_b <= 0.0_dp ) THEN
      T = InverseLogInterp( T_a, T_b, X_a, X_b, X, OS )
      RETURN
    END IF

    ! -------------------------------------------------------------------

    i_a = 1
    T_a = Ts(i_a)
    Xs_a = 10.0d0**Xs(iD:iD+1,i_a,iYp:iYp+1) - OS
    ! Calculate electron and muon contribution
    DO iL_D=1,2
      DO iL_Y=1,2
        ElectronPhotonState % t = T_a
        ElectronPhotonState % rho = Ds(iD+iL_D-1)
        ElectronPhotonState % ye = Yps(iYp+iL_Y-1) * Ye_over_Yp
        CALL ElectronPhotonEOS(HelmholtzTable, ElectronPhotonState)

        MuonState % t = T_a
        MuonState % rhoym = Ds(iD+iL_D-1) * Yps(iYp+iL_Y-1) * Ym_over_Yp
        CALL FullMuonEOS(MuonTable, MuonState)
          
        E_LeptPhot = ElectronPhotonState % e + MuonState % e 
        P_LeptPhot = ElectronPhotonState % p + MuonState % p
        S_LeptPhot = ElectronPhotonState % s + MuonState % s
        
        Xs_a(iL_D,iL_Y) = Xs_a(iL_D,iL_Y) + InputE*E_LeptPhot + InputP*P_LeptPhot + InputS*S_LeptPhot

      END DO
    END DO
    Xs_a = LOG10(Xs_a+OS)

    CALL LogInterpolateSingleVariable_2D_Custom_Point &
           ( LogD, Yp, LogDs_i, Yps_i, OS, Xs_a, X_a )

    f_a = X - X_a
    
    i_b = SizeTs
    T_b = Ts(i_b)
    Xs_b = 10.0d0**Xs(iD:iD+1,i_b,iYp:iYp+1) - OS
    ! Calculate electron and muon contribution
    DO iL_D=1,2
      DO iL_Y=1,2
        ElectronPhotonState % t = T_b
        ElectronPhotonState % rho = Ds(iD+iL_D-1)
        ElectronPhotonState % ye = Yps(iYp+iL_Y-1) * Ye_over_Yp
        CALL ElectronPhotonEOS(HelmholtzTable, ElectronPhotonState)

        MuonState % t = T_b
        MuonState % rhoym = Ds(iD+iL_D-1) * Yps(iYp+iL_Y-1) * Ym_over_Yp
        CALL FullMuonEOS(MuonTable, MuonState)

        E_LeptPhot = ElectronPhotonState % e + MuonState % e 
        P_LeptPhot = ElectronPhotonState % p + MuonState % p
        S_LeptPhot = ElectronPhotonState % s + MuonState % s
        
        Xs_b(iL_D,iL_Y) = Xs_b(iL_D,iL_Y) + InputE*E_LeptPhot + InputP*P_LeptPhot + InputS*S_LeptPhot

      END DO
    END DO
    Xs_b = LOG10(Xs_b+OS)

    CALL LogInterpolateSingleVariable_2D_Custom_Point &
           ( LogD, Yp, LogDs_i, Yps_i, OS, Xs_b, X_b )

    f_b = X - X_b

    IF ( f_a * f_b <= 0.0_dp ) THEN

      i_c = MIN( MAX( i_a + 1, iT ), i_b - 1 )
      DO WHILE ( i_b > i_a + 1 )

        T_c = Ts(i_c)
        Xs_c = 10.0d0**Xs(iD:iD+1,i_c,iYp:iYp+1) - OS
        ! Calculate electron and muon contribution
        DO iL_D=1,2
          DO iL_Y=1,2
            ElectronPhotonState % t = T_c
            ElectronPhotonState % rho = Ds(iD+iL_D-1)
            ElectronPhotonState % ye = Yps(iYp+iL_Y-1) * Ye_over_Yp
            CALL ElectronPhotonEOS(HelmholtzTable, ElectronPhotonState)

            MuonState % t = T_c
            MuonState % rhoym = Ds(iD+iL_D-1) * Yps(iYp+iL_Y-1) * Ym_over_Yp
            CALL FullMuonEOS(MuonTable, MuonState)

            E_LeptPhot = ElectronPhotonState % e + MuonState % e 
            P_LeptPhot = ElectronPhotonState % p + MuonState % p
            S_LeptPhot = ElectronPhotonState % s + MuonState % s
            
            Xs_c(iL_D,iL_Y) = Xs_c(iL_D,iL_Y) + InputE*E_LeptPhot + InputP*P_LeptPhot + InputS*S_LeptPhot

          END DO
        END DO
        Xs_c = LOG10(Xs_c+OS)
    
        CALL LogInterpolateSingleVariable_2D_Custom_Point &
               ( LogD, Yp, LogDs_i, Yps_i, OS, Xs_c, X_c )

        f_c = X - X_c

        IF ( f_a * f_c <= 0.0_dp ) THEN
          i_b = i_c
          T_b = T_c
          X_b = X_c
          f_b = f_c
        ELSE
          i_a = i_c
          T_a = T_c
          X_a = X_c
          f_a = f_c
        END IF

        i_c = MAX( i_a + 1, ( i_a + i_b ) / 2 )

      END DO

    ELSE

      ! --- Pick Root Nearest Initial Guess ---
      d_i = SizeTs
      i_c = MIN( MAX( 1, iT ), SizeTs - 1 )
      T_c = T_a
      X_c = X_a
      DO i = 2, SizeTs - 1

        T_i = Ts(i)
        Xs_i = 10.0d0**Xs(iD:iD+1,i,iYp:iYp+1) - OS
        ! Calculate electron and muon contribution
        DO iL_D=1,2
          DO iL_Y=1,2
            ElectronPhotonState % t = T_i
            ElectronPhotonState % rho = Ds(iD+iL_D-1)
            ElectronPhotonState % ye = Yps(iYp+iL_Y-1) * Ye_over_Yp
            CALL ElectronPhotonEOS(HelmholtzTable, ElectronPhotonState)

            MuonState % t = T_i
            MuonState % rhoym = Ds(iD+iL_D-1) * Yps(iYp+iL_Y-1) * Ym_over_Yp
            CALL FullMuonEOS(MuonTable, MuonState)
  
            E_LeptPhot = ElectronPhotonState % e + MuonState % e 
            P_LeptPhot = ElectronPhotonState % p + MuonState % p
            S_LeptPhot = ElectronPhotonState % s + MuonState % s
            
            Xs_i(iL_D,iL_Y) = Xs_i(iL_D,iL_Y) + InputE*E_LeptPhot + InputP*P_LeptPhot + InputS*S_LeptPhot

          END DO
        END DO
        Xs_i = LOG10(Xs_i+OS)

        CALL LogInterpolateSingleVariable_2D_Custom_Point &
               ( LogD, Yp, LogDs_i, Yps_i, OS, Xs_i, X_i )

        f_c = X - X_c
        f_b = X - X_i

        d_c = ABS( i - i_c )

        IF ( f_c * f_b <= 0.0_dp .AND. d_c < d_i ) THEN
          d_i = d_c
          i_a = i - 1
          T_a = T_c
          X_a = X_c
          i_b = i
          T_b = T_i
          X_b = X_i
        ELSE
          T_c = T_i
          X_c = X_i
        END IF

      END DO

      IF ( d_i >= SizeTs ) THEN
        Error = 13
        WRITE(*,*) X, X_a, X_b, i_a, i_b
      ENDIF
    END IF

    IF ( Error .NE. 0 ) THEN
      T = 0.0_dp
    ELSE
      T = InverseLogInterp( T_a, T_b, X_a, X_b, X, OS )
    END IF

  END SUBROUTINE ComputeTemperatureWith_DXYpYl_Guess
  
  SUBROUTINE ComputeTemperatureWith_DXYpYl_NoGuess &
    ( D, X, Ye, Ym, Ds, Ts, Yps, Xs, InputE, InputP, InputS, OS, T, &
    Error )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(dp), INTENT(in)  :: D     , X     , Ye, Ym
    REAL(dp), INTENT(in)  :: Ds(1:), Ts(1:), Yps(1:)
    REAL(dp), INTENT(in)  :: Xs(1:,1:,1:)
    REAL(dp), INTENT(in)  :: OS
    ! CHARACTER(len=*), INTENT(in) :: X_name
    REAL(dp), INTENT(in)  :: InputE, InputP, InputS
    REAL(dp), INTENT(out) :: T
    INTEGER,  INTENT(out) :: Error

    REAL(dp) :: Yp
    INTEGER  :: iD, iT, iYp, iL_D, iL_Y
    INTEGER  :: SizeDs, SizeTs, SizeYps
    INTEGER  :: d_c, d_i
    INTEGER  :: i_a, i_b, i_c, i
    REAL(dp) :: T_a, T_b, T_c, T_i
    REAL(dp) :: X_a, X_b, X_c, X_i
    REAL(dp) :: f_a, f_b, f_c
    REAL(dp) :: LogD
    REAL(dp) :: LogDs_i(2), Yps_i(2)
    REAL(dp) :: Xs_a(2,2), Xs_b(2,2), Xs_c(2,2), Xs_i(2,2)
    
    REAL(dp) :: tBegin, tEnd
    REAL(dp) :: E_LeptPhot, P_LeptPhot, S_LeptPhot
    
    LOGICAL :: InvertP, InvertE, InvertS
    REAL(dp) :: Is_this_P, Is_this_E, Is_this_S
    
    ! Electron and Muon quantities
    TYPE(ElectronPhotonStateType) :: ElectronPhotonState
    TYPE(MuonStateType) :: MuonState
    
    ! Make sure that Yp = Ye + Ym also at the table level
    REAL(dp) :: Ye_over_Yp, Ym_over_Yp
    
    Yp = Ye + Ym

    Ye_over_Yp = Ye/Yp
    Ym_over_Yp = Ym/Yp

    ! -------------------------------------------------------------------
    Error = 0

    LogD = LOG10( D )

    SizeDs = SIZE( Ds )
    SizeTs = SIZE( Ts )
    SizeYps = SIZE( Yps )

    iD = Index1D_Log( D, Ds )
    iYp = Index1D_Lin( Yp, Yps )

    iD = MIN( MAX( 1, iD ), SizeDs - 1 )
    iYp = MIN( MAX( 1, iYp ), SizeYps - 1 )

    LogDs_i = LOG10( Ds(iD:iD+1) )
    Yps_i = Yps(iYp:iYp+1)

    ! -------------------------------------------------------------------

    i_a = 1
    T_a = Ts(i_a)
    Xs_a = 10.0d0**Xs(iD:iD+1,i_a,iYp:iYp+1) - OS
    ! Calculate electron and muon contribution
    DO iL_D=1,2
      DO iL_Y=1,2
        ElectronPhotonState % t = T_a
        ElectronPhotonState % rho = Ds(iD+iL_D-1)
        ElectronPhotonState % ye = Yps(iYp+iL_Y-1) * Ye_over_Yp
        CALL ElectronPhotonEOS(HelmholtzTable, ElectronPhotonState)
          
        MuonState % t = T_a
        MuonState % rhoym = Ds(iD+iL_D-1) * Yps(iYp+iL_Y-1) * Ym_over_Yp       
        CALL FullMuonEOS(MuonTable, MuonState)
          
        E_LeptPhot = ElectronPhotonState % e + MuonState % e 
        P_LeptPhot = ElectronPhotonState % p + MuonState % p
        S_LeptPhot = ElectronPhotonState % s + MuonState % s
        
        Xs_a(iL_D,iL_Y) = Xs_a(iL_D,iL_Y) + InputE*E_LeptPhot + InputP*P_LeptPhot + InputS*S_LeptPhot
        
      END DO
    END DO
    Xs_a = LOG10(Xs_a+OS)

    CALL LogInterpolateSingleVariable_2D_Custom_Point &
           ( LogD, Yp, LogDs_i, Yps_i, OS, Xs_a, X_a )

    f_a = X - X_a

    i_b = SizeTs
    T_b = Ts(i_b)
    Xs_b = 10.0d0**Xs(iD:iD+1,i_b,iYp:iYp+1) - OS
    ! Calculate electron and muon contribution
    DO iL_D=1,2
      DO iL_Y=1,2
        ElectronPhotonState % t = T_b
        ElectronPhotonState % rho = Ds(iD+iL_D-1)
        ElectronPhotonState % ye = Yps(iYp+iL_Y-1) * Ye_over_Yp       
        CALL ElectronPhotonEOS(HelmholtzTable, ElectronPhotonState)

        MuonState % t = T_b
        MuonState % rhoym = Ds(iD+iL_D-1) * Yps(iYp+iL_Y-1) * Ym_over_Yp        
        CALL FullMuonEOS(MuonTable, MuonState)
          
        E_LeptPhot = ElectronPhotonState % e + MuonState % e 
        P_LeptPhot = ElectronPhotonState % p + MuonState % p
        S_LeptPhot = ElectronPhotonState % s + MuonState % s
        
        Xs_b(iL_D,iL_Y) = Xs_b(iL_D,iL_Y) + InputE*E_LeptPhot + InputP*P_LeptPhot + InputS*S_LeptPhot

      END DO
    END DO
    Xs_b = LOG10(Xs_b+OS)

    CALL LogInterpolateSingleVariable_2D_Custom_Point &
           ( LogD, Yp, LogDs_i, Yps_i, OS, Xs_b, X_b )

    f_b = X - X_b

    IF ( f_a * f_b <= 0.0_dp ) THEN

      DO WHILE ( i_b > i_a + 1 )

        i_c = MAX( i_a + 1, ( i_a + i_b ) / 2 )
        T_c = Ts(i_c)
        Xs_c = 10.0d0**Xs(iD:iD+1,i_c,iYp:iYp+1) - OS
        ! Calculate electron and muon contribution
        DO iL_D=1,2
          DO iL_Y=1,2
            ElectronPhotonState % t = T_c
            ElectronPhotonState % rho = Ds(iD+iL_D-1)
            ElectronPhotonState % ye = Yps(iYp+iL_Y-1) * Ye_over_Yp
            CALL ElectronPhotonEOS(HelmholtzTable, ElectronPhotonState)

            MuonState % t = T_c
            MuonState % rhoym = Ds(iD+iL_D-1) * Yps(iYp+iL_Y-1) * Ym_over_Yp
            CALL FullMuonEOS(MuonTable, MuonState)
              
            E_LeptPhot = ElectronPhotonState % e + MuonState % e 
            P_LeptPhot = ElectronPhotonState % p + MuonState % p
            S_LeptPhot = ElectronPhotonState % s + MuonState % s
            
            Xs_c(iL_D,iL_Y) = Xs_c(iL_D,iL_Y) + InputE*E_LeptPhot + InputP*P_LeptPhot + InputS*S_LeptPhot

          END DO
        END DO
        Xs_c = LOG10(Xs_c+OS)

        CALL LogInterpolateSingleVariable_2D_Custom_Point &
               ( LogD, Yp, LogDs_i, Yps_i, OS, Xs_c, X_c )

        f_c = X - X_c

        IF ( f_a * f_c <= 0.0_dp ) THEN
          i_b = i_c
          T_b = T_c
          X_b = X_c
          f_b = f_c
        ELSE
          i_a = i_c
          T_a = T_c
          X_a = X_c
          f_a = f_c
        END IF

      END DO

    ELSE

      ! --- Pick Highest Temperature Root ---
      DO i = SizeTs - 1, 2, -1

        T_i = Ts(i)
        Xs_i = 10.0d0**Xs(iD:iD+1,i,iYp:iYp+1) - OS
        ! Calculate electron and muon contribution
        DO iL_D=1,2
          DO iL_Y=1,2
            ElectronPhotonState % t = T_i
            ElectronPhotonState % rho = Ds(iD+iL_D-1)
            ElectronPhotonState % ye = Yps(iYp+iL_Y-1) * Ye_over_Yp
            CALL ElectronPhotonEOS(HelmholtzTable, ElectronPhotonState)

            MuonState % t = T_i
            MuonState % rhoym = Ds(iD+iL_D-1) * Yps(iYp+iL_Y-1) * Ym_over_Yp
            CALL FullMuonEOS(MuonTable, MuonState)
              
            E_LeptPhot = ElectronPhotonState % e + MuonState % e 
            P_LeptPhot = ElectronPhotonState % p + MuonState % p
            S_LeptPhot = ElectronPhotonState % s + MuonState % s
            
            Xs_i(iL_D,iL_Y) = Xs_i(iL_D,iL_Y) + InputE*E_LeptPhot + InputP*P_LeptPhot + InputS*S_LeptPhot

          END DO
        END DO
        Xs_i = LOG10(Xs_i+OS)
        
        CALL LogInterpolateSingleVariable_2D_Custom_Point &
               ( LogD, Yp, LogDs_i, Yps_i, OS, Xs_i, X_i )

        f_a = X - X_i
        f_b = X - X_b

        IF ( f_a * f_b <= 0.0_dp ) THEN
          i_a = i
          T_a = T_i
          X_a = X_i
          EXIT
        ELSE
          i_b = i
          T_b = T_i
          X_b = X_i
        END IF

      END DO

      f_a = X - X_a
      f_b = X - X_b

      IF ( f_a * f_b > 0.0_dp ) Error = 13

    END IF

    IF( Error .NE. 0 )THEN
      T = 0.0_dp
    ELSE
      T = InverseLogInterp( T_a, T_b, X_a, X_b, X, OS )
    END IF

  END SUBROUTINE ComputeTemperatureWith_DXYpYl_NoGuess

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
    
    REAL(DP) :: InputE, InputP, InputS

    InputE = 1.0d0
    InputP  = 0.0d0
    InputS  = 0.0d0
    
    T = 0.0_dp
    Error = CheckInputError( D, E, Ye, Ym, MinE, MaxE )
    IF ( Error == 0 ) THEN
      CALL ComputeTemperatureWith_DXYpYl_Guess &
             ( D, E, Ye, Ym, Ds, Ts, Yps, Es, &
             InputE, InputP, InputS, OS, T, T_Guess, Error )
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
    
    REAL(DP) :: InputE, InputP, InputS

    InputE = 1.0d0
    InputP  = 0.0d0
    InputS  = 0.0d0

    T = 0.0_dp
    Error = CheckInputError( D, E, Ye, Ym, MinE, MaxE )
    IF ( Error == 0 ) THEN
      CALL ComputeTemperatureWith_DXYpYl_Guess &
             ( D, E, Ye, Ym, Ds, Ts, Yps, Es, &
             InputE, InputP, InputS, Os, T, T_Guess, Error )
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
    
    REAL(DP) :: InputE, InputP, InputS

    InputE = 1.0d0
    InputP  = 0.0d0
    InputS  = 0.0d0
    
    T = 0.0_dp
    Error = CheckInputError( D, E, Ye, Ym, MinE, MaxE )
    IF ( Error == 0 ) THEN
      CALL ComputeTemperatureWith_DXYpYl_NoGuess &
             ( D, E, Ye, Ym, Ds, Ts, Yps, Es, &
             InputE, InputP, InputS, Os, T, Error )
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
    
    REAL(DP) :: InputE, InputP, InputS

    InputE = 1.0d0
    InputP  = 0.0d0
    InputS  = 0.0d0
    
    T = 0.0_dp
    Error = CheckInputError( D, E, Ye, Ym, MinE, MaxE )
    IF ( Error == 0 ) THEN
      CALL ComputeTemperatureWith_DXYpYl_NoGuess &
             ( D, E, Ye, Ym, Ds, Ts, Yps, Es, &
             InputE, InputP, InputS, Os, T, Error )
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
    
    REAL(DP) :: InputE, InputP, InputS

    InputE = 0.0d0
    InputP  = 1.0d0
    InputS  = 0.0d0
    
    T = 0.0_dp
    Error = CheckInputError( D, P, Ye, Ym, MinP, MaxP )
    IF ( Error == 0 ) THEN
      CALL ComputeTemperatureWith_DXYpYl_Guess &
             ( D, P, Ye, Ym, Ds, Ts, Yps, Ps, &
             InputE, InputP, InputS, Os, T, T_Guess, Error )
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
    
    REAL(DP) :: InputE, InputP, InputS

    InputE = 0.0d0
    InputP  = 1.0d0
    InputS  = 0.0d0
    
    T = 0.0_dp
    Error = CheckInputError( D, P, Ye, Ym, MinP, MaxP )
    IF ( Error == 0 ) THEN
      CALL ComputeTemperatureWith_DXYpYl_Guess &
             ( D, P, Ye, Ym, Ds, Ts, Yps, Ps, &
             InputE, InputP, InputS, Os, T, T_Guess, Error )
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
    
    REAL(DP) :: InputE, InputP, InputS

    InputE = 0.0d0
    InputP  = 1.0d0
    InputS  = 0.0d0
    
    T = 0.0_dp
    Error = CheckInputError( D, P, Ye, Ym, MinP, MaxP )
    IF ( Error == 0 ) THEN
      CALL ComputeTemperatureWith_DXYpYl_NoGuess &
             ( D, P, Ye, Ym, Ds, Ts, Yps, Ps, &
             InputE, InputP, InputS, Os, T, Error )
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
    
    REAL(DP) :: InputE, InputP, InputS

    InputE = 0.0d0
    InputP  = 1.0d0
    InputS  = 0.0d0
    
    T = 0.0_dp
    Error = CheckInputError( D, P, Ye, Ym, MinP, MaxP )
    IF ( Error == 0 ) THEN
      CALL ComputeTemperatureWith_DXYpYl_NoGuess &
             ( D, P, Ye, Ym, Ds, Ts, Yps, Ps, &
             InputE, InputP, InputS, Os, T, Error )
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
    
    REAL(DP) :: InputE, InputP, InputS

    InputE = 0.0d0
    InputP  = 0.0d0
    InputS  = 1.0d0
    
    T = 0.0_dp
    Error = CheckInputError( D, S, Ye, Ym, MinS, MaxS )
    IF ( Error == 0 ) THEN
      CALL ComputeTemperatureWith_DXYpYl_Guess &
             ( D, S, Ye, Ym, Ds, Ts, Yps, Ss, &
             InputE, InputP, InputS, OS, T, T_Guess, Error )
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
    
    REAL(DP) :: InputE, InputP, InputS

    InputE = 0.0d0
    InputP  = 0.0d0
    InputS  = 1.0d0
    
    T = 0.0_dp
    Error = CheckInputError( D, S, Ye, Ym, MinS, MaxS )
    IF ( Error == 0 ) THEN
      CALL ComputeTemperatureWith_DXYpYl_Guess &
             ( D, S, Ye, Ym, Ds, Ts, Yps, Ss, &
             InputE, InputP, InputS, Os, T, T_Guess, Error )
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
    
    REAL(DP) :: InputE, InputP, InputS

    InputE = 0.0d0
    InputP  = 0.0d0
    InputS  = 1.0d0
    
    T = 0.0_dp
    Error = CheckInputError( D, S, Ye, Ym, MinS, MaxS )
    IF ( Error == 0 ) THEN
      CALL ComputeTemperatureWith_DXYpYl_NoGuess &
             ( D, S, Ye, Ym, Ds, Ts, Yps, Ss, &
             InputE, InputP, InputS, Os, T, Error )
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
    
    REAL(DP) :: InputE, InputP, InputS

    InputE = 0.0d0
    InputP  = 0.0d0
    InputS  = 1.0d0
    
    T = 0.0_dp
    Error = CheckInputError( D, S, Ye, Ym, MinS, MaxS )
    IF ( Error == 0 ) THEN
      CALL ComputeTemperatureWith_DXYpYl_NoGuess &
             ( D, S, Ye, Ym, Ds, Ts, Yps, Ss, &
             InputE, InputP, InputS, Os, T, Error )
    END IF

  END SUBROUTINE ComputeTemperatureWith_DSYpYl_Single_NoGuess_NoError
  
END MODULE wlEOSComponentsCombinedInversionModule
