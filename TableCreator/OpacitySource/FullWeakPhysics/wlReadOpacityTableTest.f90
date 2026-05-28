PROGRAM wlReadOpacityTableTest

  USE wlIOModuleHDF, ONLY: &
    InitializeHDF, &
    FinalizeHDF
  USE wlOpacityTableModule, ONLY: &
    OpacityTableType, &
    DescribeOpacityTable
  USE wlOpacityTableIOModuleHDF, ONLY: &
    ReadOpacityTableHDF
  
  IMPLICIT NONE

  TYPE(OpacityTableType), TARGET :: OpacityTable

  CALL InitializeHDF( )

  CALL ReadOpacityTableHDF &
          ( OpacityTable, &
           FileName_EmAb_Option &
           = "wl-Op-SFHo-15-25-50-E40-EmAb.h5", &
           FileName_Iso_Option &
           = "wl-Op-SFHo-15-25-50-E40-Iso.h5", &
           FileName_NNS_Option &
           = "wl-Op-SFHo-15-25-50-E40-NNS.h5", &
           FileName_NES_Option &
           = "wl-Op-SFHo-15-25-50-E40-NES.h5", &
           FileName_Pair_Option &
           = "wl-Op-SFHo-15-25-50-E40-Pair.h5", &
           FileName_Brem_Option &
           = "wl-Op-SFHo-15-25-50-E40-Brem.h5", &
           EquationOfStateTableName_Option &
           = "wl-EOS-SFHo-15-25-50.h5", &
           Verbose_Option = .true. )

  CALL DescribeOpacityTable( OpacityTable )

  CALL FinalizeHDF( ) 

  !-- Test Interpolate_NNS_Point

  CALL random_seed ( )

  TestPoint: BLOCK

    USE wlKindModule, ONLY: &
      dp
    USE wlExtPhysicalConstantsModule, ONLY: &
      kMeV
    USE wlOpacityFieldsModule, ONLY: &
      iNu_NNS, iNuBar_NNS, &
      iNeutron_NNS, iProton_NNS

    REAL(dp) :: random
    REAL(dp) :: T, MuB
    REAL(dp), DIMENSION ( :, : ), ALLOCATABLE :: &
      Interpolated_Nu_N_0,  Interpolated_Nu_N_1,  &
      Interpolated_NuB_N_0, Interpolated_NuB_N_1, &
      Interpolated_Nu_P_0,  Interpolated_Nu_P_1,  &
      Interpolated_NuB_P_0, Interpolated_NuB_P_1
    REAL(dp), DIMENSION ( :, : ), ALLOCATABLE :: &
      Computed_Nu_N_0,  Computed_Nu_N_1,  &
      Computed_NuB_N_0, Computed_NuB_N_1, &
      Computed_Nu_P_0,  Computed_Nu_P_1,  &
      Computed_NuB_P_0, Computed_NuB_P_1

    associate &
      ( iT_TS     =>  OpacityTable % TS % Indices % iT, &
        nPointsE  =>  OpacityTable % Scat_NNS % nPoints ( 1 ) )
    associate &
      ( T_Min    =>  OpacityTable % TS % minValues ( iT_TS ), &
        T_Max    =>  OpacityTable % TS % maxValues ( iT_TS ), &
        MuB_Min  =>  OpacityTable % MuBGrid % minValue, &
        MuB_Max  =>  OpacityTable % MuBGrid % maxValue )    

    CALL random_number ( random )
    T  =  10.d0 ** ( LOG10 ( T_Min )  &
                     +  random  *  ( LOG10 ( T_Max )  -  LOG10 ( T_Min ) ) ) 

    CALL random_number ( random )
    MuB  =  MuB_Min  +  random * ( MuB_Max - MuB_Min )

    WRITE (*,*)
    WRITE (*,'(A6,ES12.6E2)') 'T   = ', T 
    WRITE (*,'(A6,ES12.6E2)') 'MuB = ', MuB 
!    WRITE (*,'(A11,I4.4)')    'nPointsE = ', nPointsE 

    allocate &
      ( Interpolated_Nu_N_0  ( nPointsE, nPointsE ), &
        Interpolated_Nu_N_1  ( nPointsE, nPointsE ), &
        Interpolated_NuB_N_0 ( nPointsE, nPointsE ), &
        Interpolated_NuB_N_1 ( nPointsE, nPointsE ), &
        Interpolated_Nu_P_0  ( nPointsE, nPointsE ), &
        Interpolated_Nu_P_1  ( nPointsE, nPointsE ), &
        Interpolated_NuB_P_0 ( nPointsE, nPointsE ), &
        Interpolated_NuB_P_1 ( nPointsE, nPointsE ) )
    allocate &
      ( Computed_Nu_N_0  ( nPointsE, nPointsE ), &
        Computed_Nu_N_1  ( nPointsE, nPointsE ), &
        Computed_NuB_N_0 ( nPointsE, nPointsE ), &
        Computed_NuB_N_1 ( nPointsE, nPointsE ), &
        Computed_Nu_P_0  ( nPointsE, nPointsE ), &
        Computed_Nu_P_1  ( nPointsE, nPointsE ), &
        Computed_NuB_P_0 ( nPointsE, nPointsE ), &
        Computed_NuB_P_1 ( nPointsE, nPointsE ) )

    CALL Interpolate_NNS_Point &
           ( T, MuB, &
             Interpolated_Nu_N_0,  Interpolated_Nu_N_1,  &
             Interpolated_NuB_N_0, Interpolated_NuB_N_1, &
             Interpolated_Nu_P_0,  Interpolated_Nu_P_1,  &
             Interpolated_NuB_P_0, Interpolated_NuB_P_1 )

    end associate !-- T_Min, etc.
    end associate !-- iT 

    WRITE (*,*)

  END BLOCK TestPoint

  
CONTAINS


  SUBROUTINE Interpolate_NNS_Point &
               ( T, MuB, &
                 Interpolated_Nu_N_0,  Interpolated_Nu_N_1,  &
                 Interpolated_NuB_N_0, Interpolated_NuB_N_1, &
                 Interpolated_Nu_P_0,  Interpolated_Nu_P_1,  &
                 Interpolated_NuB_P_0, Interpolated_NuB_P_1 )

    USE wlKindModule, ONLY: &
      dp
    USE wlInterpolationModule, ONLY: &
      GetIndexAndDelta_Log, &
      GetIndexAndDelta_Lin, &
      LinearInterp2D_2DArray_Point

    REAL(dp), INTENT(in)  :: T, MuB
    REAL(dp), DIMENSION ( :, : ), INTENT(out), TARGET :: &
      Interpolated_Nu_N_0,  Interpolated_Nu_N_1,  &
      Interpolated_NuB_N_0, Interpolated_NuB_N_1, &
      Interpolated_Nu_P_0,  Interpolated_Nu_P_1,  &
      Interpolated_NuB_P_0, Interpolated_NuB_P_1

    INTEGER  :: iOpacity, iMoment, iE, iEp, iT, iMuB
    REAL(dp) :: dLogT, dMuB
    REAL(dp), DIMENSION ( :, : ), POINTER :: Interpolated

    associate &
      ( iT_TS  =>  OpacityTable % TS % Indices % iT )
    associate &
      ( T_Values    =>  OpacityTable % TS % States ( iT_TS ) % Values, &
        MuB_Values  =>  OpacityTable % MuBGrid % Values )
  
    !-- Get Index and Delta

    CALL GetIndexAndDelta_Log ( T,     T_Values, iT,   dLogT )
    CALL GetIndexAndDelta_Lin ( MuB, MuB_Values, iMuB, dMuB )

    !-- Check Index and Delta values

    WRITE (*,*)
    WRITE (*,'(A15,I4.4)')     '          iT = ', iT 
    WRITE (*,'(A15,ES12.6E2)') 'T ( iT )     = ', T_Values ( iT )
    WRITE (*,'(A15,ES12.6E2)') 'T ( iT + 1 ) = ', T_Values ( iT + 1 )
    WRITE (*,'(A15,ES12.6E2)') '       dLogT = ', dLogT
    WRITE (*,'(A15,ES12.6E2)') 'T            = ', &
          10.d0 ** ( LOG10 ( T_Values ( iT ) )  &
                     +  LOG10 ( T_Values ( iT + 1 ) / T_Values ( iT ) )  &
                        * dLogT )

    WRITE (*,*)
    WRITE (*,'(A19,I4.4)')     '            iMuB = ', iMuB 
    WRITE (*,'(A19,ES12.6E2)') 'MuB ( iMuB )     = ', MuB_Values ( iMuB )
    WRITE (*,'(A19,ES12.6E2)') 'MuB ( iMuB + 1 ) = ', MuB_Values ( iMuB + 1 )
    WRITE (*,'(A19,ES12.6E2)') '            dMuB = ', dMuB
    WRITE (*,'(A19,ES12.6E2)') 'MuB              = ', &
          MuB_Values ( iMuB )  &
          +  ( MuB_Values ( iMuB + 1 )  -  MuB_Values ( iMuB ) )  &
             *  dMuB 

    associate & 
      ( Table_NNS  =>  OpacityTable % Scat_NNS )

    do iOpacity  =  1,  Table_NNS % nOpacities
      do iMoment  =  1,  1 !Table_NNS % nMoments
 
        select case ( iOpacity )
        case ( 1 )
          select case ( iMoment )
            case ( 1 )
              Interpolated  =>  Interpolated_Nu_N_0
            case ( 2 )
              Interpolated  =>  Interpolated_Nu_N_1
          end select !-- iMoment
        case ( 2 )
          select case ( iMoment )
            case ( 1 )
              Interpolated  =>  Interpolated_NuB_N_0
            case ( 2 )
              Interpolated  =>  Interpolated_NuB_N_1
          end select !-- iMoment
        case ( 3 )
          select case ( iMoment )
            case ( 1 )
              Interpolated  =>  Interpolated_Nu_P_0
            case ( 2 )
              Interpolated  =>  Interpolated_Nu_P_1
          end select !-- iMoment
        case ( 4 )
          select case ( iMoment )
            case ( 1 )
              Interpolated  =>  Interpolated_NuB_P_0
            case ( 2 )
              Interpolated  =>  Interpolated_NuB_P_1
          end select !-- iMoment
        end select !-- iOpacity

        iEp = 15
        iE  = 25

        associate &
          ( Table       =>  OpacityTable % Scat_NNS % Kernel ( iOpacity ) &
                              % Values ( iEp, iE, iMoment, :, : ), &
            Offset      =>  OpacityTable % Scat_NNS &
                              % Offsets ( iOpacity, iMoment ) )

        !-- Interpolate

        CALL LinearInterp2D_2DArray_Point &
               ( iT, iMuB, dLogT, dMuB, Offset, Table, &
                 Interpolated ( iEp, iE ) )

        !-- Check Interpolated

        WRITE (*,*)
        WRITE (*,'(A15,ES12.6E2)') 'Table 00    = ', &
          10.d0 ** ( Table ( iT,   iMuB   ) )  -  Offset    
        WRITE (*,'(A15,ES12.6E2)') 'Table 10    = ', &
          10.d0 ** ( Table ( iT+1, iMuB   ) )  -  Offset    
        WRITE (*,'(A15,ES12.6E2)') 'Table 01    = ', &
          10.d0 ** ( Table ( iT,   iMuB+1 ) )  -  Offset    
        WRITE (*,'(A15,ES12.6E2)') 'Table 11    = ', &
          10.d0 ** ( Table ( iT+1, iMuB+1 ) )  -  Offset    
        WRITE (*,'(A15,ES12.6E2)') 'Interpolated = ', Interpolated ( iEp, iE )    
        end associate !-- Table, etc.

      end do !-- iMoment
    end do !-- iOpacity

    end associate !-- Table_NNS
    end associate !-- T_Values, etc.
    end associate !-- iT_TS

  END SUBROUTINE Interpolate_NNS_Point


  SUBROUTINE Compute_NNS_Point &
               ( T, MuB, &
                 Computed_Nu_N_0,  Computed_Nu_N_1,  &
                 Computed_NuB_N_0, Computed_NuB_N_1, &
                 Computed_Nu_P_0,  Computed_Nu_P_1,  &
                 Computed_NuB_P_0, Computed_NuB_P_1 )

    USE wlKindModule, ONLY: &
      dp

    REAL(dp), INTENT(in)  :: T, MuB
    REAL(dp), DIMENSION ( :, : ), INTENT(out) :: &
      Computed_Nu_N_0,  Computed_Nu_N_1,  &
      Computed_NuB_N_0, Computed_NuB_N_1, &
      Computed_Nu_P_0,  Computed_Nu_P_1,  &
      Computed_NuB_P_0, Computed_NuB_P_1

  END SUBROUTINE Compute_NNS_Point


END PROGRAM wlReadOpacityTableTest
