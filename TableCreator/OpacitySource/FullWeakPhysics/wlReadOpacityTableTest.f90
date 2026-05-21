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

  TYPE(OpacityTableType), target :: OpacityTable

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

    INTEGER  :: iEp, iE, iMom, iNu, iTarget
    REAL(dp) :: random
    REAL(dp) :: T, MuB
    REAL(dp) :: Interpolant

    iMom    = 1
    iNu     = iNu_NNS
    iTarget = iNeutron_NNS

    iEp  = 15
    iE   = 25

    associate &
      ( iT  =>  OpacityTable % TS % Indices % iT )
    associate &
      ( T_Min    =>  OpacityTable % TS % minValues ( iT ), &
        T_Max    =>  OpacityTable % TS % maxValues ( iT ), &
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

    CALL Interpolate_NNS_Point &
           ( T, MuB, iEp, iE, iMom, iNu, iTarget, Interpolant )

    end associate !-- T_Min, etc.
    end associate !-- iT 

    WRITE (*,*)

  END BLOCK TestPoint

  
CONTAINS


  SUBROUTINE Interpolate_NNS_Point &
               ( T, MuB, iEp, iE, iMom, iNu, iTarget, Interpolant )

    USE wlKindModule, ONLY: &
      dp
    USE wlInterpolationModule, ONLY: &
      GetIndexAndDelta_Log, &
      GetIndexAndDelta_Lin, &
      LinearInterp2D_2DArray_Point

    REAL(dp), INTENT(in)  :: T, MuB
    INTEGER,  INTENT(in)  :: iEp, iE, iMom, iNu, iTarget
    REAL(dp), INTENT(out) :: Interpolant

    INTEGER  :: iT, iMuB
    REAL(dp) :: dLogT, dMuB
    REAL(dp), dimension ( :, : ), pointer :: Offset_2D

    associate &
      ( iT  =>  OpacityTable % TS % Indices % iT )
    associate &
      ( Table       =>  OpacityTable % Scat_NNS % Phi ( iNu, iTarget ) &
                          % Values ( iEp, iE, iMom, :, : ), &
        T_Values    =>  OpacityTable % TS % States ( iT ) % Values, &
        MuB_Values  =>  OpacityTable % MuBGrid % Values )

    Offset_2D ( 1:2, 1:2 )  &
      =>  OpacityTable % Scat_NNS % Offsets( 1:4 , iMom )

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

    !-- Interpolate

    CALL LinearInterp2D_2DArray_Point &
           ( iT, iMuB, dLogT, dMuB, Offset_2D ( iNu, iTarget ), Table, &
             Interpolant )

    !-- Check Interpolant

    WRITE (*,*)
    WRITE (*,'(A14,ES12.6E2)') 'Table 00    = ', &
      10.d0 ** ( Table ( iT,   iMuB   ) )  -  Offset_2D ( iNu, iTarget )    
    WRITE (*,'(A14,ES12.6E2)') 'Table 10    = ', &
      10.d0 ** ( Table ( iT+1, iMuB   ) )  -  Offset_2D ( iNu, iTarget )    
    WRITE (*,'(A14,ES12.6E2)') 'Table 01    = ', &
      10.d0 ** ( Table ( iT,   iMuB+1 ) )  -  Offset_2D ( iNu, iTarget )    
    WRITE (*,'(A14,ES12.6E2)') 'Table 11    = ', &
      10.d0 ** ( Table ( iT+1, iMuB+1 ) )  -  Offset_2D ( iNu, iTarget )    
    WRITE (*,'(A14,ES12.6E2)') 'Interpolant = ', Interpolant    

    end associate !-- Table, etc.
    end associate !-- iT

  END SUBROUTINE Interpolate_NNS_Point


END PROGRAM wlReadOpacityTableTest
