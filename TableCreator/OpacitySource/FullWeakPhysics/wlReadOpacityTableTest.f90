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


!   !-- Test random point

!   TestRandom: BLOCK

!     USE wlKindModule, ONLY: &
!       dp

!     REAL(dp) :: random
!     REAL(dp) :: T, MuB

!     WRITE (*,*)
!     WRITE (*,*) '>>> Random ( T, MuB ) point'

!     CALL random_seed ( )

!     associate &
!       ( iT_TS     =>  OpacityTable % TS % Indices % iT )
!     associate &
!       ( T_Min     =>  OpacityTable % TS % minValues ( iT_TS ), &
!         T_Max     =>  OpacityTable % TS % maxValues ( iT_TS ), &
!         MuB_Min   =>  OpacityTable % MuBGrid % minValue, &
!         MuB_Max   =>  OpacityTable % MuBGrid % maxValue )

!     CALL random_number ( random )
!     T  =  10.d0 ** ( LOG10 ( T_Min )  &
!                      +  random  *  ( LOG10 ( T_Max )  -  LOG10 ( T_Min ) ) ) 

!     CALL random_number ( random )
!     MuB  =  MuB_Min  +  random * ( MuB_Max - MuB_Min )

!     !-- Temporary example values
! !    T   = 1.402089E+12
! !    MuB = 1.034229E+03
! !    call Test_NNS_Point ( T, MuB, MuB )
  
!     end associate !-- T_Min, etc.
!     end associate !-- iT 

!   END BLOCK TestRandom


  !-- Bruenn et al. (2020) figures

  BruennFigures: BLOCK

    USE wlKindModule, ONLY: &
      dp

    INTEGER  :: j_Rho, k_T, l_Ye 
    REAL(dp) :: Rho, T, Ye, MuN, MuP, S_tot 

!    j_Rho = 118 !-- 10^11 g cm^-3 
!    k_t = 25  !--  0.9 MeV
!    l_ye = 20  !-- 0.4

    WRITE (*,*)
    WRITE (*,*) '>>> Condition 1, Figs. 30-31, Bruenn et al. (2020)'

    j_Rho = 103  !-- 10^10 g cm^-3 
    k_T   =  38  !--  3.0 MeV
    l_Ye  =  13  !-- 0.25
    call SetFluid ( j_Rho, k_T, l_Ye, Rho, T, Ye, MuN, MuP, S_tot )
    call Test_NNS_Point ( T, MuN, MuP, S_tot )

    WRITE (*,*)
    WRITE (*,*) '>>> Condition 2, Figs. 30-31, Bruenn et al. (2020)'

    j_Rho = 133  !-- 10^12 g cm^-3 
    k_T   =  51  !-- 10.0 MeV
    l_Ye  =   5  !-- 0.1
    call SetFluid ( j_Rho, k_T, l_Ye, Rho, T, Ye, MuN, MuP, S_tot )
    call Test_NNS_Point ( T, MuN, MuP, S_tot )

    WRITE (*,*)
    WRITE (*,*) '>>> Condition 3, Figs. 30-31, Bruenn et al. (2020)'

    j_Rho = 163  !-- 10^14 g cm^-3 
    k_T   =  53  !-- 12.0 MeV
    l_Ye  =  14  !-- 0.27
    call SetFluid ( j_Rho, k_T, l_Ye, Rho, T, Ye, MuN, MuP, S_tot )
    call Test_NNS_Point ( T, MuN, MuP, S_tot )

  END BLOCK BruennFigures


CONTAINS


  SUBROUTINE SetFluid ( j_Rho, k_T, l_Ye, Rho, T, Ye, MuN, MuP, S_tot )

    USE wlKindModule, ONLY: &
      dp
    USE wlExtPhysicalConstantsModule, ONLY: &
      kMeV, kfm, dmnp, mn_wl => mn, mp_wl => mp
    USE wlExtNumericalModule, ONLY: epsilon

    INTEGER,  INTENT(in)  :: j_Rho, k_T, l_Ye
    REAL(dp), INTENT(out) :: Rho, T, Ye, MuN, MuP, S_tot

    REAL(dp) :: chem_n, chem_p

    ASSOCIATE (  &
       iRho    => OpacityTable % TS % Indices % iRho, &
       iT      => OpacityTable % TS % Indices % iT, &
       iYe     => OpacityTable % TS % Indices % iYe, &
       Indices => OpacityTable % EOSTable % DV % Indices, &
       DVOffs  => OpacityTable % EOSTable % DV % Offsets, &
       DVar    => OpacityTable % EOSTable % DV % Variables )
       ! nRho    => OpacityTable % nPointsTS(OpacityTable % TS % Indices % iRho) , &
       ! nT      => OpacityTable % nPointsTS(OpacityTable % TS % Indices % iT)   , &
       ! nYe     => OpacityTable % nPointsTS(OpacityTable % TS % Indices % iYe)  , &
 
    WRITE (*,*)

  ! WRITE (*,*)
  ! WRITE (*,*) '>>> Rho values'
  ! DO j_rho = 1, nRho
  !   rho = OpacityTable % TS % States (iRho) % Values (j_rho)
  !   WRITE (*,'(A7,I3.3,A4,ES10.3E2)') 'Rho    ', j_rho, '    ', rho
  ! END DO

  ! WRITE (*,*)
  ! WRITE (*,*) '>>> T values'
  ! DO k_t = 1, nT
  !   T = OpacityTable % TS % States (iT) % Values (k_t)
  !   TMeV = T * kMeV
  !   WRITE (*,'(A7,I3.3,A4,ES10.3E2)') 'T_MeV  ', k_t, '    ', TMeV
  ! END DO

  ! WRITE (*,*)
  ! WRITE (*,*) '>>> Ye values'
  ! DO l_ye = 1, nYe
  !   ye = OpacityTable % TS % States (iYe) % Values (l_ye)
  !   WRITE (*,'(A7,I3.3,A4,ES10.3E2)') 'Ye     ', l_ye, '    ', ye
  ! END DO

    Rho = OpacityTable % TS % States (iRho) % Values (j_Rho)
    WRITE (*,'(A14,I3.3,A4,ES10.3E2)') 'Rho (g cm^-3) ', j_Rho, '    ', Rho

    T = OpacityTable % TS % States (iT) % Values (k_t)
    WRITE (*,'(A14,I3.3,A4,ES10.3E2)') '  T (K)       ', k_t, '    ', T
    WRITE (*,'(A14,I3.3,A4,ES10.3E2)') '  T (MeV)     ', k_t, '    ', T * kMeV

    Ye = OpacityTable % TS % States (iYe) % Values (l_ye)
    WRITE (*,'(A14,I3.3,A4,ES10.3E2)') ' Ye           ', l_ye, '    ', Ye

    chem_n = 10**DVar(Indices % iNeutronChemicalPotential) % &
                   Values(j_Rho, k_T, l_Ye) &
             - DVOffs(Indices % iNeutronChemicalPotential)   &
             - epsilon
    MuN  =  chem_n + dmnp + mn_wl  !-- Convert from Chimera to absolute
    WRITE (*,'(A10,ES10.3E2)') 'MuN (MeV) ', MuN

    chem_p = 10**DVar(Indices % iProtonChemicalPotential) % &
                   Values(j_Rho, k_T, l_Ye) &
             - DVOffs(Indices % iProtonChemicalPotential)   &
             - epsilon
    MuP  =  chem_p + dmnp + mp_wl  !-- Convert from Chimera to absolute
    WRITE (*,'(A10,ES10.3E2)') 'MuP (MeV) ', MuP

    CALL nc_manybody_corrections_weaklib( Rho * kfm, T * kMeV, Ye, S_tot )
    WRITE (*,'(A10,ES10.3E2)') 'S_tot     ', S_tot

    END ASSOCIATE ! rho-T-Ye

  END SUBROUTINE SetFluid


  SUBROUTINE Test_NNS_Point ( T, MuN, MuP, S_tot )

    USE wlKindModule, ONLY: &
      dp
    USE wlExtPhysicalConstantsModule, ONLY: &
      kMeV
    USE wlOpacityFieldsModule, ONLY: &
      iNu_NNS, iNuBar_NNS, &
      iNeutron_NNS, iProton_NNS

    REAL(dp), INTENT(in) :: T, MuN, MuP, S_tot

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
      ( nPointsE  =>  OpacityTable % EnergyGrid % nPoints, &
        E_Edges => OpacityTable % EnergyGrid % Edge ) 

    ! WRITE (*,*)
    ! WRITE (*,'(A6,ES13.6E2)') 'T   = ', T 
    ! WRITE (*,'(A6,ES13.6E2)') 'MuN = ', MuN 
    ! WRITE (*,'(A6,ES13.6E2)') 'MuP = ', MuP 
    ! WRITE (*,'(A11,I4.4)')    'nPointsE = ', nPointsE 

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
           ( T, MuN, MuP, nPointsE, &
             Interpolated_Nu_N_0,  Interpolated_Nu_N_1,  &
             Interpolated_NuB_N_0, Interpolated_NuB_N_1, &
             Interpolated_Nu_P_0,  Interpolated_Nu_P_1,  &
             Interpolated_NuB_P_0, Interpolated_NuB_P_1 )
    CALL Compute_NNS_Point &
           ( T, MuN, MuP, nPointsE, &
             Computed_Nu_N_0,  Computed_Nu_N_1,  &
             Computed_NuB_N_0, Computed_NuB_N_1, &
             Computed_Nu_P_0,  Computed_Nu_P_1,  &
             Computed_NuB_P_0, Computed_NuB_P_1 )

    call TestMeanFreePath &
           ( Interpolated_Nu_N_0, Interpolated_NuB_N_0, &
             Interpolated_Nu_P_0, Interpolated_NuB_P_0, &
             Computed_Nu_N_0, Computed_NuB_N_0, &
             Computed_Nu_P_0, Computed_NuB_P_0, &
             E_Edges, S_tot, nPointsE )

    end associate !-- nPointsE

  END SUBROUTINE Test_NNS_Point


  SUBROUTINE Interpolate_NNS_Point &
               ( T, MuN, MuP, nPointsE, &
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

    REAL(dp), INTENT(in)  :: T, MuN, MuP
    INTEGER,  INTENT(in)  :: nPointsE
    REAL(dp), DIMENSION ( :, : ), INTENT(out), TARGET :: &
      Interpolated_Nu_N_0,  Interpolated_Nu_N_1,  &
      Interpolated_NuB_N_0, Interpolated_NuB_N_1, &
      Interpolated_Nu_P_0,  Interpolated_Nu_P_1,  &
      Interpolated_NuB_P_0, Interpolated_NuB_P_1

    INTEGER  :: iOpacity, iMoment, iE, iEp, iT, iMuN, iMuP, iMuB
    REAL(dp) :: dLogT, dMuN, dMuP, dMuB
    REAL(dp), DIMENSION ( :, : ), POINTER :: Interpolated

    associate &
      ( iT_TS  =>  OpacityTable % TS % Indices % iT )
    associate &
      ( T_Values    =>  OpacityTable % TS % States ( iT_TS ) % Values, &
        MuB_Values  =>  OpacityTable % MuBGrid % Values )
  
    !-- Get Index and Delta

    CALL GetIndexAndDelta_Log ( T,     T_Values, iT,   dLogT )
    CALL GetIndexAndDelta_Lin ( MuN, MuB_Values, iMuN, dMuN )
    CALL GetIndexAndDelta_Lin ( MuP, MuB_Values, iMuP, dMuP )

    ! !-- Check Index and Delta values

    ! WRITE (*,*)
    ! WRITE (*,'(A15,I4.4)')     '          iT = ', iT 
    ! WRITE (*,'(A15,ES13.6E2)') 'T ( iT )     = ', T_Values ( iT )
    ! WRITE (*,'(A15,ES13.6E2)') 'T ( iT + 1 ) = ', T_Values ( iT + 1 )
    ! WRITE (*,'(A15,ES13.6E2)') '       dLogT = ', dLogT
    ! WRITE (*,'(A15,ES13.6E2)') 'T            = ', &
    !       10.d0 ** ( LOG10 ( T_Values ( iT ) )  &
    !                  +  LOG10 ( T_Values ( iT + 1 ) / T_Values ( iT ) )  &
    !                     * dLogT )

    ! WRITE (*,*)
    ! WRITE (*,'(A19,I4.4)')     '            iMuN = ', iMuN 
    ! WRITE (*,'(A19,ES13.6E2)') 'MuB ( iMuN )     = ', MuB_Values ( iMuN )
    ! WRITE (*,'(A19,ES13.6E2)') 'MuB ( iMuN + 1 ) = ', MuB_Values ( iMuN + 1 )
    ! WRITE (*,'(A19,ES13.6E2)') '            dMuN = ', dMuN
    ! WRITE (*,'(A19,ES13.6E2)') 'MuN              = ', &
    !       MuB_Values ( iMuN )  &
    !       +  ( MuB_Values ( iMuN + 1 )  -  MuB_Values ( iMuN ) )  &
    !          *  dMuN 

    ! WRITE (*,*)
    ! WRITE (*,'(A19,I4.4)')     '            iMuP = ', iMuP 
    ! WRITE (*,'(A19,ES13.6E2)') 'MuB ( iMuP )     = ', MuB_Values ( iMuP )
    ! WRITE (*,'(A19,ES13.6E2)') 'MuB ( iMuP + 1 ) = ', MuB_Values ( iMuP + 1 )
    ! WRITE (*,'(A19,ES13.6E2)') '            dMuP = ', dMuP
    ! WRITE (*,'(A19,ES13.6E2)') 'MuP              = ', &
    !       MuB_Values ( iMuP )  &
    !       +  ( MuB_Values ( iMuP + 1 )  -  MuB_Values ( iMuP ) )  &
    !          *  dMuP 

    associate & 
      ( Table_NNS  =>  OpacityTable % Scat_NNS )

    do iOpacity  =  1,  Table_NNS % nOpacities
      do iMoment  =  1,  1 !Table_NNS % nMoments
 
        select case ( iOpacity )
        case ( 1 )
          iMuB  =  iMuN
          dMuB  =  dMuN
          select case ( iMoment )
            case ( 1 )
              Interpolated  =>  Interpolated_Nu_N_0
            case ( 2 )
              Interpolated  =>  Interpolated_Nu_N_1
          end select !-- iMoment
        case ( 2 )
          iMuB  =  iMuN
          dMuB  =  dMuN
          select case ( iMoment )
            case ( 1 )
              Interpolated  =>  Interpolated_NuB_N_0
            case ( 2 )
              Interpolated  =>  Interpolated_NuB_N_1
          end select !-- iMoment
        case ( 3 )
          iMuB  =  iMuP
          dMuB  =  dMuP
          select case ( iMoment )
            case ( 1 )
              Interpolated  =>  Interpolated_Nu_P_0
            case ( 2 )
              Interpolated  =>  Interpolated_Nu_P_1
          end select !-- iMoment
        case ( 4 )
          iMuB  =  iMuP
          dMuB  =  dMuP
          select case ( iMoment )
            case ( 1 )
              Interpolated  =>  Interpolated_NuB_P_0
            case ( 2 )
              Interpolated  =>  Interpolated_NuB_P_1
          end select !-- iMoment
        end select !-- iOpacity

        do iE = 1, nPointsE
          do iEp = 1, nPointsE

            associate &
              ( Table   =>  OpacityTable % Scat_NNS % Kernel ( iOpacity ) &
                              % Values ( iEp, iE, iMoment, :, : ), &
                Offset  =>  OpacityTable % Scat_NNS &
                              % Offsets ( iOpacity, iMoment ) )

            !-- Interpolate

            CALL LinearInterp2D_2DArray_Point &
                   ( iT, iMuB, dLogT, dMuB, Offset, Table, &
                     Interpolated ( iEp, iE ) )

            end associate !-- Table, etc.

          end do !-- iEp
        end do !-- iE

        ! !-- Spot check

        ! iEp = 5
        ! iE  = 7

        ! associate &
        !   ( Table   =>  OpacityTable % Scat_NNS % Kernel ( iOpacity ) &
        !                   % Values ( iEp, iE, iMoment, :, : ), &
        !     Offset  =>  OpacityTable % Scat_NNS &
        !                   % Offsets ( iOpacity, iMoment ) )

        ! WRITE (*,*)
        ! WRITE (*,'(A15,ES13.6E2)') 'Table 00    = ', &
        !   10.d0 ** ( Table ( iT,   iMuB   ) )  -  Offset    
        ! WRITE (*,'(A15,ES13.6E2)') 'Table 10    = ', &
        !   10.d0 ** ( Table ( iT+1, iMuB   ) )  -  Offset    
        ! WRITE (*,'(A15,ES13.6E2)') 'Table 01    = ', &
        !   10.d0 ** ( Table ( iT,   iMuB+1 ) )  -  Offset    
        ! WRITE (*,'(A15,ES13.6E2)') 'Table 11    = ', &
        !   10.d0 ** ( Table ( iT+1, iMuB+1 ) )  -  Offset    
        ! WRITE (*,'(A15,ES13.6E2)') 'Interpolated = ', Interpolated ( iEp, iE )    
        ! end associate !-- Table, etc.

      end do !-- iMoment
    end do !-- iOpacity

    end associate !-- Table_NNS
    end associate !-- T_Values, etc.
    end associate !-- iT_TS

  END SUBROUTINE Interpolate_NNS_Point


  SUBROUTINE Compute_NNS_Point &
               ( T, MuN, MuP, nPointsE, &
                 Computed_Nu_N_0,  Computed_Nu_N_1,  &
                 Computed_NuB_N_0, Computed_NuB_N_1, &
                 Computed_Nu_P_0,  Computed_Nu_P_1,  &
                 Computed_NuB_P_0, Computed_NuB_P_1 )

    USE wlKindModule, ONLY: &
      dp
    USE wlExtPhysicalConstantsModule, ONLY: &
      kMeV, dmnp, mn_wl => mn, mp_wl => mp
    USE scat_n_module_weaklib, ONLY: &
      init_quad_scat_n

    REAL(dp), INTENT(in)  :: T, MuN, MuP
    INTEGER,  INTENT(in)  :: nPointsE
    REAL(dp), DIMENSION ( :, : ), INTENT(out) :: &
      Computed_Nu_N_0,  Computed_Nu_N_1,  &
      Computed_NuB_N_0, Computed_NuB_N_1, &
      Computed_Nu_P_0,  Computed_Nu_P_1,  &
      Computed_NuB_P_0, Computed_NuB_P_1

    INTEGER, PARAMETER :: Scat_weak_magnetism &
                          = 1
    REAL(DP), PARAMETER :: Scat_ga_strange &
!                               = -0.1d0
                           = 0.0d0
    REAL(dp) :: TMev, chem_n, chem_p
    REAL(dp), DIMENSION(nPointsE, nPointsE) :: phi0_nu_n,  phi1_nu_n, &
                                               phi0_nub_n, phi1_nub_n, &
                                               phi0_nu_p,  phi1_nu_p, &
                                               phi0_nub_p, phi1_nub_p

    CALL init_quad_scat_n
    CALL load_polylog_weaklib

    !-- convert to Chimera conventions

    TMeV    =  T * kMeV
    chem_n  =  MuN - dmnp - mn_wl
    chem_p  =  MuP - dmnp - mp_wl

    ! WRITE (*,*)
    ! WRITE (*,'(A9,ES13.6E2)') 'TMeV   = ', TMeV
    ! WRITE (*,'(A9,ES13.6E2)') 'chem_n = ', chem_n
    ! WRITE (*,'(A9,ES13.6E2)') 'chem_p = ', chem_p

    CALL scatnrgn_weaklib &
         ( nPointsE, &
           OpacityTable % EnergyGrid % Values, &
           OpacityTable % EnergyGrid % Edge, &
           TMeV, chem_n, chem_p, Scat_weak_magnetism, Scat_ga_strange, &
           phi0_nu_n, phi1_nu_n, phi0_nub_n, phi1_nub_n, &
           phi0_nu_p, phi1_nu_p, phi0_nub_p, phi1_nub_p )

    Computed_Nu_N_0  =  0.5_DP * TRANSPOSE(phi0_nu_n(:,:))  
            ! phi0_nu_n was saved as phi0_nu_n(e,ep)

    Computed_Nu_N_1  =  1.5_DP * TRANSPOSE(phi1_nu_n(:,:))  
            ! phi1_nu_n was saved as phi1_nu_n(e,ep)

    Computed_NuB_N_0  =  0.5_DP * TRANSPOSE(phi0_nub_n(:,:))  
            ! phi0_nub_n was saved as phi0_nub_n(e,ep)

    Computed_NuB_N_1  =  1.5_DP * TRANSPOSE(phi1_nub_n(:,:))  
            ! phi1_nub_n was saved as phi1_nub_n(e,ep)

    Computed_Nu_P_0  =  0.5_DP * TRANSPOSE(phi0_nu_p(:,:))  
            ! phi0_nu_p was saved as phi0_nu_p(e,ep)

    Computed_Nu_P_1  =  1.5_DP * TRANSPOSE(phi1_nu_p(:,:))  
            ! phi1_nu_p was saved as phi1_nu_p(e,ep)

    Computed_NuB_P_0  =  0.5_DP * TRANSPOSE(phi0_nub_p(:,:))  
            ! phi0_nub_p was saved as phi0_nub_p(e,ep)

    Computed_NuB_P_1  =  1.5_DP * TRANSPOSE(phi1_nub_p(:,:))  
            ! phi1_nub_p was saved as phi1_nub_p(e,ep)

    ! WRITE (*,*)
    ! WRITE (*,'(A15,ES13.6E2)') 'Computed = ', Computed_Nu_N_0  ( 5, 7 )  
    ! WRITE (*,'(A15,ES13.6E2)') 'Computed = ', Computed_NuB_N_0 ( 5, 7 )  
    ! WRITE (*,'(A15,ES13.6E2)') 'Computed = ', Computed_Nu_P_0  ( 5, 7 )  
    ! WRITE (*,'(A15,ES13.6E2)') 'Computed = ', Computed_NuB_P_0 ( 5, 7 )  

  END SUBROUTINE Compute_NNS_Point


  SUBROUTINE TestMeanFreePath &
               ( Interpolated_Nu_N_0, Interpolated_NuB_N_0, &
                 Interpolated_Nu_P_0, Interpolated_NuB_P_0, &
                 Computed_Nu_N_0, Computed_NuB_N_0, &
                 Computed_Nu_P_0, Computed_NuB_P_0, &
                 E_Edges, S_tot, nPointsE )

    USE wlKindModule, ONLY: &
      dp

    REAL(dp), DIMENSION ( :, : ), INTENT(in) :: &
      Interpolated_Nu_N_0, Interpolated_NuB_N_0, &
      Interpolated_Nu_P_0, Interpolated_NuB_P_0, &
      Computed_Nu_N_0, Computed_NuB_N_0, &
      Computed_Nu_P_0, Computed_NuB_P_0
    REAL(dp), DIMENSION ( : ), INTENT(in) :: &
      E_Edges
    REAL(dp), INTENT(in) :: S_tot
    INTEGER, INTENT(in) :: nPointsE

    INTEGER  :: iEp, iE
    REAL(dp) :: TwoPi
    REAL(dp), DIMENSION( nPointsE ) :: E_Vol, &
                                       invMFP_Nu_I,  invMFP_Nu_C, &
                                       invMFP_NuB_I, invMFP_NuB_C
    
    TwoPi  =  2.0d0 * acos ( -1.0d0 )

    DO iE = 1, nPointsE
      E_Vol ( iE )  =  ( E_Edges ( iE + 1 ) ** 3  -  E_Edges ( iE ) ** 3 ) &
                       / 3.d0 
    END DO !-- iE

    invMFP_Nu_I   =  0.d0
    invMFP_Nu_C   =  0.d0
    invMFP_NuB_I  =  0.d0
    invMFP_NuB_C  =  0.d0
    DO iE = 1, nPointsE
      DO iEp = 1, nPointsE
        invMFP_Nu_I ( iE )  &
          =  invMFP_Nu_I ( iE )  &
             +  E_Vol ( iEp )  *  (    Interpolated_Nu_N_0 ( iEp, iE )  &
                                    +  Interpolated_Nu_P_0 ( iEp, iE ) )  
        invMFP_Nu_C ( iE )  &
          =  invMFP_Nu_C ( iE )  &
             +  E_Vol ( iEp )  *  (    Computed_Nu_N_0 ( iEp, iE )  &
                                    +  Computed_Nu_P_0 ( iEp, iE ) )  
        invMFP_NuB_I ( iE )  &
          =  invMFP_NuB_I ( iE )  &
             +  E_Vol ( iEp )  *  (    Interpolated_NuB_N_0 ( iEp, iE )  &
                                    +  Interpolated_NuB_P_0 ( iEp, iE ) )  
        invMFP_NuB_C ( iE )  &
          =  invMFP_NuB_C ( iE )  &
             +  E_Vol ( iEp )  *  (    Computed_NuB_N_0 ( iEp, iE )  &
                                    +  Computed_NuB_P_0 ( iEp, iE ) )  
      END DO
      !-- Factor of 2 to undo thornado legendre moment convention
      !-- Factor of TwoPi for azimuthal integral, Bruenn et al. (2020) Eq. (364)
      !-- Many-body corrections S_tot
      invMFP_Nu_I  ( iE )  =  2.d0 * TwoPi * invMFP_Nu_I  ( iE ) * S_tot
      invMFP_Nu_C  ( iE )  =  2.d0 * TwoPi * invMFP_Nu_C  ( iE ) * S_tot
      invMFP_NuB_I ( iE )  =  2.d0 * TwoPi * invMFP_NuB_I ( iE ) * S_tot
      invMFP_NuB_C ( iE )  =  2.d0 * TwoPi * invMFP_NuB_C ( iE ) * S_tot
    END DO

    WRITE (*,*)
    WRITE (*,*) 'Nu Interpolated, Computed'
    DO iE = 1, nPointsE
      WRITE (*,'(ES13.6E2,ES13.6E2)') invMFP_Nu_I ( iE ), invMFP_Nu_C ( iE )    
    END DO

    WRITE (*,*)
    WRITE (*,*) 'NuB Interpolated, Computed'
    DO iE = 1, nPointsE
      WRITE (*,'(ES13.6E2,ES13.6E2)') invMFP_NuB_I ( iE ), invMFP_NuB_C ( iE )
    END DO

  END SUBROUTINE TestMeanFreePath


END PROGRAM wlReadOpacityTableTest
