PROGRAM wlTestGeneralKernelTwoReactions

  USE wlKindModule,         ONLY: dp
  USE wlEosConstantsModule, ONLY: &
    mp, mn, me, mmu, kmev, kmev_inv, ca, cv, pi
  USE wlIOModuleHDF, ONLY: &
    InitializeHDF, &
    FinalizeHDF
  USE wlHelmIOModuleHDF, ONLY: &
    ReadHelmholtzTableHDF
  USE wlLeptonEOSTableModule, ONLY: &
    HelmTableType
  USE wlLeptonPhotonGasEOS, ONLY: &
    LeptonGasType, LeptonGasEOS
  USE wlOpacityTableIOModuleHDF, ONLY: &
    ReadOpacityTableHDF
  USE wlGridModule, ONLY: &
    GridType, &
    AllocateGrid, &
    DescribeGrid, &
    MakeLogGrid
  USE wlGeneralLeptonScatteringModule, ONLY: &
    ProcessIndexFromReactionString
  USE wlGeneralLeptonScatteringModuleThornadoInterface, ONLY: &
    CalculateAllPhout, &
    CalculateAllRoutIntegrated, &
    InitGeneralScatteringKernels, &
    FinalizeGeneralScatteringKernels, &
    iProcessMin_Default, &
    iProcessMax_Default
  USE HDF5

  IMPLICIT NONE

  !--------- parameters for creating energy grid ----------------------------
  INTEGER, PARAMETER     :: Inte_nPointE = 80
  INTEGER, PARAMETER     :: nE1 = Inte_nPointE
  INTEGER, PARAMETER     :: nE3 = nE1 - 1
  REAL(dp)               :: Inte_Emin = 1.0d-1
  REAL(dp)               :: Inte_Emax = 3.0d02
  REAL(dp)               :: E1_arr(nE1), E3_arr(nE3)
  TYPE(GridType)         :: Inte_E

  !--- Variables for profile read ---
  REAL(dp), ALLOCATABLE :: Inte_T(:) , Inte_Rho(:), Inte_Ym(:), &
                           Inte_Ye(:), Inte_Mue(:), Inte_Mum(:), &
                           Mu2_1(:), Mu4_1(:), Mu2_2(:), Mu4_2(:)
  REAL(dp) :: buffer
  INTEGER  :: nThermoPoints
  CHARACTER(LEN=256) :: ThermoConditionsName='ThermoConditions.dat'

  !--- Variables for opacity read ---
  CHARACTER(LEN=256) :: OpTableName
  CHARACTER(LEN=256) :: Eos3DTableName
  CHARACTER(LEN=256) :: EosComposeTableName

  INTEGER                                 :: i, ii, jj

  ! Variables to calculate my general scattering       
  TYPE(HelmTableType) :: HelmTableMuons
  TYPE(HelmTableType) :: HelmTableElectrons
  TYPE(LeptonGasType) :: LeptonGasState
  CHARACTER(LEN=128)  :: process_string_1, process_string_2

  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: Phi0Out_General_1, &
                                             Phi0Out_General_2, &
                                             Phi1Out_General_1, &
                                             Phi1Out_General_2
  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: Phi0In_General_1, &
                                             Phi0In_General_2, &
                                             Phi1In_General_1, &
                                             Phi1In_General_2
  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: Rout_General_1, &
                                             Rout_General_2

  INTEGER , PARAMETER :: nL = 2
  INTEGER , PARAMETER :: iProcessMin = 1
  INTEGER , PARAMETER :: iProcessMax = 32
  REAL(dp) :: Rout_Int(iProcessMax - iProcessMin + 1)
  REAL(dp) :: Phout(iProcessMax - iProcessMin + 1, nL)
  REAL(dp) :: Phin(iProcessMax - iProcessMin + 1, nL)
  REAL(dp) :: E1, E3, Delta_mu, exponent
  INTEGER  :: iProcess, ierr, iProcess_1, iProcess_2
  INTEGER, PARAMETER :: nTheta = 24
  INTEGER(8) :: count_start, count_interm, count_end, count_rate
  REAL(dp)   :: t_General

  LOGICAL  :: Add_AntiLepton_Contribution = .TRUE.
  INTEGER  :: iProcess_nue_lep, iProcess_nueb_lep, iProcess_num_lep, iProcess_numb_lep
  INTEGER  :: iProcess_nue_alep, iProcess_nueb_alep, iProcess_num_alep, iProcess_numb_alep
  
  EosComposeTableName = "BaryonsPlusPhotonsPlusLeptonsEOS.h5"
  OpTableName         = "wl-Op-SFHo-15-25-50-E40-NES.h5"
  CALL InitializeHDF( )
  CALL ReadHelmholtzTableHDF(    &
            HelmTableElectrons,  &
            EosComposeTableName, &
            "HelmTableElectrons" )
  CALL ReadHelmholtzTableHDF(    &
            HelmTableMuons,      &
            EosComposeTableName, &
            "HelmTableMuons" )
  CALL FinalizeHDF( )

  !------------------------------------------------------
  !    Make energy grid
  !------------------------------------------------------
  CALL AllocateGrid( Inte_E, Inte_nPointE )

  Inte_E % Unit = 'MeV                  '
  Inte_E % Name = 'Intepolated Energy   '
  Inte_E % MinValue = Inte_Emin
  Inte_E % MaxValue = Inte_Emax
  Inte_E % LogInterp = 1
  Inte_E % nPoints = Inte_nPointE

  CALL MakeLogGrid &
          ( Inte_E % MinValue, Inte_E % MaxValue, &
            Inte_E % nPoints, Inte_E % Values )

  E1_arr = Inte_E % Values(:)
  E3_arr = Inte_E % Values(:nE3)

  E1_arr = [ (Inte_Emin + (Inte_Emax - Inte_Emin) * real(i - 1, dp) / real(Inte_nPointE - 1, dp), i = 1, nE1) ]
  E3_arr = E1_arr(:nE3)

  ! ------------------------------------------------------
  !    read in thermodynamic conditions
  ! ------------------------------------------------------
  OPEN(UNIT=123, FILE=trim(adjustl(ThermoConditionsName)), STATUS='OLD', ACTION='READ')
  READ(123,*) nThermoPoints
  READ(123,*)
  nThermoPoints = 1

  ALLOCATE(Inte_T  (nThermoPoints) , Inte_Rho(nThermoPoints), &
            Inte_Ye(nThermoPoints) , Inte_Mue(nThermoPoints), &
            Inte_Ym(nThermoPoints) , Inte_Mum(nThermoPoints), &
            Mu2_1  (nThermoPoints) , Mu4_1   (nThermoPoints), &
            Mu2_2  (nThermoPoints) , Mu4_2   (nThermoPoints) )

  DO i = 1, nThermoPoints
    
    READ(123,*) Inte_T(i), Inte_Rho(i), Inte_Ye(i), Inte_Ym(i), &
      buffer, buffer, buffer, buffer, buffer, buffer, buffer, buffer
    Inte_T(i) = Inte_T(i) * kmev_inv

    LeptonGasState % rho = Inte_Rho(i)
    LeptonGasState % T   = Inte_T(i)
    LeptonGasState % yL  = Inte_Ym(i)
    CALL LeptonGasEOS(HelmTableMuons, LeptonGasState)
    Inte_Mum(i) = LeptonGasState % mu

    LeptonGasState % rho = Inte_Rho(i)
    LeptonGasState % T   = Inte_T(i)
    LeptonGasState % yL  = Inte_Ye(i)
    CALL LeptonGasEOS(HelmTableElectrons, LeptonGasState)
    Inte_Mue(i) = LeptonGasState % mu

  END DO
  CLOSE(123)

  ! HARD CODE THERMODYNAMICAL STATE
  Inte_T(:)   = 15.0d0  * kmev_inv
  Inte_Mue(:) = 250.0d0
  Inte_Mum(:) = 0.0d0

  CALL InitGeneralScatteringKernels(  &
      E1_arr,                         &
      E3_arr,                         &
      nE1,                            &
      nE3,                            &
      nTheta,                         &
      iProcessMin_Option=iProcessMin, &
      iProcessMax_Option=iProcessMax  )
      
  ALLOCATE( Phi0Out_General_1 ( nE1, nE3, nThermoPoints ) )
  ALLOCATE( Phi0Out_General_2 ( nE1, nE3, nThermoPoints ) )
  ALLOCATE( Phi1Out_General_1 ( nE1, nE3, nThermoPoints ) )
  ALLOCATE( Phi1Out_General_2 ( nE1, nE3, nThermoPoints ) )

  ALLOCATE( Phi0In_General_1 ( nE1, nE3, nThermoPoints ) )
  ALLOCATE( Phi0In_General_2 ( nE1, nE3, nThermoPoints ) )
  ALLOCATE( Phi1In_General_1 ( nE1, nE3, nThermoPoints ) )
  ALLOCATE( Phi1In_General_2 ( nE1, nE3, nThermoPoints ) )

  ALLOCATE( Rout_General_1 ( nE1, nE3, nThermoPoints ) )
  ALLOCATE( Rout_General_2 ( nE1, nE3, nThermoPoints ) )

  ! ----------------------     LFC   -----------------------------
  process_string_1 = 'nu_mu      + mu+ -> nu_e       + e+  '
  Mu2_1(:) = -Inte_Mum(:); Mu4_1(:) = -Inte_Mue(:)
  CALL ProcessIndexFromReactionString( process_string_1, iProcess)
  iProcess_1 = iProcess - iProcessMin + 1
  
  process_string_2 = 'nu_e       + e+  -> nu_mu      + mu+ '
  Mu2_2(:) = -Inte_Mue(:); Mu4_2(:) = -Inte_Mum(:)
  CALL ProcessIndexFromReactionString( process_string_2, iProcess)
  iProcess_2 = iProcess - iProcessMin + 1

  ! ! ----------------------     LFE   -----------------------------
  ! process_string_1 = 'nu_bar_mu  + e+  -> nu_bar_e   + mu+ '
  ! Mu2_1(:) = -Inte_Mue(:); Mu4_1(:) = -Inte_Mum(:)
  ! CALL ProcessIndexFromReactionString( process_string_1, iProcess)
  ! iProcess_1 = iProcess - iProcessMin + 1
  
  ! process_string_2 = 'nu_bar_e   + mu+ -> nu_bar_mu  + e+  '
  ! Mu2_2(:) = -Inte_Mum(:); Mu4_2(:) = -Inte_Mue(:)
  ! CALL ProcessIndexFromReactionString( process_string_2, iProcess)
  ! iProcess_2 = iProcess - iProcessMin + 1

  ! ! ----------------------     LFE   -----------------------------
  ! process_string_1 = 'nu_mu      + e-  -> nu_e       + mu- '
  ! Mu2_1(:) = Inte_Mue(:); Mu4_1(:) = Inte_Mum(:)
  ! CALL ProcessIndexFromReactionString( process_string_1, iProcess)
  ! iProcess_1 = iProcess - iProcessMin + 1
  
  ! process_string_2 = 'nu_e       + mu- -> nu_mu      + e- '
  ! Mu2_2(:) = Inte_Mum(:); Mu4_2(:) = Inte_Mue(:)
  ! CALL ProcessIndexFromReactionString( process_string_2, iProcess)
  ! iProcess_2 = iProcess - iProcessMin + 1
  
  ! ! ----------------------     NMS   -----------------------------
  ! process_string_1 = 'nu_mu      + mu- -> nu_mu       + mu-  '
  ! Mu2_1(:) = Inte_Mum(:); Mu4_1(:) = Inte_Mum(:)
  ! CALL ProcessIndexFromReactionString( process_string_1, iProcess)
  ! iProcess_1 = iProcess - iProcessMin + 1
  
  ! process_string_2 = 'nu_mu      + mu- -> nu_mu       + mu-  '
  ! Mu2_2(:) = Inte_Mum(:); Mu4_2(:) = Inte_Mum(:)
  ! CALL ProcessIndexFromReactionString( process_string_2, iProcess)
  ! iProcess_2 = iProcess - iProcessMin + 1

  WRITE(*,*) 'Doing processes', iProcess_1, iProcess_2
  !============================================================================
  ! FIGURE 5: LOOP OVER ENERGY SCAN
  !============================================================================
  CALL SYSTEM_CLOCK(count_start, count_rate=count_rate)

#if defined(WEAKLIB_OMP)
  !$OMP PARALLEL DO PRIVATE(i, ii, jj, E1, E3, Delta_mu, exponent, Phin, Phout)
#endif
  DO i = 1, nThermoPoints
    ! Make sure you are using the exact same chemical potential
    ! WRITE(*,*) i, '/', nThermoPoints
    DO jj = 1, nE3
      DO ii = 1, nE1

        E1 = E1_arr(ii)
        E3 = E3_arr(jj)

        CALL CalculateAllPhout( ii, jj, E1, E3, Inte_T(i) * kmev, Inte_Mue(i), &
             Inte_Mum(i), Phout(:,:), nL )

        Delta_mu = Mu4_1(i) - Mu2_1(i)
        exponent = MIN( (Delta_mu + E3 - E1) / (Inte_T(i) * kmev), 500.0d0 )
        Phin = Phout * EXP(exponent)
             
        Phi0In_General_1(ii,jj,i) = Phin(iProcess_1,1)
        Phi0Out_General_1(ii,jj,i) = Phout(iProcess_1,1)

        Phi1In_General_1(ii,jj,i) = Phin(iProcess_1,2)
        Phi1Out_General_1(ii,jj,i) = Phout(iProcess_1,2)
        
        Delta_mu = Mu4_2(i) - Mu2_2(i)
        exponent = MIN( (Delta_mu + E3 - E1) / (Inte_T(i) * kmev), 500.0d0 )
        Phin = Phout * EXP(exponent)
             
        Phi0In_General_2(ii,jj,i) = Phin(iProcess_2,1)
        Phi0Out_General_2(ii,jj,i) = Phout(iProcess_2,1)

        Phi1In_General_2(ii,jj,i) = Phin(iProcess_2,2)
        Phi1Out_General_2(ii,jj,i) = Phout(iProcess_2,2)

        CALL CalculateAllRoutIntegrated( ii, jj, E1, E3, Inte_T(i) * kmev, Inte_Mue(i), &
             Inte_Mum(i), Rout_Int(:) )

        Rout_General_1(ii,jj,i) = Rout_Int(iProcess_1)
        Rout_General_2(ii,jj,i) = Rout_Int(iProcess_2)

      ENDDO
    ENDDO
  END DO
#if defined(WEAKLIB_OMP)
  !$OMP END PARALLEL DO
#endif

  CALL SYSTEM_CLOCK(count_end)
  t_General = REAL(count_end - count_start) / REAL(count_rate)

  WRITE(*,*) 't_General', t_General

  ! --- Output Routine ---
  OPEN(UNIT=66, FILE='CompareTwoReactions_PhiOut.dat', STATUS='REPLACE', ACTION='WRITE')
  WRITE(66,'(A)') '# ' // process_string_1
  WRITE(66,'(A)') '# ' // process_string_2
  WRITE(66,'(A)') '# T, Mue, Mum, E1, E3, ' // &
                   'Phi0_1, Phi1_1, Phi0_2, Phi1_2'

  ! --- Output Routine ---
  OPEN(UNIT=77, FILE='CompareTwoReactions_PhiIn.dat', STATUS='REPLACE', ACTION='WRITE')
  WRITE(77,'(A)') '# ' // process_string_1
  WRITE(77,'(A)') '# ' // process_string_2
  WRITE(77,'(A)') '# T, Mue, Mum, E1, E3, ' // &
                   'Phi0_1, Phi1_1, Phi0_2, Phi1_2'

  ! --- Output Routine ---
  OPEN(UNIT=88, FILE='CompareTwoReactions_Rout.dat', STATUS='REPLACE', ACTION='WRITE')
  WRITE(88,'(A)') '# ' // process_string_1
  WRITE(88,'(A)') '# ' // process_string_2
  WRITE(88,'(A)') '# T, Mue, Mum, E1, E3, ' // &
                   'Rout_1, Rout_2'

  ! THe order matters whatch out
  DO i = 1, nThermoPoints
    DO ii = 1, nE1
      DO jj = 1, nE3
        
        ! Retrieve the coordinate/thermodynamic variables for this row
        E1  = E1_arr(ii)
        E3  = E3_arr(jj)
        
        ! Write to CompareTwoReactions.dat (UNIT=66)
        WRITE(66, '(9(ES17.8E3, 1X))') &
          Inte_T(i) * kmev, Inte_Mue(i), Inte_Mum(i), E1, E3, &
          Phi0Out_General_1(ii,jj,i), Phi1Out_General_1(ii,jj,i), &
          Phi0Out_General_2(ii,jj,i), Phi1Out_General_2(ii,jj,i)

        WRITE(77, '(9(ES17.8E3, 1X))') &
          Inte_T(i) * kmev, Inte_Mue(i), Inte_Mum(i), E1, E3, &
          Phi0In_General_1(ii,jj,i), Phi1In_General_1(ii,jj,i), &
          Phi0In_General_2(ii,jj,i), Phi1In_General_2(ii,jj,i)

        WRITE(88, '(7(ES17.8E3, 1X))') &
          Inte_T(i) * kmev, Inte_Mue(i), Inte_Mum(i), E1, E3, &
          Rout_General_1(ii,jj,i), Rout_General_1(ii,jj,i)

      END DO
    END DO
  END DO

  ! Close the files when done
  CLOSE(66)
  CLOSE(77)
  CLOSE(88)

END PROGRAM wlTestGeneralKernelTwoReactions