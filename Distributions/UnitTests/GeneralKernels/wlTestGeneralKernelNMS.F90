PROGRAM wlTestGeneralKernelNMS

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
    InitGeneralScatteringKernels, &
    FinalizeGeneralScatteringKernels, &
    iProcessMin_Default, &
    iProcessMax_Default
  USE HDF5

  IMPLICIT NONE

  !--------- parameters for creating energy grid ----------------------------
  INTEGER, PARAMETER     :: Inte_nPointE = 40
  REAL(dp)               :: Inte_Emin = 1.0d-1
  REAL(dp)               :: Inte_Emax = 3.0d02
  TYPE(GridType)         :: Inte_E

  !--- Variables for profile read ---
  REAL(DP), ALLOCATABLE :: Inte_T(:) , Inte_Rho(:), Inte_Ym(:), &
                           Inte_Ye(:), Inte_Mue(:), Inte_Mum(:)
  REAL(DP) :: buffer
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
  CHARACTER(LEN=64)   :: process_string

  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: Phi0Out_General_nue_lep,    &
                                             Phi0Out_General_nueb_lep, &
                                             Phi0Out_General_num_lep,   &
                                             Phi0Out_General_numb_lep
  INTEGER , PARAMETER :: nL = 2
  INTEGER , PARAMETER :: iProcessMin = 13
  INTEGER , PARAMETER :: iProcessMax = 24
  REAL(DP) :: Phout(iProcessMax - iProcessMin + 1,nL)
  REAL(DP) :: Phin (iProcessMax - iProcessMin + 1,nL)
  REAL(DP) :: E1, E3, Delta_mu, exponent
  INTEGER  :: iProcess, ierr
  INTEGER, PARAMETER :: nTheta = 24
  INTEGER(8) :: count_start, count_interm, count_end, count_rate
  REAL(DP)   :: t_General

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

  ! ------------------------------------------------------
  !    read in thermodynamic conditions
  ! ------------------------------------------------------
  OPEN(UNIT=123, FILE=trim(adjustl(ThermoConditionsName)), STATUS='OLD', ACTION='READ')
  READ(123,*) nThermoPoints
  READ(123,*)
  nThermoPoints = 100

  ALLOCATE(Inte_T (nThermoPoints) , Inte_Rho(nThermoPoints) , &
            Inte_Ye(nThermoPoints) , Inte_Mue(nThermoPoints), &
            Inte_Ym(nThermoPoints) , Inte_Mum(nThermoPoints))

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

  CALL InitGeneralScatteringKernels(  &
      Inte_E % Values,                &
      Inte_E % Values,                &
      Inte_nPointE,                   &
      Inte_nPointE,                   &
      nTheta,                         &
      iProcessMin_Option=iProcessMin, &
      iProcessMax_Option=iProcessMax  )
      
  ALLOCATE( Phi0Out_General_nue_lep ( Inte_nPointE, Inte_nPointE, nThermoPoints ) )
  ALLOCATE( Phi0Out_General_nueb_lep( Inte_nPointE, Inte_nPointE, nThermoPoints ) )
  ALLOCATE( Phi0Out_General_num_lep ( Inte_nPointE, Inte_nPointE, nThermoPoints ) )
  ALLOCATE( Phi0Out_General_numb_lep( Inte_nPointE, Inte_nPointE, nThermoPoints ) )

  ! Map all indices correctly
  process_string = 'nu_e + mu- -> nu_e + mu-'
  CALL ProcessIndexFromReactionString( process_string, iProcess)
  iProcess_nue_lep = iProcess
  process_string = 'nu_bar_e + mu- -> nu_bar_e + mu-'
  CALL ProcessIndexFromReactionString( process_string, iProcess)
  iProcess_nueb_lep = iProcess
  process_string = 'nu_mu + mu- -> nu_mu + mu-'
  CALL ProcessIndexFromReactionString( process_string, iProcess)
  iProcess_num_lep = iProcess
  process_string = 'nu_bar_mu + mu- -> nu_bar_mu + mu-'
  CALL ProcessIndexFromReactionString( process_string, iProcess)
  iProcess_numb_lep = iProcess
  process_string = 'nu_e + mu+ -> nu_e + mu+'
  CALL ProcessIndexFromReactionString( process_string, iProcess)
  iProcess_nue_alep = iProcess
  process_string = 'nu_bar_e + mu+ -> nu_bar_e + mu+'
  CALL ProcessIndexFromReactionString( process_string, iProcess)
  iProcess_nueb_alep = iProcess
  process_string = 'nu_mu + mu+ -> nu_mu + mu+'
  CALL ProcessIndexFromReactionString( process_string, iProcess)
  iProcess_num_alep = iProcess
  process_string = 'nu_bar_mu + mu+ -> nu_bar_mu + mu+'
  CALL ProcessIndexFromReactionString( process_string, iProcess)
  iProcess_numb_alep = iProcess

  !============================================================================
  ! FIGURE 5: LOOP OVER ENERGY SCAN
  !============================================================================
  WRITE(*,*) 'Table done, now doing General Scattering Kernels'
  CALL SYSTEM_CLOCK(count_start, count_rate=count_rate)
  DO i = 1, nThermoPoints
    ! Make sure you are using the exact same chemical potential
    WRITE(*,*) i, '/', nThermoPoints
    DO ii = 1, Inte_nPointE
      DO jj = ii, Inte_nPointE

        E1 = Inte_E % Values(ii)
        E3 = Inte_E % Values(jj)

        CALL CalculateAllPhout( ii, jj, E1, E3, Inte_T(i) * kmev, Inte_Mue(i), &
             Inte_Mum(i), Phout(:,:), nL )

        Delta_mu = 0.0d0
        exponent = MIN( (E3 - E1 - Delta_mu) / (Inte_T(i) * kmev), 500.0d0 )
        Phin = Phout * EXP(exponent)
             
        Phi0Out_General_nue_lep(ii,jj,i) = Phin (iProcess_nue_lep,1)
        Phi0Out_General_nue_lep(jj,ii,i) = Phout(iProcess_nue_lep,1)

        Phi0Out_General_nueb_lep(ii,jj,i) = Phin (iProcess_nueb_lep,1)
        Phi0Out_General_nueb_lep(jj,ii,i) = Phout(iProcess_nueb_lep,1)

        Phi0Out_General_num_lep(ii,jj,i) = Phin (iProcess_num_lep,1)
        Phi0Out_General_num_lep(jj,ii,i) = Phout(iProcess_num_lep,1)

        Phi0Out_General_numb_lep(ii,jj,i) = Phin (iProcess_numb_lep,1)
        Phi0Out_General_numb_lep(jj,ii,i) = Phout(iProcess_numb_lep,1)

        WRITE(*,*) iProcess_nue_lep, Phout(iProcess_nue_lep,1)
        WRITE(*,*) iProcess_nueb_lep, Phout(iProcess_nueb_lep,1)
        WRITE(*,*) iProcess_num_lep, Phout(iProcess_num_lep,1)
        WRITE(*,*) iProcess_numb_lep, Phout(iProcess_numb_lep,1)
        STOP

      ENDDO
    ENDDO

  END DO
  CALL SYSTEM_CLOCK(count_end)
  t_General = REAL(count_end - count_start) / REAL(count_rate)

  WRITE(*,*) 't_General', t_General

END PROGRAM wlTestGeneralKernelNMS
