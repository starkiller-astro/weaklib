MODULE wlIOModuleCHIMERA

  USE wlKindModule, ONLY: dp
  USE HDF5
  USE wlThermoStateModule
  USE wlDependentVariablesModule
  USE wlEquationOfStateTableModule
  USE wlInterpolationModule
  USE wlExtNumericalModule, ONLY: epsilon, zero
  USE wlExtPhysicalConstantsModule
  USE e_p_eos_module
  USE EL_EOS_MODULE
  USE EOS_M4C_MODULE
  USE MAXWEL_MODULE
  USE EOS_BCK_MODULE
  USE wlExtEOSWrapperModule, ONLY: wlGetElectronEOS, wlGetFullEOS
  USE wlIOModuleHDF
  USE wlEOSIOModuleHDF
  !USE sfho_frdm_composition_module
  USE sfhx_frdm_composition_module
  !USE dd2_frdm_composition_module
  !USE iuf_roca_composition_module
  !USE fsg_roca_composition_module

  implicit none

  PUBLIC ReadChimeraProfile1D 
  PUBLIC ReadChimeraThermoProfile1D 
  PUBLIC ReadCHIMERAHDF
  PUBLIC ReadComposeTableHDF
  PUBLIC ReadSCTableHDF

CONTAINS

  SUBROUTINE ReadComposeTableHDF( EOSTable, FileName )
    
    IMPLICIT none

    TYPE(EquationOfStateTableType), INTENT(inout) :: EOSTable
    CHARACTER(len=*), INTENT(in)                  :: FileName

    INTEGER, DIMENSION(3)                         :: nPoints
    INTEGER, DIMENSION(3)                         :: nPointsBCK
    INTEGER, DIMENSION(1)                         :: nPointsTemp
    INTEGER, DIMENSION(1)                         :: TwoPoints
    INTEGER, DIMENSION(1)                         :: AllPoints
    INTEGER, DIMENSION(1)                         :: nThermo
    REAL(dp), DIMENSION(:), ALLOCATABLE           :: nb
    REAL(dp), DIMENSION(:), ALLOCATABLE           :: t
    REAL(dp), DIMENSION(:), ALLOCATABLE           :: yq
    REAL(dp), DIMENSION(:), ALLOCATABLE           :: thermo
    REAL(dp), DIMENSION(:), ALLOCATABLE           :: yi
    REAL(dp), DIMENSION(:), ALLOCATABLE           :: yav
    REAL(dp), DIMENSION(:), ALLOCATABLE           :: aav
    REAL(dp), DIMENSION(:), ALLOCATABLE           :: zav
    INTEGER                                       :: nVariables
    INTEGER                                       :: i, j, k, l
    INTEGER                                       :: iHe4
    INTEGER                                       :: BindingTableSwitch
    INTEGER                                       :: BCKTableSwitch
    INTEGER(HSIZE_T), DIMENSION(1)                :: datasize1d
    INTEGER(HID_T)                                :: file_id
    INTEGER(HID_T)                                :: group_id
    INTEGER :: hdferr
    INTEGER(HID_T)                                :: dataset_id
    REAL(dp)                                      :: chem_e      ! electron chemical potential [MeV]
    REAL(dp)                                      :: press_e     ! electron/positron/photon pressure
    REAL(dp)                                      :: entrop_e    ! electron/positron/photon entropy [kb/baryon]
    REAL(dp)                                      :: energ_e     ! electron/positron/photon energy
    REAL(dp)                                      :: press_buff     
    REAL(dp)                                      :: entrop_buff   
    REAL(dp)                                      :: energ_buff   
    REAL(dp)                                      :: pebuff   
    REAL(dp)                                      :: ppbuff   
    REAL(dp)                                      :: psbuff   
    REAL(dp)                                      :: minvar
    REAL(dp), DIMENSION(kmax)                     :: naz, xaz
    REAL(dp)                                      :: xn, xp, nn, np
    REAL(dp)                                      :: Ye_tmp
    INTEGER                                       :: sflag, count
    REAL(dp), DIMENSION(1) :: InterpolantU
    REAL(dp), DIMENSION(1) :: InterpolantP
    REAL(dp), DIMENSION(1,3) :: DerivativeU
    REAL(dp), DIMENSION(1,3) :: DerivativeP
    REAL(dp), DIMENSION(1)  :: rhobuff, tbuff, yebuff
    REAL(dp)  :: rhobuff2, tempbuff2, yebuff2
    REAL(dp), DIMENSION(1)  :: pbuff, ubuff, P1, P2, P1T, P2T, T1, T2, rho1, rho2, U1, U2
    REAL(dp)   :: press       ! pressure
    REAL(dp)   :: energ       ! internal energy
    REAL(dp)   :: entrop      ! entropy [kb/baryon]
    REAL(dp)   :: chem_n      ! free neutron chemical potential
    REAL(dp)   :: chem_p      ! free proton chemical potential
    REAL(dp)   :: xn_neut     ! free neutron fraction
    REAL(dp)   :: xn_prot     ! free proton fraction
    REAL(dp)   :: xn_alpha    ! alpha fraction
    REAL(dp)   :: xn_heavy    ! heavy fraction
    REAL(dp)   :: a_heavy     ! A for mean heavy nucleus
    REAL(dp)   :: z_heavy     ! Z for mean heavy nucleus
    REAL(dp)   :: be_heavy    ! Binding energy for mean heavy nucleus
    REAL(dp)   :: thermalenergy ! Internal energy without rest mass and constant offsets
    REAL(dp)   :: gamma1      ! Gamma from LS 
    REAL(dp)   :: pressbuff       ! pressure
    REAL(dp)   :: energbuff      ! internal energy
    REAL(dp)   :: entropbuff      ! entropy [kb/baryon]
    REAL(dp)   :: chem_nbuff      ! free neutron chemical potential
    REAL(dp)   :: chem_ebuff      ! electron chemical potential
    REAL(dp)   :: chem_pbuff      ! free proton chemical potential
    REAL(dp)   :: xn_neutbuff     ! free neutron fraction
    REAL(dp)   :: xn_protbuff     ! free proton fraction
    REAL(dp)   :: xn_alphabuff    ! alpha fraction
    REAL(dp)   :: xn_heavybuff    ! heavy fraction
    REAL(dp)   :: a_heavybuff     ! A for mean heavy nucleus
    REAL(dp)   :: z_heavybuff     ! Z for mean heavy nucleus
    REAL(dp)   :: be_heavybuff    ! Binding energy for mean heavy nucleus
    REAL(dp)   :: thermalenergybuff ! Internal energy without rest mass and constant offsets
    REAL(dp)   :: gamma1buff      ! Gamma from LS 
    REAL(dp)   :: xnbuff          ! xn buffer for light cluster read
    REAL(dp)   :: xpbuff          ! xn buffer for light cluster read
    REAL(dp)   :: Yd              ! Deuteron buffer
    REAL(dp)   :: Ytr             ! Tritium buffer
    REAL(dp)   :: Yhel            ! Helium 3 buffer
    REAL(dp)   :: xalphabuff      ! x alpha buffer for renormalization   
    REAL(dp)   :: xheavybuff      ! x heavy buffer for renormalization   
    REAL(dp)   :: totalmassfrac   ! total mass fraction for renormalization
    CHARACTER(len=1) :: flag     ! nuclear eos selection flag

    LOGICAL :: fail        ! did EoS fail to converge


    CALL OpenFileHDF( FileName, .false., file_id )

    BCKTableSwitch = 1

write(*,*) 'hdf5 table opened'

    datasize1d(1) = 1
    CALL ReadHDF( "pointsnb", nPointsTemp(:), file_id, datasize1d )
    nPoints(1) = nPointsTemp(1)

    datasize1d(1) = 1
    CALL ReadHDF( "pointst", nPointsTemp(:), file_id, datasize1d )
    nPoints(2) = nPointsTemp(1)


    datasize1d(1) = 1
    CALL ReadHDF( "pointsyq", nPointsTemp(:), file_id, datasize1d )
    nPoints(3) = nPointsTemp(1)

    CALL ReadHDF( "pointsqty", nThermo(:), file_id, datasize1d )
    nVariables = nThermo(1) + 9 ! ( (3 or 4)? mass frax, then heavy Z and heavy A)

write (*,*) 'nPoints', nPoints
write (*,*) 'nVariables', nVariables

    TwoPoints(1) = nPoints(1) * nPoints(2)
    AllPoints(1) = TwoPoints(1) * nPoints(3)

    nPointsBCK(1) = nPoints(1)
    nPointsBCK(2) = nPoints(2)
    IF ( BCKTableSwitch == 1 ) THEN
      nPointsBCK(3) = nPoints(3) + 40
      write (*,*) 'Expanded table range', nPointsBCK
    ELSE
      nPointsBCK(3) = nPoints(3)
      write (*,*) 'Standard table range', nPoints
    END IF

    IF ( BCKTableSwitch == 1 ) THEN
      CALL AllocateEquationOfStateTable( EOSTable, nPointsBCK , nVariables )
      write (*,*) 'Expanded EOSTable allocated'
    ELSE
      CALL AllocateEquationOfStateTable( EOSTable, nPoints , nVariables )
      write (*,*) 'Standard EOSTable allocated'
    END IF
 
    EOSTable % MD % IDTag = 'wl-EOS-SFHx+HiYeBCK-25-50-100, 5-11-17  '
    EOSTable % MD % TableResolution = '25 pts/dec rho, 50 pts/dec, delta ye = .01'
    EOSTable % MD % NucEOSLink = 'Nuc EOS Paper Link'
    EOSTable % MD % LeptonEOSLink = 'Lepton EOS Paper Link'
    EOSTable % MD % SourceLink = 'Table Source Link'
    EOSTable % MD % WLRevision = 'Rev'
    EOSTable % MD % TableLink = &
& 'http://eagle.phys.utk.edu/weaklib/trac/browser/External/Tables/EquationsOfState/wl-EOS-SFHo+HiYeBCK-25-50-100.h5'

    ALLOCATE( nb(0:(nPoints(1) - 1 ) ), t(0:(nPoints(2) - 1) ), yq( 0:(nPoints(3) - 1 ) ) )
    ALLOCATE( thermo(0:(nThermo(1)*AllPoints(1) - 1)) )
    ALLOCATE( yi(0:(3*AllPoints(1) - 1)), yav(0:(AllPoints(1) - 1)), &
      aav(0:(AllPoints(1) - 1)), zav( 0:(AllPoints(1) - 1) ) )
    
write (*,*) 'Buffers Allocated'
    datasize1d(1) = SIZE(nb)
    CALL ReadHDF( "nb", nb(:), file_id, datasize1d )

    datasize1d(1) = SIZE(t)
    CALL ReadHDF( "t", t(:), file_id, datasize1d )
write(*,*) 't', t

    datasize1d(1) = SIZE(yq)
    CALL ReadHDF( "yq", yq(:), file_id, datasize1d )

write (*,*) 'Buffers read'

    DO i = 0, nPoints(1) - 1
      EOSTable % TS % States(1) % Values(i+1) = nb(i) / kfm
    END DO
!write (*,*) 'rho', EOSTable % TS % States(1) % Values 

    DO i = 0, nPoints(2) - 1
      EOSTable % TS % States(2) % Values(i+1) = t(i) / kmev
    END DO
!write (*,*) 'T'

    DO i = 0, nPoints(3) - 1
      EOSTable % TS % States(3) % Values(i+1) = yq(i) 
    END DO
    
    IF ( BCKTableSwitch == 1 ) THEN
      DO i = nPoints(3) + 1, nPointsBCK(3)
        EOSTable % TS % States(3) % Values(i) = (1.0d-02)*i
write (*,*) 'ye', EOSTable % TS % States(3) % Values(i), i
      END DO
    END IF
!-----------------------------------------------------------
! Now that we've written the TS data, we need to fill in the
! additional TS metadata, like names, min/max, etc. 
!-----------------------------------------------------------

    EOSTable % TS % Names(1:3) = (/'Density                         ',&
                                   'Temperature                     ',&
                                   'Electron Fraction               '/)

    EOSTable % TS % Indices % iRho = 1
    EOSTable % TS % Indices % iT   = 2
    EOSTable % TS % Indices % iYe  = 3

write(*,*), "Allocate Independent Variable Units "

    EOSTable % TS % Units(1:3) = (/'Grams per cm^3                  ', &
                                   'K                               ', &
                                   '                                '/)

    EOSTable % TS % minValues(1) = EOSTable % TS % States(1) % Values(1)
    EOSTable % TS % minValues(2) = EOSTable % TS % States(2) % Values(1)
    EOSTable % TS % minValues(3) = EOSTable % TS % States(3) % Values(1)
    EOSTable % TS % maxValues(1) = EOSTable % TS % States(1) % Values(nPoints(1))
    EOSTable % TS % maxValues(2) = EOSTable % TS % States(2) % Values(nPoints(2))
    !EOSTable % TS % maxValues(3) = EOSTable % TS % States(3) % Values(nPoints(3))
    EOSTable % TS % maxValues(3) = EOSTable % TS % States(3) % Values(nPointsBCK(3))

    EOSTable % TS % LogInterp(1:3) =  (/1, 1, 0/)

write (*,*) 'TS filled'
    !CALL ReadDependentVariablesHDF( EOSTable % DV, file_id )

!-----------------------------------------------------------
! Now that we've written the TS, we need to read and write the DV
! In addition to the DV data, we need to fill in the
! additional DV metadata, like names, min/max, etc. 
! We also need to fill in the DVID.
!-----------------------------------------------------------

write(*,*), "Allocate Names "
    EOSTable % DV % Names(1:15) = (/'Pressure                        ', &
                                    'Entropy Per Baryon              ', &
                                    'Internal Energy Density         ', &
                                    'Electron Chemical Potential     ', &
                                    'Proton Chemical Potential       ', &
                                    'Neutron Chemical Potential      ', &
                                    'Proton Mass Fraction            ', &
                                    'Neutron Mass Fraction           ', &
                                    'Alpha Mass Fraction             ', &
                                    'Heavy Mass Fraction             ', &
                                    'Heavy Charge Number             ', &
                                    'Heavy Mass Number               ', &
                                    'Heavy Binding Energy            ', &
                                    'Thermal Energy                  ', &
                                    'Gamma1                          '/)

write(*,*), "Set Dependent Variable Identifier Indicies "

    EOSTable % DV % Indices % iPressure = 1
    EOSTable % DV % Indices % iEntropyPerBaryon = 2
    EOSTable % DV % Indices % iInternalEnergyDensity = 3
    EOSTable % DV % Indices % iElectronChemicalPotential = 4
    EOSTable % DV % Indices % iProtonChemicalPotential = 5
    EOSTable % DV % Indices % iNeutronChemicalPotential = 6
    EOSTable % DV % Indices % iProtonMassFraction = 7
    EOSTable % DV % Indices % iNeutronMassFraction = 8
    EOSTable % DV % Indices % iAlphaMassFraction = 9
    EOSTable % DV % Indices % iHeavyMassFraction = 10
    EOSTable % DV % Indices % iHeavyChargeNumber = 11
    EOSTable % DV % Indices % iHeavyMassNumber = 12
    EOSTable % DV % Indices % iHeavyBindingEnergy = 13
    EOSTable % DV % Indices % iThermalEnergy = 14
    EOSTable % DV % Indices % iGamma1 = 15

write(*,*), "Allocate Dependent Variable Units "
    EOSTable % DV % Units(1:15) = (/'Dynes per cm^2                  ', &
                                    'k_b per baryon                  ', &
                                    'erg per gram                    ', &
                                    'MeV                             ', &
                                    'MeV                             ', &
                                    'MeV                             ', &
                                    '                                ', &
                                    '                                ', &
                                    '                                ', &
                                    '                                ', &
                                    '                                ', &
                                    '                                ', &
                                    'MeV                             ', &
                                    'MeV                             ', &
                                    '                                '/)

write(*,*), "Allocate Dependent Variable Logical"
    EOSTable % DV % Repaired(:,:,:) = 0

    datasize1d(1) = SIZE(thermo)
    CALL ReadHDF( "thermo", thermo(:), file_id, datasize1d )
write(*,*), "thermo read"

    TwoPoints(1) = nPoints(1) * nPoints(2)
    AllPoints(1) = TwoPoints(1) * nPoints(3)

    DO k = 1, nPoints(3)
      DO j = 1, nPoints(2)
        DO i = 0, nPoints(1) - 1
          EOSTable % DV % Variables(1) % Values(i+1,j,k) &
            = thermo(i + nPoints(1)*(j-1) + TwoPoints(1)*(k-1)) * 1.6022e33 ! this is kp, or ergmev/cm3fm3 
        END DO
      END DO
    END DO
write (*,*) 'pressure(1,1,1)', EOSTable % DV % Variables(1) % Values(1,1,1)

    DO k = 1, nPoints(3)
      DO j = 1, nPoints(2)
        DO i = 0, nPoints(1) - 1
          EOSTable % DV % Variables(2) % Values(i+1,j,k) &
            = thermo((i + nPoints(1)*(j-1) + TwoPoints(1)*(k-1)) + AllPoints(1) )
        END DO
      END DO
    END DO
write (*,*) 'entropy(1,1,1)', EOSTable % DV % Variables(2) % Values(1,1,1)

    DO k = 1, nPoints(3)
      DO j = 1, nPoints(2)
        DO i = 0, nPoints(1) - 1
          EOSTable % DV % Variables(3) % Values(i+1,j,k) &
            = ergmev * mn * ( thermo((i + nPoints(1)*(j-1) + TwoPoints(1)*(k-1)) &
              + 4*AllPoints(1) ) + 8.9d0/mn ) * avn ! need to multiply by baryons per gram, avn 
        END DO
      END DO
    END DO
write (*,*) 'energy(1,1,1)', EOSTable % DV % Variables(3) % Values(1,1,1)

    DO k = 1, nPoints(3)
      DO j = 1, nPoints(2)
        DO i = 0, nPoints(1) - 1
          CALL wlGetElectronEOS( EOSTable % TS % States(1) % Values(i+1), & 
                                 EOSTable % TS % States(2) % Values(j),   &
                                 EOSTable % TS % States(3) % Values(k),   &
                                 press_e,  &
                                 entrop_e, &
                                 energ_e,  & 
                                 chem_e )

          
          pebuff = asig*( (kmev * EOSTable % TS % States(2) % Values(j) )**4 ) * ku & 
           / ( kfm * EOSTable % TS % States(1) % Values(i+1) )

          energ_buff = EOSTable % DV % Variables(3) % Values(i+1,j,k) 
          EOSTable % DV % Variables(3) % Values(i+1,j,k) &
            = energ_buff + energ_e + pebuff 

          ppbuff = EOSTable % TS % States(1) % Values(i+1) * pebuff / ( 3 ) 
          press_buff = EOSTable % DV % Variables(1) % Values(i+1,j,k)
          EOSTable % DV % Variables(1) % Values(i+1,j,k) &
            = press_buff + press_e + ppbuff

          psbuff = (4/3) * (pebuff/ku) / ( kmev * EOSTable % TS % States(2) % Values(j) )

          entrop_buff = EOSTable % DV % Variables(2) % Values(i+1,j,k)
          EOSTable % DV % Variables(2) % Values(i+1,j,k) &
            = entrop_buff + entrop_e + psbuff

          EOSTable % DV % Variables(4) % Values(i+1,j,k) &
            = chem_e 
        END DO
      END DO
    END DO
write (*,*) 'electron chem pot(1,1,1)', EOSTable % DV % Variables(4) % Values(1,1,1)

    DO k = 1, nPoints(3)
      DO j = 1, nPoints(2)
        DO i = 0, nPoints(1) - 1
          EOSTable % DV % Variables(6) % Values(i+1,j,k) &
            != ( thermo((i + nPoints(1)*(j-1) + TwoPoints(1)*(k-1)) + 2*AllPoints(1) ) + 1 ) * mn!- dmnp
            = ( thermo((i + nPoints(1)*(j-1) + TwoPoints(1)*(k-1)) + 2*AllPoints(1) ) ) - dmnp
        END DO
      END DO
    END DO
write (*,*) 'neutron chem pot(1,1,1)', EOSTable % DV % Variables(6) % Values(1,1,1)

    DO k = 1, nPoints(3)
      DO j = 1, nPoints(2)
        DO i = 0, nPoints(1) - 1
          EOSTable % DV % Variables(5) % Values(i+1,j,k) &
            = thermo((i + nPoints(1)*(j-1) + TwoPoints(1)*(k-1)) + 3*AllPoints(1) ) &
              + EOSTable % DV % Variables(6) % Values(i+1,j,k) + dmnp
        END DO
      END DO
    END DO
write (*,*) 'proton chem pot(1,1,1)', EOSTable % DV % Variables(5) % Values(1,1,1)
write(*,*), "thermo DV's filled" 

    datasize1d(1) = SIZE(yi)
    CALL ReadHDF( "yi", yi(:), file_id, datasize1d )
write(*,*), "yi read"
    DO k = 1, nPoints(3)
      DO j = 1, nPoints(2)
        DO i = 0, nPoints(1) - 1
          EOSTable % DV % Variables(8) % Values(i+1,j,k) &
            = MAX(yi(i + nPoints(1)*(j-1) + TwoPoints(1)*(k-1)),1.e-31)
        END DO
      END DO
    END DO
write (*,*) 'neutron mass frac(1,1,1)', EOSTable % DV % Variables(8) % Values(1,1,1)

    DO k = 1, nPoints(3)
      DO j = 1, nPoints(2)
        DO i = 0, nPoints(1) - 1
          EOSTable % DV % Variables(7) % Values(i+1,j,k) &
            = MAX(yi((i + nPoints(1)*(j-1) + TwoPoints(1)*(k-1)) + AllPoints(1) ), 1.e-31)
        END DO
      END DO
    END DO
write (*,*) 'proton mass frac(1,1,1)', EOSTable % DV % Variables(7) % Values(1,1,1)

    DO k = 1, nPoints(3)
      DO j = 1, nPoints(2)
        DO i = 0, nPoints(1) - 1
          EOSTable % DV % Variables(9) % Values(i+1,j,k) &
            = MAX(4*yi((i + nPoints(1)*(j-1) + TwoPoints(1)*(k-1)) + 2*AllPoints(1) ), 1.e-31 )
        END DO
      END DO
    END DO
write (*,*) 'alpha mass frac(1,1,1)', EOSTable % DV % Variables(9) % Values(1,1,1)


    DO k = 1, nPoints(3)
      DO j = 1, nPoints(2)
        DO i = 0, nPoints(1) - 1
            xpbuff = EOSTable % DV % Variables(7) % Values(i+1,j,k) ! xp
            xnbuff = EOSTable % DV % Variables(8) % Values(i+1,j,k) ! xn
            Yd = MAX(yi((i + nPoints(1)*(j-1) + TwoPoints(1)*(k-1)) + 3*AllPoints(1) ), 1.e-31 ) ! deuterons
            Ytr= MAX(yi((i + nPoints(1)*(j-1) + TwoPoints(1)*(k-1)) + 4*AllPoints(1) ), 1.e-31 ) ! tritons
            Yhel= MAX(yi((i + nPoints(1)*(j-1) + TwoPoints(1)*(k-1)) + 5*AllPoints(1) ), 1.e-31 ) ! helium3
            EOSTable % DV % Variables(7) % Values(i+1,j,k) = xpbuff + Yd + Ytr + 2*Yhel 
            EOSTable % DV % Variables(8) % Values(i+1,j,k) = xnbuff + Yd + 2*Ytr + Yhel 
        END DO
      END DO
    END DO
write (*,*) 'light cluster mass fractions added'


write(*,*), "yi DV's filled"

    datasize1d(1) = SIZE(zav)
    CALL ReadHDF( "zav", zav(:), file_id, datasize1d )
write(*,*), "zav read"
    DO k = 1, nPoints(3)
      DO j = 1, nPoints(2)
        DO i = 0, nPoints(1) - 1
          EOSTable % DV % Variables(11) % Values(i+1,j,k) = MAX(zav(i + nPoints(1)*(j-1) + TwoPoints(1)*(k-1)), 1.e-31)
        END DO
      END DO
    END DO
write (*,*) 'heavy charge number(1,1,1)', EOSTable % DV % Variables(11) % Values(1,1,1)
write(*,*), "zav DV's filled"

    datasize1d(1) = SIZE(aav)
    CALL ReadHDF( "aav", aav(:), file_id, datasize1d )
write(*,*), "aav read"
    DO k = 1, nPoints(3)
      DO j = 1, nPoints(2)
        DO i = 0, nPoints(1) - 1
          EOSTable % DV % Variables(12) % Values(i+1,j,k) = MAX(aav(i + nPoints(1)*(j-1) + TwoPoints(1)*(k-1)) + epsilon, 1.e-31)
        END DO
      END DO
    END DO
write (*,*) 'heavy mass number(1,1,1)', EOSTable % DV % Variables(12) % Values(1,1,1)
write(*,*), "aav DV filled"

    datasize1d(1) = SIZE(yav)
    CALL ReadHDF( "yav", yav(:), file_id, datasize1d )
write(*,*), "yav read"
    DO k = 1, nPoints(3)
      DO j = 1, nPoints(2)
        DO i = 0, nPoints(1) - 1
          EOSTable % DV % Variables(10) % Values(i+1,j,k) &
            = MAX(yav(i + nPoints(1)*(j-1) + TwoPoints(1)*(k-1)) & 
            * EOSTable % DV % Variables(12) % Values(i+1,j,k) + epsilon, 1.e-31 )
        END DO
      END DO
    END DO
write (*,*) 'heavy mass frac(1,1,1)', EOSTable % DV % Variables(10) % Values(1,1,1)
write(*,*), "yav DV filled"

    DO k = 1, nPoints(3)
      DO j = 1, nPoints(2)
        DO i = 0, nPoints(1) - 1
            xpbuff = EOSTable % DV % Variables(7) % Values(i+1,j,k) ! xp
            xnbuff = EOSTable % DV % Variables(8) % Values(i+1,j,k) ! xn
            xalphabuff = EOSTable % DV % Variables(9) % Values(i+1,j,k) ! xalpha
            xheavybuff = EOSTable % DV % Variables(10) % Values(i+1,j,k) ! xheavy
             
            totalmassfrac = xpbuff + xnbuff + xalphabuff + xheavybuff

            EOSTable % DV % Variables(7) % Values(i+1,j,k) = xpbuff/totalmassfrac
            EOSTable % DV % Variables(8) % Values(i+1,j,k) = xnbuff/totalmassfrac
            EOSTable % DV % Variables(9) % Values(i+1,j,k) = xalphabuff/totalmassfrac
            EOSTable % DV % Variables(10) % Values(i+1,j,k) = xheavybuff/totalmassfrac
        END DO
      END DO
    END DO

  IF ( BCKTableSwitch == 1 ) THEN 
 write(*,*) "starting bck loop"

    DO k = nPointsBCK(3), nPoints(3) + 1, -1
      DO j = 1, nPoints(2)
        DO i = 1, nPoints(1)

          flag = 'B'
          rhobuff2 = EOSTable % TS % States(1) % Values(i)
          tempbuff2 = EOSTable % TS % States(2) % Values(j)
          yebuff2 = EOSTable % TS % States(3) % Values(k)

          IF ( rhobuff2 > 1.0e11 ) THEN
            EOSTable % DV % Variables(1) % Values(i,j,k) = 0.
            EOSTable % DV % Variables(2) % Values(i,j,k) = 0.
            EOSTable % DV % Variables(3) % Values(i,j,k) = 0.
            EOSTable % DV % Variables(4) % Values(i,j,k) = 0.
            EOSTable % DV % Variables(5) % Values(i,j,k) = 0.
            EOSTable % DV % Variables(6) % Values(i,j,k) = 0.
            EOSTable % DV % Variables(7) % Values(i,j,k) = 0.
            EOSTable % DV % Variables(8) % Values(i,j,k) = 0.
            EOSTable % DV % Variables(9) % Values(i,j,k) = 0.
            EOSTable % DV % Variables(10) % Values(i,j,k) = 0.
            EOSTable % DV % Variables(11) % Values(i,j,k) = 0.
            EOSTable % DV % Variables(12) % Values(i,j,k) = 0.
            EOSTable % DV % Variables(13) % Values(i,j,k) = 0.
            EOSTable % DV % Variables(14) % Values(i,j,k) = 0.
            EOSTable % DV % Variables(15) % Values(i,j,k) = 0.
            EOSTable % DV % Repaired(i,j,k) = 2
          ELSE

            CALL wlGetFullEOS( rhobuff2, tempbuff2, yebuff2, flag, fail, press, energ,       &
              & entrop, chem_n, chem_p, chem_e, xn_neut, xn_prot, xn_alpha, xn_heavy, &
              & a_heavy, z_heavy, be_heavy, thermalenergy, gamma1, i, j, k )

              EOSTable % DV % Variables(1) % Values(i,j,k) = press
              EOSTable % DV % Variables(2) % Values(i,j,k) = entrop
              EOSTable % DV % Variables(3) % Values(i,j,k) = energ
              EOSTable % DV % Variables(4) % Values(i,j,k) = chem_e
              EOSTable % DV % Variables(5) % Values(i,j,k) = chem_p
              EOSTable % DV % Variables(6) % Values(i,j,k) = chem_n
              EOSTable % DV % Variables(7) % Values(i,j,k) = xn_prot
              EOSTable % DV % Variables(8) % Values(i,j,k) = xn_neut
              EOSTable % DV % Variables(9) % Values(i,j,k) = xn_alpha
              EOSTable % DV % Variables(10) % Values(i,j,k) = xn_heavy
              EOSTable % DV % Variables(11) % Values(i,j,k) = z_heavy
              EOSTable % DV % Variables(12) % Values(i,j,k) = a_heavy
              EOSTable % DV % Variables(13) % Values(i,j,k) = be_heavy
              EOSTable % DV % Variables(14) % Values(i,j,k) = thermalenergy
          !    EOSTable % DV % Variables(15) % Values(i,j,k) = gamma1
          END IF
        END DO
      END DO
    END DO

    DO k = nPointsBCK(3), nPoints(3) + 1, -1
      DO j = 1, nPoints(2)
        DO i = 1, nPoints(1)
          IF ( rhobuff2 > 1.0e11 ) THEN
            EOSTable % DV % Variables(15) % Values(i,j,k) = 0. 
          ELSE
              rhobuff(1) = EOSTable % TS % States(1) % Values(i)
              tbuff(1)   = EOSTable % TS % States(2) % Values(j)
              pbuff(1)   = EOSTable % DV % Variables(1) % Values(i,j,k)
              ubuff(1)   = EOSTable % DV % Variables(3) % Values(i,j,k)
              P1(1) = LOG10( EOSTable % DV % Variables(1) % Values(i-1,j,k) )
              P2(1)= LOG10( EOSTable % DV % Variables(1) % Values(i+1,j,k) )
              P1T(1) = LOG10( EOSTable % DV % Variables(1) % Values(i,j-1,k) )
              P2T(1)= LOG10( EOSTable % DV % Variables(1) % Values(i,j+1,k) )
              U1(1) = LOG10( EOSTable % DV % Variables(3) % Values(i,j-1,k) )
              U2(1) = LOG10( EOSTable % DV % Variables(3) % Values(i,j+1,k) )
              rho1(1) = LOG10( EOSTable % TS % States(1) % Values(i-1) )
              rho2(1) = LOG10( EOSTable % TS % States(1) % Values(i+1) )
              T1(1) = LOG10( EOSTable % TS % States(2) % Values(j-1) )
              T2(1) = LOG10( EOSTable % TS % States(2) % Values(j+1) )
              DerivativeP(1,3) = 0.d0
              DerivativeU(1,3) = 0.d0
              DerivativeU(1,1) = 0.d0

              IF ( i == 1 ) THEN
              P1(1) = LOG10( EOSTable % DV % Variables(1) % Values(i,j,k) )
              rho1(1) = LOG10( EOSTable % TS % States(1) % Values(i) )
              END IF

              !IF ( i == nPoints(1) ) THEN
              IF ( EOSTable % TS % States(1) % Values(i+1) > 1.0e11 ) THEN
              P2(1) = LOG10( EOSTable % DV % Variables(1) % Values(i,j,k) )
              rho2(1) = LOG10( EOSTable % TS % States(1) % Values(i) )
              END IF
    
              IF( j == 1 ) THEN
              P1T(1) = LOG10( EOSTable % DV % Variables(1) % Values(i,j,k) )
              U1(1) = LOG10( EOSTable % DV % Variables(3) % Values(i,j,k) )
              T1(1) = LOG10( EOSTable % TS % States(2) % Values(j) )
              END IF

              IF( j == nPoints(2) ) THEN
              P2T(1)= LOG10( EOSTable % DV % Variables(1) % Values(i,j,k) )
              U2(1) = LOG10( EOSTable % DV % Variables(3) % Values(i,j,k) )
              T2(1) = LOG10( EOSTable % TS % States(2) % Values(j) )
              END IF

              DerivativeP(1,1) = (P2(1) - P1(1))/(rho2(1) - rho1(1)) * (pbuff(1)/rhobuff(1))
              DerivativeP(1,2) = (P2T(1) - P1T(1))/(T2(1) - T1(1)) * (pbuff(1)/tbuff(1))
              DerivativeU(1,2) = (U2(1) - U1(1))/(T2(1) - T1(1)) * (ubuff(1)/tbuff(1))

              EOSTable % DV % Variables(15) % Values(i,j,k) &
              = MAX( ( rhobuff(1) * DerivativeP(1,1) + ( tbuff(1)/(rhobuff(1) + epsilon) ) &
                    * ((DerivativeP(1,2))**2)/(DerivativeU(1,2)) )/( pbuff(1) + epsilon ), 1.e-31 )

          END IF
      END DO
    END DO
  END DO
write(*,*) 'BCK loop finished'
END IF

write (*,*) 'nPoints', nPoints
    DO k = 1, nPoints(3) 
      DO j = 1, nPoints(2)
        DO i = 1, nPoints(1)

          rhobuff(1) = EOSTable % TS % States(1) % Values(i)
          tbuff(1)   = EOSTable % TS % States(2) % Values(j)
          pbuff(1)   = EOSTable % DV % Variables(1) % Values(i,j,k)
          ubuff(1)   = EOSTable % DV % Variables(3) % Values(i,j,k)
          P1(1) = LOG10( EOSTable % DV % Variables(1) % Values(i-1,j,k) )
          P2(1)= LOG10( EOSTable % DV % Variables(1) % Values(i+1,j,k) )
          P1T(1) = LOG10( EOSTable % DV % Variables(1) % Values(i,j-1,k) )
          P2T(1)= LOG10( EOSTable % DV % Variables(1) % Values(i,j+1,k) )
          U1(1) = LOG10( EOSTable % DV % Variables(3) % Values(i,j-1,k) )
          U2(1) = LOG10( EOSTable % DV % Variables(3) % Values(i,j+1,k) )
          rho1(1) = LOG10( EOSTable % TS % States(1) % Values(i-1) )
          rho2(1) = LOG10( EOSTable % TS % States(1) % Values(i+1) )
          T1(1) = LOG10( EOSTable % TS % States(2) % Values(j-1) )
          T2(1) = LOG10( EOSTable % TS % States(2) % Values(j+1) )
          DerivativeP(1,3) = 0.d0
          DerivativeU(1,3) = 0.d0
          DerivativeU(1,1) = 0.d0
          
          IF ( i == 1 ) THEN
          P1(1) = LOG10( EOSTable % DV % Variables(1) % Values(i,j,k) )
          rho1(1) = LOG10( EOSTable % TS % States(1) % Values(i) )
          END IF

          IF ( i == nPoints(1) ) THEN
          P2(1) = LOG10( EOSTable % DV % Variables(1) % Values(i,j,k) )
          rho2(1) = LOG10( EOSTable % TS % States(1) % Values(i) )
          END IF

          IF( j == 1 ) THEN
          P1T(1) = LOG10( EOSTable % DV % Variables(1) % Values(i,j,k) )
          U1(1) = LOG10( EOSTable % DV % Variables(3) % Values(i,j,k) )
          T1(1) = LOG10( EOSTable % TS % States(2) % Values(j) )
          END IF

          IF( j == nPoints(2) ) THEN
          P2T(1)= LOG10( EOSTable % DV % Variables(1) % Values(i,j,k) )
          U2(1) = LOG10( EOSTable % DV % Variables(3) % Values(i,j,k) )
          T2(1) = LOG10( EOSTable % TS % States(2) % Values(j) )
          END IF

          DerivativeP(1,1) = (P2(1) - P1(1))/(rho2(1) - rho1(1)) * (pbuff(1)/rhobuff(1))
          DerivativeP(1,2) = (P2T(1) - P1T(1))/(T2(1) - T1(1)) * (pbuff(1)/tbuff(1))
          DerivativeU(1,2) = (U2(1) - U1(1))/(T2(1) - T1(1)) * (ubuff(1)/tbuff(1))

          EOSTable % DV % Variables(15) % Values(i,j,k) &
          = MAX( ( rhobuff(1) * DerivativeP(1,1) + ( tbuff(1)/(rhobuff(1) + epsilon) ) &
                * ((DerivativeP(1,2))**2)/(DerivativeU(1,2)) )/( pbuff(1) + epsilon ), 1.e-31 )

        END DO
      END DO
    END DO

    CALL compdata_readin
!do  i = 1, 3162
!write(*,*) i, az(i, 1), az(i,2), bind(i) 
!end do
    BindingTableSwitch = 1 
write(*,*) 'starting binding table loops'
    IF ( BindingTableSwitch == 0 ) THEN 

    DO k = 1, nPoints(3)
      DO j = 1, nPoints(2)
        DO i = 1, nPoints(1)
          EOSTable % DV % Variables(13) % Values(i,j,k) = 0.0d0 
        END DO
      END DO
    END DO

    ELSE IF ( BindingTableSwitch == 2 ) THEN 

write(*,*) 'starting roca-maza loops'

!------------------------------------
! Start Roca-maza (IUFSU, FSUGold) Binding Energy Read Block
!------------------------------------
    DO k = 1, nPoints(3)
      DO j = 1, nPoints(2)
        DO i = 1, nPoints(1)
          !CALL sub_dist_grid(j,k,i,xaz,xn,xp,naz,nn,np)

          CALL sub_dist_interpol(t(j-1),yq(k-1),nb(i-1),xaz,xn,xp,naz,nn,np,sflag)
             IF ( ( SUM( az(:,1)*naz(:)/nb(i-1) ) - az(1515,1)*naz(1515)/nb(i-1) &
                 - az(1514,1)*naz(1514)/nb(i-1) - az(1513,1)*naz(1513)/nb(i-1) &
                 - az(1512,1)*naz(1512)/nb(i-1) ) .le. 0.0d0 ) THEN
              !EOSTable % DV % Variables(10) % Values(i,j,k) = 0.0d0
              EOSTable % DV % Variables(11) % Values(i,j,k) = 1.8d0
              EOSTable % DV % Variables(12) % Values(i,j,k) = 4.0d0
              EOSTable % DV % Variables(13) % Values(i,j,k) = 0.0d0
            ELSE

            EOSTable % DV % Variables(13) % Values(i,j,k) &
            = ( SUM( -bind(:) * xaz(:) / ( az(:,1) ) )       & !+ bind(iHe4) * xaz(iHe4) / 4.d0   &     
               + bind(1515) * xaz(1515) / az(1515,1)    &
               + bind(1514) * xaz(1514) / az(1514,1)    &
               + bind(1513) * xaz(1513) / az(1513,1)    &
               + bind(1512) * xaz(1512) / az(1512,1) )  &
               / ( SUM( az(:,1)*naz(:)/nb(i-1) ) - az(1515,1)*naz(1515)/nb(i-1) &
               - az(1514,1)*naz(1514)/nb(i-1) - az(1513,1)*naz(1513)/nb(i-1) &
               - az(1512,1)*naz(1512)/nb(i-1) )
          END IF
        END DO
      END DO
    END DO

!------------------------------------
! END Roca-maza (IUFSU, FSU-Gold) Binding Energy Read Block
!------------------------------------

    ELSE IF ( BindingTableSwitch == 1 ) THEN

write(*,*) "starting sfho frdm loop"
!write(*,*) "starting sfhx frdm loop"

!------------------------------------
! Start FDRM (SFHo/x, DD2) Binding Energy Read Block
!------------------------------------
    DO k = 1, nPoints(3)
      DO j = 1, nPoints(2) - 1
        DO i = 1, nPoints(1)
      !  DO i = 1, nPoints(1)
      !DO j = 1, nPoints(2)
          CALL sub_dist_interpol(t(j-1),yq(k-1),nb(i-1),xaz,xn,xp,naz,nn,np,sflag)
      !    CALL sub_dist_grid(j,k,i,xaz,xn,xp,naz,nn,np)
             IF ( ( SUM( az(:,1)*naz(:)/nb(i-1) ) - az(8077,1)*naz(8077)/nb(i-1) &
                 - az(8076,1)*naz(8076)/nb(i-1) - az(8075,1)*naz(8075)/nb(i-1) &
                 - az(8074,1)*naz(8074)/nb(i-1) ) .le. 0.0d0 ) THEN
              !EOSTable % DV % Variables(10) % Values(i,j,k) = 0.0d0
              EOSTable % DV % Variables(11) % Values(i,j,k) = 1.8d0
              EOSTable % DV % Variables(12) % Values(i,j,k) = 4.0d0
              EOSTable % DV % Variables(13) % Values(i,j,k) = 0.0d0
            ELSE

            EOSTable % DV % Variables(13) % Values(i,j,k) &
            = ( SUM( -bind(:) * xaz(:) / ( az(:,1) ) )       & !+ bind(iHe4) * xaz(iHe4) / 4.d0   &     
               + bind(8077) * xaz(8077) / az(8077,1)    &
               + bind(8076) * xaz(8076) / az(8076,1)    &
               + bind(8075) * xaz(8075) / az(8075,1)    &
               + bind(8074) * xaz(8074) / az(8074,1) )  &
               / ( SUM( az(:,1)*naz(:)/nb(i-1) ) - az(8077,1)*naz(8077)/nb(i-1) &
               - az(8076,1)*naz(8076)/nb(i-1) - az(8075,1)*naz(8075)/nb(i-1) &
               - az(8074,1)*naz(8074)/nb(i-1) )
          END IF
        END DO
write(*,*) i-1,j,k,t(j-1),"Binding Energy=", EOSTable % DV % Variables(13) % Values(i-1,j,k) 
      END DO
    END DO

    DO k = 1, nPoints(3)
        DO i = 1, nPoints(1) 
      !    CALL sub_dist_interpol(t(j-1),yq(k-1),nb(i-1),xaz,xn,xp,naz,nn,np,sflag)
           j = nPoints(2) !81
          CALL sub_dist_grid(j,k,i,xaz,xn,xp,naz,nn,np)
             IF ( ( SUM( az(:,1)*naz(:)/nb(i-1) ) - az(8077,1)*naz(8077)/nb(i-1) &
                 - az(8076,1)*naz(8076)/nb(i-1) - az(8075,1)*naz(8075)/nb(i-1) &
                 - az(8074,1)*naz(8074)/nb(i-1) ) .le. 0.0d0 ) THEN
              !EOSTable % DV % Variables(10) % Values(i,j,k) = 0.0d0
              EOSTable % DV % Variables(11) % Values(i,j,k) = 1.8d0
              EOSTable % DV % Variables(12) % Values(i,j,k) = 4.0d0
              EOSTable % DV % Variables(13) % Values(i,j,k) = 0.0d0
            ELSE

            EOSTable % DV % Variables(13) % Values(i,j,k) &
            = ( SUM( -bind(:) * xaz(:) / ( az(:,1) ) )       & !+ bind(iHe4) * xaz(iHe4) / 4.d0   &     
               + bind(8077) * xaz(8077) / az(8077,1)    &
               + bind(8076) * xaz(8076) / az(8076,1)    &
               + bind(8075) * xaz(8075) / az(8075,1)    &
               + bind(8074) * xaz(8074) / az(8074,1) )  &
               / ( SUM( az(:,1)*naz(:)/nb(i-1) ) - az(8077,1)*naz(8077)/nb(i-1) &
               - az(8076,1)*naz(8076)/nb(i-1) - az(8075,1)*naz(8075)/nb(i-1) &
               - az(8074,1)*naz(8074)/nb(i-1) )
          END IF
      END DO
write(*,*) i-1,nPoints(2),k,t(nPoints(2)-1),"Binding Energy max=", EOSTable % DV % Variables(13) % Values(i-1,npoints(2),k)
    END DO


write(*,*) 'binding energy block complete' 
!------------------------------------
! End FDRM (SFHo/x, DD2) Binding Energy Read Block
!------------------------------------

    END IF

    DO k = 1, nPoints(3)
      DO j = 1, nPoints(2)
        DO i = 1, nPoints(1)

          EOSTable % DV % Variables(14) % Values(i,j,k) =                    &
            & EOSTable % DV % Variables(3) % Values(i,j,k)                   &
            & - EOSTable % DV % Variables(13) % Values(i,j,k)                & 
            & + ku * ( dmnp * EOSTable % DV % Variables(7) % Values(i,j,k)   &
            & + 7.075 * EOSTable % DV % Variables(9) % Values(i,j,k)         & 
            & + 1.5d0 * kmev * EOSTable % TS % States(2) % Values(j)         &
            & * ( EOSTable % DV % Variables(10) % Values(i,j,k)              &
            & / EOSTable % DV % Variables(12) % Values(i,j,k) )              & 
            & - ( EOSTable % TS % States(3) % Values(k) * me ) )
        END DO
      END DO
    END DO

IF ( BCKTableSwitch == 1) THEN
write(*,*) "starting sfho-bck merge loop"
!
    k = 60!nPoints(3) 
    DO j = 1, nPoints(2)
      DO i = 1, nPoints(1)

      flag = 'B'
       rhobuff2 = EOSTable % TS % States(1) % Values(i) 
       tempbuff2 = EOSTable % TS % States(2) % Values(j) 
       yebuff2 = EOSTable % TS % States(3) % Values(k) 


          IF ( rhobuff2 > 1.0e11 ) THEN
            CYCLE
          ELSE
              pressbuff = EOSTable % DV % Variables(1) % Values(i,j,k) 
              entropbuff = EOSTable % DV % Variables(2) % Values(i,j,k) 
              energbuff = EOSTable % DV % Variables(3) % Values(i,j,k) 
              chem_ebuff = EOSTable % DV % Variables(4) % Values(i,j,k) 
              chem_pbuff = EOSTable % DV % Variables(5) % Values(i,j,k) 
              chem_nbuff = EOSTable % DV % Variables(6) % Values(i,j,k) 
              xn_protbuff = EOSTable % DV % Variables(7) % Values(i,j,k) 
              xn_neutbuff = EOSTable % DV % Variables(8) % Values(i,j,k) 
              xn_alphabuff = EOSTable % DV % Variables(9) % Values(i,j,k) 
              xn_heavybuff = EOSTable % DV % Variables(10) % Values(i,j,k) 
              z_heavybuff = EOSTable % DV % Variables(11) % Values(i,j,k) 
              a_heavybuff = EOSTable % DV % Variables(12) % Values(i,j,k) 
              be_heavybuff = EOSTable % DV % Variables(13) % Values(i,j,k) 
              thermalenergybuff = EOSTable % DV % Variables(14) % Values(i,j,k) 

            CALL wlGetFullEOS( rhobuff2, tempbuff2, yebuff2, flag, fail, press, energ,       &
              & entrop, chem_n, chem_p, chem_e, xn_neut, xn_prot, xn_alpha, xn_heavy, &
              & a_heavy, z_heavy, be_heavy, thermalenergy, gamma1, i, j, k )

              EOSTable % DV % Variables(1) % Values(i,j,k) &
              = (pressbuff + press)/2
              EOSTable % DV % Variables(2) % Values(i,j,k) &
              = (entropbuff + entrop)/2
              EOSTable % DV % Variables(3) % Values(i,j,k) &               
              = (energbuff + energ)/2
              EOSTable % DV % Variables(4) % Values(i,j,k) &               
              = (chem_ebuff + chem_e)/2
              EOSTable % DV % Variables(5) % Values(i,j,k) &               
              = (chem_pbuff + chem_p)/2
              EOSTable % DV % Variables(6) % Values(i,j,k) &               
              = (chem_nbuff + chem_n)/2
              EOSTable % DV % Variables(7) % Values(i,j,k) &               
              =(xn_protbuff + xn_prot)/2
              EOSTable % DV % Variables(8) % Values(i,j,k) &               
              =(xn_neutbuff + xn_neut)/2
              EOSTable % DV % Variables(9) % Values(i,j,k) &               
              =(xn_alphabuff + xn_alpha)/2
              EOSTable % DV % Variables(10) % Values(i,j,k) &               
              =(xn_heavybuff + xn_heavy)/2
              EOSTable % DV % Variables(11) % Values(i,j,k) &               
              =(z_heavybuff + z_heavy)/2
              EOSTable % DV % Variables(12) % Values(i,j,k) &               
              =(a_heavybuff + a_heavy)/2
              EOSTable % DV % Variables(13) % Values(i,j,k) &               
              =(be_heavybuff + be_heavy)/2
              EOSTable % DV % Variables(14) % Values(i,j,k) &               
              =(thermalenergybuff + thermalenergy)/2
write(*,*) 'P-SFHo, P-BCK', pressbuff, press

            END IF
          END DO
        END DO
  

    k = 60!nPoints(3) 
    DO j = 1, nPoints(2)
      DO i = 1, nPoints(1)

       rhobuff2 = EOSTable % TS % States(1) % Values(i) 


          IF ( rhobuff2 > 1.0e11 ) THEN
            CYCLE
          ELSE
              rhobuff(1) = EOSTable % TS % States(1) % Values(i)
              tbuff(1)   = EOSTable % TS % States(2) % Values(j)
              pbuff(1)   = EOSTable % DV % Variables(1) % Values(i,j,k)
              ubuff(1)   = EOSTable % DV % Variables(3) % Values(i,j,k)
              P1(1) = LOG10( EOSTable % DV % Variables(1) % Values(i-1,j,k) )
              P2(1)= LOG10( EOSTable % DV % Variables(1) % Values(i+1,j,k) )
              P1T(1) = LOG10( EOSTable % DV % Variables(1) % Values(i,j-1,k) )
              P2T(1)= LOG10( EOSTable % DV % Variables(1) % Values(i,j+1,k) )
              U1(1) = LOG10( EOSTable % DV % Variables(3) % Values(i,j-1,k) )
              U2(1) = LOG10( EOSTable % DV % Variables(3) % Values(i,j+1,k) )
              rho1(1) = LOG10( EOSTable % TS % States(1) % Values(i-1) )
              rho2(1) = LOG10( EOSTable % TS % States(1) % Values(i+1) )
              T1(1) = LOG10( EOSTable % TS % States(2) % Values(j-1) )
              T2(1) = LOG10( EOSTable % TS % States(2) % Values(j+1) )
              DerivativeP(1,3) = 0.d0
              DerivativeU(1,3) = 0.d0
              DerivativeU(1,1) = 0.d0

              IF ( i == 1 ) THEN
              P1(1) = LOG10( EOSTable % DV % Variables(1) % Values(i,j,k) )
              rho1(1) = LOG10( EOSTable % TS % States(1) % Values(i) )
              END IF

              !IF ( i == nPoints(1) ) THEN
              IF ( EOSTable % TS % States(1) % Values(i+1) > 1.0e11 ) THEN
              P2(1) = LOG10( EOSTable % DV % Variables(1) % Values(i,j,k) )
              rho2(1) = LOG10( EOSTable % TS % States(1) % Values(i) )
              END IF

              IF( j == 1 ) THEN
              P1T(1) = LOG10( EOSTable % DV % Variables(1) % Values(i,j,k) )
              U1(1) = LOG10( EOSTable % DV % Variables(3) % Values(i,j,k) )
              T1(1) = LOG10( EOSTable % TS % States(2) % Values(j) )
              END IF

              IF( j == nPoints(2) ) THEN
              P2T(1)= LOG10( EOSTable % DV % Variables(1) % Values(i,j,k) )
              U2(1) = LOG10( EOSTable % DV % Variables(3) % Values(i,j,k) )
              T2(1) = LOG10( EOSTable % TS % States(2) % Values(j) )
              END IF

              DerivativeP(1,1) = (P2(1) - P1(1))/(rho2(1) - rho1(1)) * (pbuff(1)/rhobuff(1))
              DerivativeP(1,2) = (P2T(1) - P1T(1))/(T2(1) - T1(1)) * (pbuff(1)/tbuff(1))
              DerivativeU(1,2) = (U2(1) - U1(1))/(T2(1) - T1(1)) * (ubuff(1)/tbuff(1))

              EOSTable % DV % Variables(15) % Values(i,j,k) &
              = MAX( ( rhobuff(1) * DerivativeP(1,1) + ( tbuff(1)/(rhobuff(1) + epsilon) ) &
                    * ((DerivativeP(1,2))**2)/(DerivativeU(1,2)) )/( pbuff(1) + epsilon ), 1.e-31 )
          END IF
      END DO
    END DO

!write(*,*) i,j,k,"Energy=", EOSTable % DV % Variables(3) % Values(i,j,k) 


write(*,*) 'SFHo-BCK merge loop finished'
END IF

    DO l = 1, EOSTable % nVariables
      EOSTable % DV % minValues(l) = MINVAL( EOSTable % DV % Variables(l) % Values )
      EOSTable % DV % maxValues(l) = MAXVAL( EOSTable % DV % Variables(l) % Values )
    END DO

  DO l = 1, 15 !EOSTable % nVariables
    WRITE (*,*) EOSTable % DV % Names(l)
    minvar = MINVAL( EOSTable % DV % Variables(l) % Values )
    WRITE (*,*) "minvar=", minvar
    EOSTable % DV % Offsets(l) = -2.d0 * MIN( 0.d0, minvar )
    WRITE (*,*) "Offset=", EOSTable % DV % Offsets(l)
    EOSTable % DV % Variables(l) % Values &
      = LOG10( EOSTable % DV % Variables(l) % Values &
               + EOSTable % DV % Offsets(l) + epsilon )

  END DO

    CALL WriteEquationOfStateTableHDF( EOSTable )

    CALL DescribeEquationOfStateTable( EOSTable )

    !CALL CloseFileHDF( file_id )

  END SUBROUTINE ReadComposeTableHDF

  SUBROUTINE ReadSCTableHDF( EOSTable, FileName )
    
    IMPLICIT none

    TYPE(EquationOfStateTableType), INTENT(inout) :: EOSTable
    CHARACTER(len=*), INTENT(in)                  :: FileName

    INTEGER, DIMENSION(3)                         :: nPoints
    INTEGER, DIMENSION(1)                         :: nPointsTemp
    INTEGER, DIMENSION(1)                         :: TwoPoints
    INTEGER, DIMENSION(1)                         :: AllPoints
    INTEGER, DIMENSION(1)                         :: nThermo
    REAL(dp), DIMENSION(1)                        :: energyshift
    REAL(dp), DIMENSION(:), ALLOCATABLE           :: logrho
    REAL(dp), DIMENSION(:), ALLOCATABLE           :: logtemp
    REAL(dp), DIMENSION(:), ALLOCATABLE           :: ye
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE       :: logpress
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE       :: logenergy
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE       :: entropy
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE       :: mu_e
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE       :: mu_n
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE       :: mu_p
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE       :: Xa
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE       :: Xpp
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE       :: Xnn
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE       :: Xh
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE       :: Zbar
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE       :: Abar
    REAL(dp), DIMENSION(:), ALLOCATABLE           :: thermo
    REAL(dp), DIMENSION(:), ALLOCATABLE           :: yi
    REAL(dp), DIMENSION(:), ALLOCATABLE           :: yav
    REAL(dp), DIMENSION(:), ALLOCATABLE           :: aav
    REAL(dp), DIMENSION(:), ALLOCATABLE           :: zav
    INTEGER                                       :: nVariables
    INTEGER                                       :: i, j, k, l, TestUnit1, TestUnit2
    INTEGER                                       :: iHe4
    INTEGER(HSIZE_T), DIMENSION(1)                :: datasize1d
    INTEGER(HSIZE_T), DIMENSION(3)                :: datasize3d
    INTEGER(HID_T)                                :: file_id
    INTEGER(HID_T)                                :: group_id
    INTEGER :: hdferr
    INTEGER(HID_T)                                :: dataset_id
    REAL(dp)                                      :: chem_e      ! electron chemical potential [MeV]
    REAL(dp)                                      :: press_e     ! electron/positron/photon pressure
    REAL(dp)                                      :: entrop_e    ! electron/positron/photon entropy [kb/baryon]
    REAL(dp)                                      :: energ_e     ! electron/positron/photon energy
    REAL(dp)                                      :: press_buff     
    REAL(dp)                                      :: entrop_buff   
    REAL(dp)                                      :: energ_buff   
    REAL(dp)                                      :: pebuff   
    REAL(dp)                                      :: ppbuff   
    REAL(dp)                                      :: psbuff   
    REAL(dp)                                      :: minvar
    REAL(dp), DIMENSION(kmax)                     :: naz, xaz
    REAL(dp)                                      :: xn, xp, nn, np
    INTEGER :: sflag
    REAL(dp), DIMENSION(1) :: InterpolantU
    REAL(dp), DIMENSION(1) :: InterpolantP
    REAL(dp), DIMENSION(1,3) :: DerivativeU
    REAL(dp), DIMENSION(1,3) :: DerivativeP
    REAL(dp), DIMENSION(1)  :: rhobuff, tbuff, yebuff
    REAL(dp), DIMENSION(1)  :: pbuff, ubuff, P1, P2, P1T, P2T, T1, T2, rho1, rho2, U1, U2


    OPEN( newunit = TestUnit1, FILE="BindingEnergyTable.d")
    OPEN( newunit = TestUnit2, FILE="HeavyFracCheck.d")

    CALL OpenFileHDF( FileName, .false., file_id )

write(*,*) 'hdf5 table opened'

    datasize1d(1) = 1
    CALL ReadHDF( "pointsrho", nPointsTemp(:), file_id, datasize1d )
    nPoints(1) = nPointsTemp(1)

    datasize1d(1) = 1
    CALL ReadHDF( "pointstemp", nPointsTemp(:), file_id, datasize1d )
    nPoints(2) = nPointsTemp(1)

    datasize1d(1) = 1
    CALL ReadHDF( "pointsye", nPointsTemp(:), file_id, datasize1d )
    nPoints(3) = nPointsTemp(1)

    nVariables = 15!nThermo(1) + 9 ! ( (3 or 4)? mass frax, then heavy Z and heavy A)

write (*,*) 'nPoints', nPoints
write (*,*) 'nVariables', nVariables

    CALL AllocateEquationOfStateTable( EOSTable, nPoints , nVariables )
write (*,*) 'EOSTable allocated'

    ALLOCATE( logrho(0:(nPoints(1) - 1 ) ), logtemp(0:(nPoints(2) - 1) ), ye( 0:(nPoints(3) - 1 ) ) )
    
write (*,*) 'Buffers Allocated'
    datasize1d(1) = SIZE(logrho)
    CALL ReadHDF( "logrho", logrho(:), file_id, datasize1d )

    datasize1d(1) = SIZE(logtemp)
    CALL ReadHDF( "logtemp", logtemp(:), file_id, datasize1d )

    datasize1d(1) = SIZE(ye)
    CALL ReadHDF( "ye", ye(:), file_id, datasize1d )

write (*,*) 'Buffers read'

    DO i = 0, nPoints(1) - 1
      EOSTable % TS % States(1) % Values(i+1) = 10**logrho(i)
    END DO
write (*,*) 'rho', EOSTable % TS % States(1) % Values 

    DO i = 0, nPoints(2) - 1
      EOSTable % TS % States(2) % Values(i+1) = 10**logtemp(i) / kmev
    END DO
write (*,*) 'T', EOSTable % TS % States(2) % Values 

    DO i = 0, nPoints(3) - 1
      EOSTable % TS % States(3) % Values(i+1) = ye(i) 
    END DO
!write (*,*) 'ye'
write (*,*) 'TS Loaded'
!-----------------------------------------------------------
! Now that we've written the TS data, we need to fill in the
! additional TS metadata, like names, min/max, etc. 
!-----------------------------------------------------------

    EOSTable % TS % Names(1:3) = (/'Density                         ',&
                                   'Temperature                     ',&
                                   'Electron Fraction               '/)

    EOSTable % TS % Indices % iRho = 1
    EOSTable % TS % Indices % iT   = 2
    EOSTable % TS % Indices % iYe  = 3

write(*,*), "Allocate Independent Variable Units "

    EOSTable % TS % Units(1:3) = (/'Grams per cm^3                  ', &
                                   'K                               ', &
                                   '                                '/)

    EOSTable % TS % minValues(1) = EOSTable % TS % States(1) % Values(1)
    EOSTable % TS % minValues(2) = EOSTable % TS % States(2) % Values(1)
    EOSTable % TS % minValues(3) = EOSTable % TS % States(3) % Values(1)
    EOSTable % TS % maxValues(1) = EOSTable % TS % States(1) % Values(nPoints(1))
    EOSTable % TS % maxValues(2) = EOSTable % TS % States(2) % Values(nPoints(2))
    EOSTable % TS % maxValues(3) = EOSTable % TS % States(3) % Values(nPoints(3))

    EOSTable % TS % LogInterp(1:3) =  (/1, 1, 0/)

write (*,*) 'TS filled'

!-----------------------------------------------------------
! Now that we've written the TS, we need to read and write the DV
! In addition to the DV data, we need to fill in the
! additional DV metadata, like names, min/max, etc. 
! We also need to fill in the DVID.
!-----------------------------------------------------------

write(*,*), "Allocate Names "
    EOSTable % DV % Names(1:15) = (/'Pressure                        ', &
                                    'Entropy Per Baryon              ', &
                                    'Internal Energy Density         ', &
                                    'Electron Chemical Potential     ', &
                                    'Proton Chemical Potential       ', &
                                    'Neutron Chemical Potential      ', &
                                    'Proton Mass Fraction            ', &
                                    'Neutron Mass Fraction           ', &
                                    'Alpha Mass Fraction             ', &
                                    'Heavy Mass Fraction             ', &
                                    'Heavy Charge Number             ', &
                                    'Heavy Mass Number               ', &
                                    'Heavy Binding Energy            ', &
                                    'Thermal Energy                  ', &
                                    'Gamma1                          '/)

write(*,*), "Set Dependent Variable Identifier Indicies "

    EOSTable % DV % Indices % iPressure = 1
    EOSTable % DV % Indices % iEntropyPerBaryon = 2
    EOSTable % DV % Indices % iInternalEnergyDensity = 3
    EOSTable % DV % Indices % iElectronChemicalPotential = 4
    EOSTable % DV % Indices % iProtonChemicalPotential = 5
    EOSTable % DV % Indices % iNeutronChemicalPotential = 6
    EOSTable % DV % Indices % iProtonMassFraction = 7
    EOSTable % DV % Indices % iNeutronMassFraction = 8
    EOSTable % DV % Indices % iAlphaMassFraction = 9
    EOSTable % DV % Indices % iHeavyMassFraction = 10
    EOSTable % DV % Indices % iHeavyChargeNumber = 11
    EOSTable % DV % Indices % iHeavyMassNumber = 12
    EOSTable % DV % Indices % iHeavyBindingEnergy = 13
    EOSTable % DV % Indices % iThermalEnergy = 14
    EOSTable % DV % Indices % iGamma1 = 15

write(*,*), "Allocate Dependent Variable Units "
    EOSTable % DV % Units(1:15) = (/'Dynes per cm^2                  ', &
                                    'k_b per baryon                  ', &
                                    'erg per gram                    ', &
                                    'MeV                             ', &
                                    'MeV                             ', &
                                    'MeV                             ', &
                                    '                                ', &
                                    '                                ', &
                                    '                                ', &
                                    '                                ', &
                                    '                                ', &
                                    '                                ', &
                                    'MeV                             ', &
                                    'MeV                             ', &
                                    '                                '/)

write(*,*), "Allocate Dependent Variable Logical"
    EOSTable % DV % Repaired(:,:,:) = 0
write(*,*), "repaired step"

ALLOCATE( logpress(0:(nPoints(1)-1), 0:(nPoints(2)-1), 0:(nPoints(3)-1) ) )
write(*,*), "logpress allocated"

    datasize3d = SHAPE( logpress )
write(*,*), "datasize shaped"
    CALL ReadHDF( "logpress", logpress(:,:,:), &
                              file_id, datasize3d )
write(*,*), "logpress read from table"

    DO k = 1, nPoints(3)
      DO j = 1, nPoints(2)
        DO i = 1, nPoints(1)
          EOSTable % DV % Variables(1) % Values(i,j,k) = 10**logpress(i-1,j-1,k-1)
        END DO
      END DO
    END DO
write (*,*) 'pressure(1,1,1)', EOSTable % DV % Variables(1) % Values(1,1,1)

ALLOCATE( entropy(0:(nPoints(1)-1), 0:(nPoints(2)-1), 0:(nPoints(3)-1) ) )
write(*,*), "entropy allocated"

    datasize3d = SHAPE( entropy )
write(*,*), "datasize shaped"
    CALL ReadHDF( "entropy", entropy(:,:,:), &
                              file_id, datasize3d )
    DO k = 1, nPoints(3)
      DO j = 1, nPoints(2)
        DO i = 1, nPoints(1)
          EOSTable % DV % Variables(2) % Values(i,j,k) &
            = entropy(i-1,j-1,k-1)
        END DO
      END DO
    END DO
write (*,*) 'entropy(1,1,1)', EOSTable % DV % Variables(2) % Values(1,1,1)

ALLOCATE( logenergy(0:(nPoints(1)-1), 0:(nPoints(2)-1), 0:(nPoints(3)-1) ) )
write(*,*), "logenergy allocated"

    datasize1d(1) = 1
    CALL ReadHDF( "energy_shift", energyshift(:), file_id, datasize1d )

    datasize3d = SHAPE( logenergy )
write(*,*), "datasize shaped"
    CALL ReadHDF( "logenergy", logenergy(:,:,:), &
                              file_id, datasize3d )
write(*,*), "logenergy read from table"

    DO k = 1, nPoints(3)
      DO j = 1, nPoints(2)
        DO i = 1, nPoints(1)
          EOSTable % DV % Variables(3) % Values(i,j,k) = 10**logenergy(i-1,j-1,k-1) - energyshift(1)
        END DO
      END DO
    END DO
write (*,*) 'energy(1,1,1)', EOSTable % DV % Variables(3) % Values(1,1,1)

ALLOCATE( mu_e(0:(nPoints(1)-1), 0:(nPoints(2)-1), 0:(nPoints(3)-1) ) )
write(*,*), "mu_e allocated"

    datasize3d = SHAPE( mu_e )
write(*,*), "datasize shaped"
    CALL ReadHDF( "mu_e", mu_e(:,:,:), &
                              file_id, datasize3d )
write(*,*), "mu_e read from table"


    DO k = 1, nPoints(3)
      DO j = 1, nPoints(2)
        DO i = 1, nPoints(1)
          EOSTable % DV % Variables(4) % Values(i,j,k) &
            = mu_e(i-1,j-1,k-1)
        END DO
      END DO
    END DO
write (*,*) 'electron chem pot(1,1,1)', EOSTable % DV % Variables(4) % Values(1,1,1)

ALLOCATE( mu_n(0:(nPoints(1)-1), 0:(nPoints(2)-1), 0:(nPoints(3)-1) ) )
write(*,*), "mu_n allocated"

    datasize3d = SHAPE( mu_n )
write(*,*), "datasize shaped"
    CALL ReadHDF( "mu_n", mu_n(:,:,:), &
                              file_id, datasize3d )
write(*,*), "mu_n read from table"

    DO k = 1, nPoints(3)
      DO j = 1, nPoints(2)
        DO i = 1, nPoints(1)
          EOSTable % DV % Variables(6) % Values(i,j,k) = mu_n(i-1,j-1,k-1)
        END DO
      END DO
    END DO
write (*,*) 'neutron chem pot(1,1,1)', EOSTable % DV % Variables(6) % Values(1,1,1)

ALLOCATE( mu_p(0:(nPoints(1)-1), 0:(nPoints(2)-1), 0:(nPoints(3)-1) ) )
write(*,*), "mu_p allocated"

    datasize3d = SHAPE( mu_p )
write(*,*), "datasize shaped"
    CALL ReadHDF( "mu_p", mu_p(:,:,:), &
                              file_id, datasize3d )
write(*,*), "mu_p read from table"

    DO k = 1, nPoints(3)
      DO j = 1, nPoints(2)
        DO i = 1, nPoints(1) 
          EOSTable % DV % Variables(6) % Values(i,j,k) = mu_p(i-1,j-1,k-1)    
        END DO
      END DO
    END DO
    
write (*,*) 'proton chem pot(1,1,1)', EOSTable % DV % Variables(5) % Values(1,1,1)

ALLOCATE( Xnn(0:(nPoints(1)-1), 0:(nPoints(2)-1), 0:(nPoints(3)-1) ) )
write(*,*), "Xnn allocated"

    datasize3d = SHAPE( Xnn )
write(*,*), "datasize shaped"
    CALL ReadHDF( "Xn", Xnn(:,:,:), &
                              file_id, datasize3d )
write(*,*), "Xn read from table"
    DO k = 1, nPoints(3)
      DO j = 1, nPoints(2)
        DO i = 1, nPoints(1)
          EOSTable % DV % Variables(8) % Values(i,j,k) &
            = Xnn(i-1,j-1,k-1)
        END DO
      END DO
    END DO
write (*,*) 'neutron mass frac(1,1,1)', EOSTable % DV % Variables(8) % Values(1,1,1)

ALLOCATE( Xpp(0:(nPoints(1)-1), 0:(nPoints(2)-1), 0:(nPoints(3)-1) ) )
write(*,*), "Xpp allocated"

    datasize3d = SHAPE( Xpp ) 
write(*,*), "datasize shaped"
    CALL ReadHDF( "Xp", Xpp(:,:,:), & 
                              file_id, datasize3d )
write(*,*), "Xpp read from table"
    DO k = 1, nPoints(3)
      DO j = 1, nPoints(2)
        DO i = 1, nPoints(1) 
          EOSTable % DV % Variables(7) % Values(i,j,k) &
            = Xpp(i-1,j-1,k-1)
        END DO
      END DO
    END DO
write (*,*) 'proton mass frac(1,1,1)', EOSTable % DV % Variables(7) % Values(1,1,1)

ALLOCATE( Xa(0:(nPoints(1)-1), 0:(nPoints(2)-1), 0:(nPoints(3)-1) ) )
write(*,*), "Xa allocated"

    datasize3d = SHAPE( Xa ) 
write(*,*), "datasize shaped"
    CALL ReadHDF( "Xa", Xa(:,:,:), & 
                              file_id, datasize3d )
write(*,*), "Xa read from table"
    DO k = 1, nPoints(3)
      DO j = 1, nPoints(2)
        DO i = 1, nPoints(1) 
          EOSTable % DV % Variables(9) % Values(i,j,k) &
            = Xa(i-1,j-1,k-1)
        END DO
      END DO
    END DO
write (*,*) 'alpha mass frac(1,1,1)', EOSTable % DV % Variables(9) % Values(1,1,1)
write(*,*), "yi DV's filled"

ALLOCATE( Zbar(0:(nPoints(1)-1), 0:(nPoints(2)-1), 0:(nPoints(3)-1) ) )
write(*,*), "Zbar allocated"

    datasize3d = SHAPE( Zbar ) 
write(*,*), "datasize shaped"
    CALL ReadHDF( "Zbar", Zbar(:,:,:), & 
                              file_id, datasize3d )
write(*,*), "Zbar read from table"
    DO k = 1, nPoints(3)
      DO j = 1, nPoints(2)
        DO i = 1, nPoints(1) 
          EOSTable % DV % Variables(11) % Values(i,j,k) &
            = Zbar(i-1,j-1,k-1)
        END DO
      END DO
    END DO
write (*,*) 'heavy charge number(1,1,1)', EOSTable % DV % Variables(11) % Values(1,1,1)
write(*,*), "zav DV's filled"

ALLOCATE( Abar(0:(nPoints(1)-1), 0:(nPoints(2)-1), 0:(nPoints(3)-1) ) )
write(*,*), "Abar allocated"

    datasize3d = SHAPE( Abar ) 
write(*,*), "datasize shaped"
    CALL ReadHDF( "Abar", Abar(:,:,:), & 
                              file_id, datasize3d )
write(*,*), "Abar read from table"
    DO k = 1, nPoints(3)
      DO j = 1, nPoints(2)
        DO i = 1, nPoints(1) 
          EOSTable % DV % Variables(12) % Values(i,j,k) &
            = Abar(i-1,j-1,k-1)
        END DO
      END DO
    END DO
write (*,*) 'heavy mass number(1,1,1)', EOSTable % DV % Variables(12) % Values(1,1,1)

ALLOCATE( Xh(0:(nPoints(1)-1), 0:(nPoints(2)-1), 0:(nPoints(3)-1) ) )
write(*,*), "Xh allocated"

    datasize3d = SHAPE( Xh ) 
write(*,*), "datasize shaped"
    CALL ReadHDF( "Xh", Xh(:,:,:), & 
                              file_id, datasize3d )
write(*,*), "Xh read from table"
    DO k = 1, nPoints(3)
      DO j = 1, nPoints(2)
        DO i = 1, nPoints(1) 
          EOSTable % DV % Variables(10) % Values(i,j,k) &
            = Xh(i-1,j-1,k-1)
        END DO
      END DO
    END DO
write (*,*) 'heavy mass frac(1,1,1)', EOSTable % DV % Variables(10) % Values(1,1,1)
write(*,*), "yav DV filled"

write (*,*) 'nPoints', nPoints
    DO k = 1, nPoints(3) 
      DO j = 1, nPoints(2)
        DO i = 1, nPoints(1)

          rhobuff(1) = EOSTable % TS % States(1) % Values(i)
          tbuff(1)   = EOSTable % TS % States(2) % Values(j)
          pbuff(1)   = EOSTable % DV % Variables(1) % Values(i,j,k)
          ubuff(1)   = EOSTable % DV % Variables(3) % Values(i,j,k)
          P1(1) = LOG10( EOSTable % DV % Variables(1) % Values(i-1,j,k) )
          P2(1)= LOG10( EOSTable % DV % Variables(1) % Values(i+1,j,k) )
          P1T(1) = LOG10( EOSTable % DV % Variables(1) % Values(i,j-1,k) )
          P2T(1)= LOG10( EOSTable % DV % Variables(1) % Values(i,j+1,k) )
          U1(1) = LOG10( EOSTable % DV % Variables(3) % Values(i,j-1,k) )
          U2(1) = LOG10( EOSTable % DV % Variables(3) % Values(i,j+1,k) )
          rho1(1) = LOG10( EOSTable % TS % States(1) % Values(i-1) )
          rho2(1) = LOG10( EOSTable % TS % States(1) % Values(i+1) )
          T1(1) = LOG10( EOSTable % TS % States(2) % Values(j-1) )
          T2(1) = LOG10( EOSTable % TS % States(2) % Values(j+1) )
          DerivativeP(1,3) = 0.d0
          DerivativeU(1,3) = 0.d0
          DerivativeU(1,1) = 0.d0
          
          IF ( i == 1 ) THEN
          P1(1) = LOG10( EOSTable % DV % Variables(1) % Values(i,j,k) )
          rho1(1) = LOG10( EOSTable % TS % States(1) % Values(i) )
          END IF

          IF ( i == nPoints(1) ) THEN
          P2(1) = LOG10( EOSTable % DV % Variables(1) % Values(i,j,k) )
          rho2(1) = LOG10( EOSTable % TS % States(1) % Values(i) )
          END IF

          IF( j == 1 ) THEN
          P1T(1) = LOG10( EOSTable % DV % Variables(1) % Values(i,j,k) )
          U1(1) = LOG10( EOSTable % DV % Variables(3) % Values(i,j,k) )
          T1(1) = LOG10( EOSTable % TS % States(2) % Values(j) )
          END IF

          IF( j == nPoints(2) ) THEN
          P2T(1)= LOG10( EOSTable % DV % Variables(1) % Values(i,j,k) )
          U2(1) = LOG10( EOSTable % DV % Variables(3) % Values(i,j,k) )
          T2(1) = LOG10( EOSTable % TS % States(2) % Values(j) )
          END IF

          DerivativeP(1,1) = (P2(1) - P1(1))/(rho2(1) - rho1(1)) * (pbuff(1)/rhobuff(1))
          DerivativeP(1,2) = (P2T(1) - P1T(1))/(T2(1) - T1(1)) * (pbuff(1)/tbuff(1))
          DerivativeU(1,2) = (U2(1) - U1(1))/(T2(1) - T1(1)) * (ubuff(1)/tbuff(1))

          EOSTable % DV % Variables(15) % Values(i,j,k) &
          = MAX( ( rhobuff(1) * DerivativeP(1,1) + ( tbuff(1)/(rhobuff(1) + epsilon) ) &
                * ((DerivativeP(1,2))**2)/(DerivativeU(1,2)) )/( pbuff(1) + epsilon ), 1.e-31 )

        END DO
      END DO
    END DO

    CALL compdata_readin

    DO k = 1, nPoints(3)
write(*,*) "binding energy: Ye=", EOSTable % TS % States(3) % Values(k)
      DO j = 1, nPoints(2)
        DO i = 1, nPoints(1)
           CALL sub_dist_interpol((10**logtemp(j-1)),ye(k-1),(10**logrho(i-1) * kfm),xaz,xn,xp,naz,nn,np,sflag)
             IF ( ( SUM( az(:,1)*naz(:)/(10**logrho(i-1) * kfm ) ) - az(8077,1)*naz(8077)/( 10**logrho(i-1) * kfm ) &
                 - az(8076,1)*naz(8076)/(10**logrho(i-1) * kfm ) - az(8075,1)*naz(8075)/(10**logrho(i-1) * kfm ) &
                 - az(8074,1)*naz(8074)/(10**logrho(i-1) * kfm ) ) .le. 0.0d0 ) THEN
              EOSTable % DV % Variables(11) % Values(i,j,k) = 1.8d0
              EOSTable % DV % Variables(12) % Values(i,j,k) = 4.0d0
              EOSTable % DV % Variables(13) % Values(i,j,k) = 0.0d0
            ELSE

            EOSTable % DV % Variables(13) % Values(i,j,k) &
            = ( SUM( -bind(:) * xaz(:) / ( az(:,1) ) )       & !+ bind(iHe4) * xaz(iHe4) / 4.d0   &     
               + bind(8077) * xaz(8077) / az(8077,1)    &
               + bind(8076) * xaz(8076) / az(8076,1)    &
               + bind(8075) * xaz(8075) / az(8075,1)    &
               + bind(8074) * xaz(8074) / az(8074,1) )  &
               / ( SUM( az(:,1)*naz(:)/(10**logrho(i-1) * kfm ) ) - az(8077,1)*naz(8077)/(10**logrho(i-1) * kfm ) &
               - az(8076,1)*naz(8076)/(10**logrho(i-1) * kfm ) - az(8075,1)*naz(8075)/(10**logrho(i-1) * kfm ) &
               - az(8074,1)*naz(8074)/(10**logrho(i-1) * kfm ) )
          END IF
        END DO
      END DO
    END DO

    DO k = 1, nPoints(3)
      DO j = 1, nPoints(2)
        DO i = 1, nPoints(1)

          EOSTable % DV % Variables(14) % Values(i,j,k) =                    &
            & EOSTable % DV % Variables(3) % Values(i,j,k)                   &
            & - EOSTable % DV % Variables(13) % Values(i,j,k)                & ! BINDING ENERGY SHOULD BE OUTSIDE PARENTHESES 
            & + ku * ( dmnp * EOSTable % DV % Variables(7) % Values(i,j,k)   &
            & + 7.075 * EOSTable % DV % Variables(9) % Values(i,j,k)         & 
            !& - EOSTable % DV % Variables(13) % Values(i,j,k)                & ! BINDING ENERGY SHOULD BE OUTSIDE PARENTHESES 
            & + 1.5d0 * kmev * EOSTable % TS % States(2) % Values(j)         &
            & * ( EOSTable % DV % Variables(10) % Values(i,j,k)              &
            & / EOSTable % DV % Variables(12) % Values(i,j,k) )              & 
            & - ( EOSTable % TS % States(3) % Values(k) * me ) )
        END DO
      END DO
    END DO


  DO l = 1, 15 !EOSTable % nVariables
    WRITE (*,*) EOSTable % DV % Names(l)
    minvar = MINVAL( EOSTable % DV % Variables(l) % Values )
    WRITE (*,*) "minvar=", minvar
    EOSTable % DV % Offsets(l) = -2.d0 * MIN( 0.d0, minvar )
    WRITE (*,*) "Offset=", EOSTable % DV % Offsets(l)
    EOSTable % DV % Variables(l) % Values &
      = LOG10( EOSTable % DV % Variables(l) % Values &
               + EOSTable % DV % Offsets(l) + epsilon )

  END DO

    CALL WriteEquationOfStateTableHDF( EOSTable )

    CALL DescribeEquationOfStateTable( EOSTable )

  END SUBROUTINE ReadSCTableHDF

  SUBROUTINE ReadChimeraProfile1D( FileName, MaxZone, r, rho, T, Ye, p, s, &
                                     e_internal, xn, xp, xhe, xa, chem_n,  &
                                     chem_p, chem_e, a_heavy, z_heavy,     &
                                     be_heavy, u, rstmss, vsound, gamma1,  &
                                     SkipLinesOption )

    CHARACTER(len=*), INTENT(in)  :: FileName
    INTEGER, INTENT(in)           :: MaxZone
    INTEGER, INTENT(in), OPTIONAL :: SkipLinesOption
    REAL(dp), DIMENSION(:), INTENT(out) :: r, rho, T, Ye, p, s, e_internal, xn, &
                                             xp, xhe, xa, chem_n, chem_p, chem_e, &
                                             a_heavy, z_heavy, be_heavy, u, &
                                             rstmss, vsound, gamma1  

    INTEGER :: i
    INTEGER :: j, FileUnit, SkipLines, istate
    LOGICAL, DIMENSION(MaxZone) :: L, shock
    CHARACTER(len=1), DIMENSION(MaxZone) :: nnse, EOS
    CHARACTER(LEN=255) :: stanza_name   ! Name of stanza headers
    REAL(dp), DIMENSION(MaxZone) :: v, w, dr, dummy14, flat, lum, &
                                      Tmev, ah, dum1, dum2, dum3, &
                                      dum4, dum5, dum6, dum7,     &
                                      dum8, dum9, dum10, dum11,   &
                                      dum12, dum13

    257 FORMAT (1x,i4,6es11.3,a1,es11.3,L3,es11.3,es14.6,4(es11.3),es10.2,a1,es11.3,1x,a1)
    258 FORMAT (1x,i4,10es12.4)
    259 FORMAT (1x,i4,8es11.4)
    260 FORMAT (1x,i4,11es11.4)
    261 FORMAT (1x,i4,14es11.4)
    262 FORMAT (1x,i4,10es11.4)
    263 FORMAT (1x,i4,18es11.4)
    264 FORMAT (1x,i4,7es11.4)

    SkipLines = 0
    IF ( Present ( SkipLinesOption ) ) &
      SkipLines = SkipLinesOption

    OPEN( newunit = FileUnit, FILE = FileName, status ="old", iostat=istate )

    DO i = 1, SkipLines
      READ(FileUnit,*)
    END DO

    DO
      IF ( istate /= 0 ) EXIT

      READ(FileUnit,'(a)',iostat=istate) stanza_name

      IF ( TRIM(ADJUSTL(stanza_name)) == 'Model configuration' ) THEN
        PRINT*, TRIM(ADJUSTL(stanza_name))
        DO i = 1, 4 
          READ(FileUnit,*)
        END DO

        DO i = MaxZone, 1, -1
          READ(FileUnit, 257) j, u(i), v(i), w(i), r(i), dr(i), dummy14(i), shock(i), &
                              flat(i), L(i), lum(i), rstmss(i), rho(i), Tmev(i),   &
                              T(i), s(i), ah(i), nnse(i), Ye(i), EOS(i)
        END DO
      END IF

      IF ( TRIM(ADJUSTL(stanza_name)) == 'Energy and shock ram pressure data' ) THEN
        PRINT*, TRIM(ADJUSTL(stanza_name))
        DO i = 1, 4 
          READ(FileUnit,*)
        END DO

        DO i = MaxZone, 1, -1
          READ(FileUnit, 258) j, dum1(i), dum2(i), dum3(i), dum4(i), p(i), &
                              dum5(i), dum6(i), dum7(i), dum8(i), dum9(i) 
        END DO
      END IF

      IF ( TRIM(ADJUSTL(stanza_name)) == 'EVH1 Energy Data' ) THEN
        PRINT*, TRIM(ADJUSTL(stanza_name))
        DO i = 1, 4 
          READ(FileUnit,*)
        END DO

        DO i = MaxZone, 1, -1
          READ(FileUnit, 259) j, e_internal(i), dum1(i), dum2(i), dum3(i), dum4(i), &
                              dum5(i), dum6(i), dum7(i)
        END DO
      END IF

      IF ( TRIM(ADJUSTL(stanza_name)) == 'Hydrodynamic data' ) THEN
        PRINT*, TRIM(ADJUSTL(stanza_name))
        DO i = 1, 4 
          READ(FileUnit,*)
        END DO

        DO i = MaxZone, 1, -1
          READ(FileUnit, 261) j, dum1(i), dum2(i), dum3(i), dum4(i), &
                              dum5(i), vsound(i), dum6(i), dum7(i), gamma1(i), &
                              dum9(i), dum10(i), dum11(i), dum12(i), dum13(i)
        END DO
      END IF

      IF ( TRIM(ADJUSTL(stanza_name)) == 'Chemical potential data' ) THEN
        PRINT*, TRIM(ADJUSTL(stanza_name))
        DO i = 1, 4 
          READ(FileUnit,*)
        END DO

        DO i = MaxZone, 1, -1
          READ(FileUnit, 262) j, chem_n(i), chem_p(i), dum1(i), chem_e(i), &
                              dum2(i), dum3(i), dum4(i), dum5(i), dum6(i), &
                              dum7(i)
        END DO
      END IF

      IF ( TRIM(ADJUSTL(stanza_name)) == 'Comp. & Mass Data' ) THEN
        PRINT*, TRIM(ADJUSTL(stanza_name))
        DO i = 1, 8 
          READ(FileUnit,*)
        END DO

        DO i = MaxZone, 1, -1
          READ(FileUnit, 263) j, xn(i), xp(i), xa(i), xhe(i), a_heavy(i), &
                              dum1(i), z_heavy(i), dum2(i), dum3(i), dum4(i), &
                              dum5(i), dum6(i), dum7(i), dum12(i), dum8(i),&
                              dum9(i), dum10(i), dum11(i)
        END DO
      END IF

      IF ( TRIM(ADJUSTL(stanza_name)) == 'Auxilary Nucleus Comp. Data and Burn Data' ) THEN
        PRINT*, TRIM(ADJUSTL(stanza_name))
        DO i = 1, 8 
          READ(FileUnit,*)
        END DO

        DO i = MaxZone, 1, -1
          READ(FileUnit, 262) j, dum1(i), dum2(i), be_heavy(i), dum3(i), dum4(i), &
                              dum5(i), dum6(i), dum7(i), dum8(i), dum9(i)
        END DO
        EXIT
      END IF
      
    END DO

    CLOSE(FileUnit)

  END SUBROUTINE ReadChimeraProfile1D

SUBROUTINE ReadChimeraThermoProfile1D( FileName, MaxZone, r, rho, T, Ye, s, &
                                     SkipLinesOption )

    CHARACTER(len=*), INTENT(in)  :: FileName
    INTEGER, INTENT(in)           :: MaxZone
    INTEGER, INTENT(in), OPTIONAL :: SkipLinesOption
    REAL(dp), DIMENSION(:), INTENT(out) :: r, rho, T, Ye, s
    INTEGER :: i, j, FileUnit, SkipLines, istate
    LOGICAL, DIMENSION(MaxZone) :: L, shock
    CHARACTER(len=1), DIMENSION(MaxZone) :: nnse, EOS
    CHARACTER(LEN=255) :: stanza_name   ! Name of stanza headers
    REAL(dp), DIMENSION(MaxZone) :: u, v, w, dr, dummy14, flat, lum, &
                                      Tmev, ah, dum1, dum2, dum3, &
                                      dum4, dum5, dum6, dum7,     &
                                      dum8, dum9, dum10, dum11,   &
                                      dum12, dum13, rstmss

    257 FORMAT (1x,i4,6es11.3,a1,es11.3,L3,es11.3,es14.6,4(es11.3),es10.2,a1,es11.3,1x,a1)
    258 FORMAT (1x,i4,10es12.4)
    259 FORMAT (1x,i4,8es11.4)
    260 FORMAT (1x,i4,11es11.4)
    261 FORMAT (1x,i4,14es11.4)
    262 FORMAT (1x,i4,10es11.4)
    263 FORMAT (1x,i4,18es11.4)
    264 FORMAT (1x,i4,7es11.4)

    SkipLines = 0
    IF ( Present ( SkipLinesOption ) ) &
      SkipLines = SkipLinesOption

    OPEN( newunit = FileUnit, FILE = FileName, status ="old", iostat=istate )

    DO i = 1, SkipLines
      READ(FileUnit,*)
    END DO

    DO
      !IF ( istate /= 0 ) EXIT
      READ(FileUnit,'(a)',iostat=istate) stanza_name
      IF ( TRIM(ADJUSTL(stanza_name)) == 'Model configuration' ) THEN
        PRINT*, TRIM(ADJUSTL(stanza_name))
        DO i = 1, 4
          READ(FileUnit,*)
        END DO

        DO i = MaxZone, 1, -1
          READ(FileUnit, 257) j, u(i), v(i), w(i), r(i), dr(i), dummy14(i), shock(i), &
                              flat(i), L(i), lum(i), rstmss(i), rho(i), Tmev(i),   &
                              T(i), s(i), ah(i), nnse(i), Ye(i), EOS(i)
        END DO
        EXIT
      END IF
    END DO

    CLOSE(FileUnit)

  END SUBROUTINE ReadChimeraThermoProfile1D

  SUBROUTINE ReadCHIMERAHDF( Rho, T, Ye, E_Int, Entropy, NSE, imax, nx, ny, &
                             nz, FileName)

    CHARACTER(len=*), INTENT(in)                :: FileName
    INTEGER, INTENT(out) :: imax, nx, ny, nz
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(out) :: Rho
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(out) :: T
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(out) :: Ye
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(out) :: E_Int
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(out) :: Entropy
    INTEGER, DIMENSION(:,:,:), ALLOCATABLE, INTENT(out) :: NSE
    INTEGER, DIMENSION(2) :: indices

    INTEGER(HSIZE_T), DIMENSION(1)              :: datasize1d
    INTEGER(HSIZE_T), DIMENSION(3)              :: datasize3d
    INTEGER(HID_T)                              :: file_id
    INTEGER(HID_T)                              :: group_id

    CALL OpenFileHDF( FileName, .false., file_id )

    CALL OpenGroupHDF( "mesh", .false., file_id, group_id )

    datasize1d(1) = 2
    CALL read_1d_slab_int('radial_index_bound', indices, group_id, &
           datasize1d)
    imax = indices(2)
    nx = imax + 2

    CALL read_1d_slab_int('theta_index_bound', indices, group_id, &
           datasize1d)
    ny = indices(2)

    CALL read_1d_slab_int('phi_index_bound', indices, group_id, &
           datasize1d)
    nz = indices(2)

    CALL CloseGroupHDF( group_id )

    ALLOCATE( Rho( nx, ny, nz ), T( nx, ny, nz ), Ye( nx, ny, nz ),           &
              E_Int( nx, ny, nz ), Entropy( nx, ny, nz ), NSE( nx + 1, ny, nz ) )

    CALL OpenGroupHDF( "fluid", .false., file_id, group_id )

    datasize3d = (/nx,ny,nz/)
    CALL ReadHDF( "rho_c", Rho(:,:,:), &
                              group_id, datasize3d )

    CALL ReadHDF( "t_c", T(:,:,:), &
                              group_id, datasize3d )

    CALL ReadHDF( "ye_c", Ye(:,:,:), &
                              group_id, datasize3d )

    CALL ReadHDF( "e_int", E_Int(:,:,:), &
                              group_id, datasize3d )

    CALL ReadHDF( "entropy", Entropy(:,:,:), &
                              group_id, datasize3d )

    CALL CloseGroupHDF( group_id )

    CALL OpenGroupHDF( "abundance", .false., file_id, group_id )

    datasize3d = SHAPE(NSE)
    CALL ReadHDF( "nse_c", NSE(:,:,:), &
                              group_id, datasize3d )

    CALL CloseGroupHDF( group_id )

  END SUBROUTINE ReadCHIMERAHDF

END MODULE wlIOModuleCHIMERA
