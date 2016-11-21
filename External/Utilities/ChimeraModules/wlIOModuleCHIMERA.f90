MODULE wlIOModuleCHIMERA

  USE wlKindModule, ONLY: dp
  USE HDF5
  USE wlThermoStateModule
  USE wlDependentVariablesModule
  USE wlEquationOfStateTableModule
!  USE wlInterpolationModule
  USE wlExtNumericalModule, ONLY: epsilon, zero
  USE wlExtPhysicalConstantsModule
  USE e_p_eos_module
  USE EL_EOS_MODULE
  USE EOS_M4C_MODULE
  USE MAXWEL_MODULE
  USE EOS_BCK_MODULE
  USE wlExtEOSWrapperModule, ONLY: wlGetElectronEOS
  USE wlIOModuleHDF
  USE wlEOSIOModuleHDF
  !USE sfhx_frdm_composition_module

  implicit none

  PUBLIC ReadChimeraProfile1D 
  PUBLIC ReadCHIMERAHDF
  PUBLIC ReadComposeTableHDF

CONTAINS

  SUBROUTINE ReadComposeTableHDF( EOSTable, FileName )
    TYPE(EquationOfStateTableType), INTENT(inout) :: EOSTable
    CHARACTER(len=*), INTENT(in)                  :: FileName

    INTEGER, DIMENSION(3)                         :: nPoints
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
    !REAL(dp), DIMENSION(0:325)                    :: nb
    !REAL(dp), DIMENSION(0:80)                     :: t
    !REAL(dp), DIMENSION(0:59)                     :: yq
    !REAL(dp), DIMENSION(0:14259239)               :: thermo !  nPoints(1) x 2 x 3 x # of thermo quantities, nThermo(1) 
    !REAL(dp), DIMENSION(0:4753079)                :: yi     ! ! nPoints(1) x 2 x 3 x 3
    !REAL(dp), DIMENSION(0:1584359)                :: yav ! nPoints(1) x 2 x 3
    !REAL(dp), DIMENSION(0:1584359)                :: aav ! nPoints(1) x 2 x 3
    !REAL(dp), DIMENSION(0:1584359)                :: zav ! nPoints(1) x 2 x 3
    INTEGER                                       :: nVariables
    INTEGER                                       :: i, j, k, l
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
    REAL(dp)                                      :: minvar
    REAL(dp), DIMENSION(8140)           :: mass, bind, degen
    REAL(dp), DIMENSION(81,0:60,326,8)           :: comp
    INTEGER, DIMENSION(8140,2) :: az

    CALL OpenFileHDF( FileName, .false., file_id )

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

    CALL AllocateEquationOfStateTable( EOSTable, nPoints , nVariables )
write (*,*) 'EOSTable allocated'

    ALLOCATE( nb(0:(nPoints(1) - 1 ) ), t(0:(nPoints(2) - 1) ), yq( 0:(nPoints(3) - 1 ) ) )
    ALLOCATE( thermo(0:(nThermo(1)*AllPoints(1) - 1)) )
    ALLOCATE( yi(0:(3*AllPoints(1) - 1)), yav(0:(AllPoints(1) - 1)), aav(0:(AllPoints(1) - 1)), zav( 0:(AllPoints(1) - 1) ) )
    
write (*,*) 'Buffers Allocated'
    datasize1d(1) = SIZE(nb)
    CALL ReadHDF( "nb", nb(:), file_id, datasize1d )

    datasize1d(1) = SIZE(t)
    CALL ReadHDF( "t", t(:), file_id, datasize1d )

    datasize1d(1) = SIZE(yq)
    CALL ReadHDF( "yq", yq(:), file_id, datasize1d )

write (*,*) 'Buffers read'

    DO i = 0, nPoints(1) - 1
      EOSTable % TS % States(1) % Values(i+1) = nb(i) / kfm
    END DO
write (*,*) 'rho', EOSTable % TS % States(1) % Values 

    DO i = 0, nPoints(2) - 1
      EOSTable % TS % States(2) % Values(i+1) = t(i) / kmev
    END DO
write (*,*) 'T'

    DO i = 0, nPoints(3) - 1
      EOSTable % TS % States(3) % Values(i+1) = yq(i) 
    END DO
write (*,*) 'ye'
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
            != ergmev * mn * ( thermo((i + nPoints(1)*(j-1) + TwoPoints(1)*(k-1)) + 4*AllPoints(1) ) + 1 ) * avn ! need to multiply by baryons per gram, avn 
            = ergmev * mn * ( thermo((i + nPoints(1)*(j-1) + TwoPoints(1)*(k-1)) + 4*AllPoints(1) ) + 8.9d0/mn ) * avn ! need to multiply by baryons per gram, avn 
        END DO
      END DO
    END DO
write (*,*) 'energy(1,1,1)', EOSTable % DV % Variables(3) % Values(1,1,1)

    DO k = 1, nPoints(3)
      DO j = 1, nPoints(2)
        DO i = 0, nPoints(1) - 1
!          EOSTable % DV % Variables(4) % Values(i+1,j,k) &
!            = thermo((i + nPoints(1)*(j-1) + TwoPoints(1)*(k-1)) + 4*AllPoints(1) ) * mn
          CALL wlGetElectronEOS( EOSTable % TS % States(1) % Values(i+1), & 
                                 EOSTable % TS % States(2) % Values(j),   &
                                 EOSTable % TS % States(3) % Values(k),   &
                                 press_e, &
                                 entrop_e,&
                                 energ_e, &
                                 chem_e )

! Add corrections to pressure, internal energy, entropy)  

          press_buff = EOSTable % DV % Variables(1) % Values(i+1,j,k)
          EOSTable % DV % Variables(1) % Values(i+1,j,k) &
            = press_buff + press_e

          entrop_buff = EOSTable % DV % Variables(2) % Values(i+1,j,k)
          EOSTable % DV % Variables(2) % Values(i+1,j,k) &
            = entrop_buff + entrop_e

          energ_buff = EOSTable % DV % Variables(3) % Values(i+1,j,k) 
          EOSTable % DV % Variables(3) % Values(i+1,j,k) &
            = energ_buff + energ_e

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
            != ( thermo((i + nPoints(1)*(j-1) + TwoPoints(1)*(k-1)) + 2*AllPoints(1) ) + 1 ) * mn/1000d0 - dmnp
            = ( thermo((i + nPoints(1)*(j-1) + TwoPoints(1)*(k-1)) + 2*AllPoints(1) ) ) - dmnp
        END DO
      END DO
    END DO
write (*,*) 'neutron chem pot(1,1,1)', EOSTable % DV % Variables(6) % Values(1,1,1)

    DO k = 1, nPoints(3)
      DO j = 1, nPoints(2)
        DO i = 0, nPoints(1) - 1
          EOSTable % DV % Variables(5) % Values(i+1,j,k) &
            != thermo((i + nPoints(1)*(j-1) + TwoPoints(1)*(k-1)) + 3*AllPoints(1) ) * mn/1000d0 &
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
          EOSTable % DV % Variables(8) % Values(i+1,j,k) = yi(i + nPoints(1)*(j-1) + TwoPoints(1)*(k-1))
        END DO
      END DO
    END DO
write (*,*) 'neutron mass frac(1,1,1)', EOSTable % DV % Variables(8) % Values(1,1,1)

    DO k = 1, nPoints(3)
      DO j = 1, nPoints(2)
        DO i = 0, nPoints(1) - 1
          EOSTable % DV % Variables(7) % Values(i+1,j,k) &
            = yi((i + nPoints(1)*(j-1) + TwoPoints(1)*(k-1)) + AllPoints(1) )
        END DO
      END DO
    END DO
write (*,*) 'proton mass frac(1,1,1)', EOSTable % DV % Variables(7) % Values(1,1,1)

    DO k = 1, nPoints(3)
      DO j = 1, nPoints(2)
        DO i = 0, nPoints(1) - 1
          EOSTable % DV % Variables(9) % Values(i+1,j,k) &
            = yi((i + nPoints(1)*(j-1) + TwoPoints(1)*(k-1)) + 2*AllPoints(1) )
        END DO
      END DO
    END DO
write (*,*) 'alpha mass frac(1,1,1)', EOSTable % DV % Variables(9) % Values(1,1,1)
write(*,*), "yi DV's filled"

    datasize1d(1) = SIZE(yav)
    CALL ReadHDF( "yav", yav(:), file_id, datasize1d )
write(*,*), "yav read"
    DO k = 1, nPoints(3)
      DO j = 1, nPoints(2)
        DO i = 0, nPoints(1) - 1
          EOSTable % DV % Variables(10) % Values(i+1,j,k) = yav(i + nPoints(1)*(j-1) + TwoPoints(1)*(k-1))
        END DO
      END DO
    END DO
write (*,*) 'heavy mass frac(1,1,1)', EOSTable % DV % Variables(10) % Values(1,1,1)
write(*,*), "yav DV filled"

    datasize1d(1) = SIZE(zav)
    CALL ReadHDF( "zav", zav(:), file_id, datasize1d )
write(*,*), "zav read"
    DO k = 1, nPoints(3)
      DO j = 1, nPoints(2)
        DO i = 0, nPoints(1) - 1
          EOSTable % DV % Variables(11) % Values(i+1,j,k) = zav(i + nPoints(1)*(j-1) + TwoPoints(1)*(k-1))
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
          EOSTable % DV % Variables(12) % Values(i+1,j,k) = aav(i + nPoints(1)*(j-1) + TwoPoints(1)*(k-1))
        END DO
      END DO
    END DO
write (*,*) 'heavy mass number(1,1,1)', EOSTable % DV % Variables(12) % Values(1,1,1)
write(*,*), "aav DV filled"

    DO k = 1, nPoints(3)
      DO j = 1, nPoints(2)
        DO i = 1, nPoints(1)
          EOSTable % DV % Variables(13) % Values(i,j,k) = 8.755831d0
        END DO
      END DO
    END DO

    DO k = 1, nPoints(3)
      DO j = 1, nPoints(2)
        DO i = 1, nPoints(1)

     !REAL(dp), DIMENSION(8140)           :: mass, bind, degen
     !REAL(dp), DIMENSION(81,0:60,326,8)           :: comp
     !INTEGER, DIMENSION(8140,2) :: az

        !open(10,file='sfhx_frdm_comp_v1.03.bin',form='unformatted',
     !& status='old')
        !read(10) az,mass,bind,comp
        !write(10) az,mass,bind,comp
        !close(10)

        !return

          !EOSTable % DV % Variables(14) % Values(i,j,k) = ku * ( UTOT - EU + ee             &
!&            + dmnp * x_proton + 7.075 * x_alpha - BUNUC + 1.5d0 * tmev * XH/A - ye * me )

          EOSTable % DV % Variables(14) % Values(i,j,k) =                    &
            & EOSTable % DV % Variables(3) % Values(i,j,k)                   &
            & + ku * ( dmnp * EOSTable % DV % Variables(7) % Values(i,j,k)   &
            & + 7.075 * EOSTable % DV % Variables(9) % Values(i,j,k)         &
            & - EOSTable % DV % Variables(13) % Values(i,j,k)                &  
            & + 1.5d0 * EOSTable % TS % States(2) % Values(j)                &
            & * ( EOSTable % DV % Variables(10) % Values(i,j,k)              &
            & / EOSTable % DV % Variables(12) % Values(i,j,k) )              & 
            & - ( EOSTable % TS % States(3) % Values(k) * me ) )
        END DO
      END DO
    END DO

    DO k = 1, nPoints(3)
      DO j = 1, nPoints(2)
        DO i = 0, nPoints(1) - 1
          EOSTable % DV % Variables(15) % Values(i+1,j,k) &
            = thermo((i + nPoints(1)*(j-1) + TwoPoints(1)*(k-1)) + 5*AllPoints(1) ) ! *(neutron mass for units)
        END DO
      END DO
    END DO
write (*,*) 'gamma(1,1,1)', EOSTable % DV % Variables(15) % Values(1,1,1)

  DO l = 1, EOSTable % nVariables
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

  SUBROUTINE ReadChimeraProfile1D( FileName, MaxZone, r, rho, T, Ye, p, s, &
                                     e_internal, xn, xp, xhe, xa, chem_n,  &
                                     chem_p, chem_e, a_heavy, z_heavy,     &
                                     be_heavy, u, rstmss, vsound,          &
                                     SkipLinesOption )

    CHARACTER(len=*), INTENT(in)  :: FileName
    INTEGER, INTENT(in)           :: MaxZone
    INTEGER, INTENT(in), OPTIONAL :: SkipLinesOption
    REAL(dp), DIMENSION(:), INTENT(out) :: r, rho, T, Ye, p, s, e_internal, xn, &
                                             xp, xhe, xa, chem_n, chem_p, chem_e, &
                                             a_heavy, z_heavy, be_heavy, u, &
                                             rstmss, vsound  

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
                              dum5(i), vsound(i), dum6(i), dum7(i), dum8(i), &
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
      END IF

    END DO

    CLOSE(FileUnit)

  END SUBROUTINE ReadChimeraProfile1D

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
