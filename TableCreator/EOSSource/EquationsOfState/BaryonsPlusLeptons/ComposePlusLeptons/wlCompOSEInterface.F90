MODULE wlCompOSEInterface
    
    USE wlKindModule, ONLY: dp
    USE wlEosConstantsModule, ONLY: cvel, ergmev, cm3fm3, kmev_inv, rmu, mn, me, mp
    USE HDF5
    USE wlIOModuleHDF
    
    IMPLICIT NONE
    PRIVATE
    
    PUBLIC :: ReadnPointsFromCompOSE, ReadCompOSETable, ReadCompOSEHDFTable
    
    REAL(dp), PUBLIC, DIMENSION(:,:,:,:), ALLOCATABLE :: EOSCompOSE
    
    INTEGER(KIND=4) :: h5err !< hdf5 error code                                                                                                    
    ! output baryon density (fm^-3), temperature (kelvin), and proton (or charge) fraction	
    REAL(dp),         DIMENSION(:), ALLOCATABLE :: nbCompOSE
    REAL(dp), PUBLIC, DIMENSION(:), ALLOCATABLE :: RhoCompOSE
    REAL(dp), PUBLIC, DIMENSION(:), ALLOCATABLE :: TempCompOSE
    REAL(dp), PUBLIC, DIMENSION(:), ALLOCATABLE :: YpCompOSE
    
    INTEGER, PARAMETER :: nVariablesCompOSE             = 17
    INTEGER, PARAMETER :: iPressCompOSE                 = 1
    INTEGER, PARAMETER :: iEntropyCompOSE               = 2
    INTEGER, PARAMETER :: iInternalEnergyDensityCompOSE = 3
    INTEGER, PARAMETER :: iElectronChemPotCompOSE       = 4
    INTEGER, PARAMETER :: iProtonChemPotCompOSE         = 5
    INTEGER, PARAMETER :: iNeutronChemPotCompOSE        = 6
    INTEGER, PARAMETER :: iProtonMassFractionCompOSE    = 7
    INTEGER, PARAMETER :: iNeutronMassFractionCompOSE   = 8
    INTEGER, PARAMETER :: iAlphaMassFractionCompOSE     = 9
    INTEGER, PARAMETER :: iHeavyMassFractionCompOSE     = 10
    INTEGER, PARAMETER :: iHeavyChargeNumberCompOSE     = 11
    INTEGER, PARAMETER :: iHeavyMassNumberCompOSE       = 12
    INTEGER, PARAMETER :: iHeavyBindingEnergyCompOSE    = 13
    INTEGER, PARAMETER :: iThermalEnergyCompOSE         = 14
    INTEGER, PARAMETER :: iGammaCompose                 = 15
    INTEGER, PARAMETER :: iProtonSelfEnergyCompOSE      = 16
    INTEGER, PARAMETER :: iNeutronSelfEnergyCompOSE     = 17
    
    ! These are used for the HDF5 thermo table (i.e. Table 7.1 from Compose manual v. 3.0)
    INTEGER, PARAMETER :: iThermoPressure = 1
    INTEGER, PARAMETER :: iThermoEntropy  = 2
    INTEGER, PARAMETER :: iThermoMub      = 3
    INTEGER, PARAMETER :: iThermoMuq      = 4
    INTEGER, PARAMETER :: iThermoMul      = 5
    INTEGER, PARAMETER :: iThermoFreeene  = 6
    INTEGER, PARAMETER :: iThermoEps      = 7
    INTEGER, PARAMETER :: iThermoEnth     = 8
    INTEGER, PARAMETER :: iThermoFreeenth = 9
    INTEGER, PARAMETER :: iThermoDpdnb    = 10
    INTEGER, PARAMETER :: iThermoDpdeps   = 11
    INTEGER, PARAMETER :: iThermoCs2      = 12
    INTEGER, PARAMETER :: iThermoCv       = 13
    INTEGER, PARAMETER :: iThermoCp       = 14
    INTEGER, PARAMETER :: iThermoGamma    = 15
    INTEGER, PARAMETER :: iThermoAlphaP   = 16
    INTEGER, PARAMETER :: iThermoBetaV    = 17
    INTEGER, PARAMETER :: iThermoKappaT   = 18
    INTEGER, PARAMETER :: iThermoKappaS   = 19
    
    ! These are indices for the specific EOS at hand. Might have to be changed depending on the EOS you are 
    ! reading. This is for the SFHo. Notice that technically CompOSE provides a recipe to assign indices,
    ! but be careful becasue not all the tables follow that recipe, so double check each time, this is why we
    ! hard code the indices. Information can be found in the eos.pdf file from compose
    INTEGER, PARAMETER :: iCompoNeutron = 10
    INTEGER, PARAMETER :: iCompoProton = 11
    INTEGER, PARAMETER :: iCompoHe4 = 4002
    INTEGER, PARAMETER :: iCompoHe3 = 3002
    INTEGER, PARAMETER :: iCompoH3 = 3001
    INTEGER, PARAMETER :: iCompoH2 = 2001
    INTEGER, PARAMETER :: iCompoHeavy = 999
    INTEGER, PARAMETER :: iMicroNeutSelfEnergy = 10051
    INTEGER, PARAMETER :: iMicroProtSelfEnergy = 11051

    CONTAINS
    
    SUBROUTINE ReadnPointsFromCompOSE( CompOSEFilePath, nRho, nTemp, nYp )
        
        ! input path of the directory with the CompOSE output
        CHARACTER(len=128), INTENT(IN) :: CompOSEFilePath
        
        ! Output variables
        INTEGER, INTENT(INOUT) :: nRho, nTemp, nYp
        
        ! get density grid
        OPEN(123, FILE=TRIM(ADJUSTL(CompOSEFilePath))//'/eos.nb', & 
        FORM = "formatted", ACTION = 'read')
        READ(123,*) 
        READ(123,*) nRho
        
        CLOSE(123)
        
        ! get temperature grid
        OPEN(123, FILE=TRIM(ADJUSTL(CompOSEFilePath))//'/eos.t', & 
        FORM = "formatted", ACTION = 'read')
        READ(123,*)
        READ(123,*) nTemp
        
        CLOSE(123)
        
        ! get yp grid
        OPEN(123, FILE=TRIM(ADJUSTL(CompOSEFilePath))//'/eos.yq', & 
        FORM = "formatted", ACTION = 'read')
        READ(123,*)
        READ(123,*) nYp
        
        CLOSE(123)
        
    END SUBROUTINE ReadnPointsFromCompOSE
    
    SUBROUTINE ReadCompOSETable( CompOSEFilePath, nRho, nTemp, nYp, &
          AllocateEOS_Optional )
        
        ! input path of the directory with the CompOSE output
        CHARACTER(len=128), INTENT(IN) :: CompOSEFilePath
        
        ! output number of points in the table
        INTEGER, INTENT(IN) :: nRho, nTemp, nYp
        
        LOGICAL, INTENT(IN), OPTIONAL :: AllocateEOS_Optional

        ! Local Variables
        LOGICAL  :: AllocateEOS
        INTEGER  :: iLepton, iRho, iT, iYp, i_tot
        REAL(dp) :: NeutronMass, ProtonMass
        REAL(dp) :: Q1, Q2, Q3, Q4, Q5, Q6, Q7, TotalInternalEnergy
        REAL(dp) :: Xp, Xn, Xa, Xh, Abar, Zbar, neut_self_ene, prot_self_ene
        
        IF ( PRESENT(AllocateEOS_Optional) ) THEN
          AllocateEOS = AllocateEOS_Optional
        ELSE
          AllocateEOS = .TRUE.
        END IF

        ! Get proton and neutron mass
        OPEN(123, FILE=TRIM(ADJUSTL(CompOSEFilePath))//'/eos.thermo', & 
        FORM = "formatted", ACTION = 'read')
        READ(123,*) NeutronMass, ProtonMass, iLepton
        CLOSE(123)
        
        ! Uncomment if you DO NOT want to read the lepton component
        IF (iLepton .ne. 0.0_dp) THEN
            WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!! CAREFUL !!!!!!!!!!!!!!!!!!!!!!!!!!'
            WRITE(*,*) 'The table you are trying to read already contains leptons!'
            WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!! CAREFUL !!!!!!!!!!!!!!!!!!!!!!!!!!'
        END IF
        
        ! Get density. DECIDE HOW TO CALCULATE THE MASS DENSITY SINCE IT IS NOT OBVIOUS,
        ! Here we use the a.m.u.
        IF (AllocateEOS) THEN
          ALLOCATE( nbCompOSE( nRho ) )
          ALLOCATE( RhoCompOSE( nRho ) )
        ENDIF
        
        OPEN(123, FILE=TRIM(ADJUSTL(CompOSEFilePath))//'/eos.nb', & 
        FORM = "formatted", ACTION = 'read')
        READ(123,*)
        READ(123,*)
        DO iRho=1,nRho
            READ(123,*) nbCompOSE(iRho)
        END DO
        CLOSE(123)	
        ! Convert to grams per cm^3
        RhoCompOSE(:) = nbCompOSE(:)*rmu/cm3fm3
        
        ! get temperature
        OPEN(123, FILE=TRIM(ADJUSTL(CompOSEFilePath))//'/eos.t', & 
        FORM = "formatted", ACTION = 'read')
        READ(123,*)
        READ(123,*)
        
        IF (AllocateEOS) THEN
          ALLOCATE( TempCompOSE( nTemp ) )
        ENDIF
        
        DO iT=1,nTemp
            READ(123,*) TempCompOSE(iT)
        END DO
        
        CLOSE(123)
        ! Convert to kelvin
        TempCompOSE(:) = TempCompOSE(:)*kmev_inv
        
        ! get yp grid
        OPEN(123, FILE=TRIM(ADJUSTL(CompOSEFilePath))//'/eos.yq', & 
        FORM = "formatted", ACTION = 'read')
        READ(123,*)
        READ(123,*)
        
        IF (AllocateEOS) THEN
          ALLOCATE( YpCompOSE( nYp ) )
        ENDIF

        DO iYp=1,nYp
            READ(123,*) YpCompOSE(iYp)
        END DO
        
        CLOSE(123)
        
        ! Allocate table
        IF (AllocateEOS) THEN
          ALLOCATE( EOSCompOSE(nRho,nTemp,nYp,nVariablesCompOSE) )
        ENDIF
        
        ! get thermal state
        OPEN(123, FILE=TRIM(ADJUSTL(CompOSEFilePath))//'/eos.thermo', & 
        FORM = "formatted", ACTION = 'read')
        READ(123,*) NeutronMass, ProtonMass, iLepton
        
        ! Read CompOSE EOS and convert to weaklib units
        ! for the internal energy you have to be careful. You have to convert baryon density to mass density. It seems to me that
        ! the best way of doing this is to set a fixed conversion factor, i.e. the nucleon mass, i.e. the a.m.u.
        ! Notice that Q7 technically uses the neutron mass, so be careful with that
        DO i_tot=1,nTemp*nRho*nYp
            
            READ(123,*) iT, iRho, iYp, Q1, Q2, Q3, Q4, Q5, Q6, Q7
            
            EOSCompOSE(iRho,iT,iYp,iPressCompOSE) = Q1 * nbCompOSE(iRho) * ergmev / cm3fm3
            EOSCompOSE(iRho,iT,iYp,iEntropyCompOSE) = Q2
            ! This is how I used to do it
            ! TotalInternalEnergy = (1.0_dp + Q7) * NeutronMass
            ! EOSCompOSE(iRho,iT,iYp,iInternalEnergyDensityCompOSE) = &
            ! TotalInternalEnergy*nbCompOSE(iRho)/cm3fm3 * ergmev / RhoCompOSE(iRho) - cvel**2.0d0
            
            ! This is how it should be done I think
            !EOSCompOSE(iRho,iT,iYp,iInternalEnergyDensityCompOSE) = Q7 * ergmev / rmu
            EOSCompOSE(iRho,iT,iYp,iInternalEnergyDensityCompOSE) = &
              (Q7 + 1.0_dp) * NeutronMass * ergmev / rmu - cvel**2.0_dp

            ! now handle the chemical potentials. This is the relativistic expression including the 
            ! nucleon mass. Be careful because some opacities need the NR version, or shifted to a reference
            ! mass, e.g. both mup and mun use the neutron mass...
            EOSCompOSE(iRho,iT,iYp,iNeutronChemPotCompOSE) = (Q3 + 1.0_dp)*NeutronMass + NeutronMass
            EOSCompOSE(iRho,iT,iYp,iProtonChemPotCompOSE) = (Q3 + 1.0_dp)*NeutronMass + Q4*NeutronMass + ProtonMass

            IF (Q4 == 0.0d0) THEN
                WRITE(*,*) iT, iRho, iYp
            ENDIF
            
        END DO
        
        CLOSE(123)
        
        ! get composition and microscopic quantities from custom made table 
        ! (see python script to convert eos.compo into eos.compomicro.wl
        OPEN(123, FILE=TRIM(ADJUSTL(CompOSEFilePath))//'/eos.compomicro.wl', & 
        FORM = "formatted", ACTION = 'read')
        
        ! Here you prrobably need to decide what to DO with the light nuclei
        DO i_tot=1,nTemp*nRho*nYp
            
            READ(123,*) iT, iRho, iYp, Xp, Xn, Xa, Xh, Abar, Zbar, neut_self_ene, prot_self_ene
            
            EOSCompOSE(iRho,iT,iYp,iProtonMassFractionCompOSE)  = Xp
            EOSCompOSE(iRho,iT,iYp,iNeutronMassFractionCompOSE) = Xn
            EOSCompOSE(iRho,iT,iYp,iAlphaMassFractionCompOSE)   = Xa
            EOSCompOSE(iRho,iT,iYp,iHeavyMassFractionCompOSE)   = Xh
            EOSCompOSE(iRho,iT,iYp,iHeavyMassNumberCompOSE)     = Abar
            EOSCompOSE(iRho,iT,iYp,iHeavyChargeNumberCompOSE)   = Zbar
            EOSCompOSE(iRho,iT,iYp,iNeutronSelfEnergyCompOSE)   = neut_self_ene
            EOSCompOSE(iRho,iT,iYp,iProtonSelfEnergyCompOSE)    = prot_self_ene
            
            ! Finally set other things to zero for now, CompOSE does not tell you what they are
            EOSCompOSE(iRho,iT,iYp,iHeavyBindingEnergyCompOSE) = 0.0_dp
            EOSCompOSE(iRho,iT,iYp,iThermalEnergyCompOSE)      = 0.0_dp
            
        END DO
        
        CLOSE(123)

    END SUBROUTINE ReadCompOSETable
    
    SUBROUTINE ReadCompOSEHDFTable( FileName, nRho, nTemp, nYp, iLepton, &
          AllocateEOS_Optional )
        
        ! HDF5 variables
        INTEGER(HID_T) :: file_id
        INTEGER(HID_T) :: group_id
        INTEGER(HSIZE_T) :: datasize4d(4)
        INTEGER(HSIZE_T) :: datasize3d(3)
        INTEGER(HSIZE_T) :: datasize1d(1)
        REAL(dp), DIMENSION(1) :: buffer
        
        ! input path of the directory with the CompOSE output
        CHARACTER(len=128), INTENT(IN) :: FileName
        
        ! output number of points in the table
        INTEGER, INTENT(IN)    :: iLepton 
        INTEGER, INTENT(INOUT) :: nRho, nTemp, nYp

        LOGICAL, INTENT(IN), OPTIONAL :: AllocateEOS_Optional
        
        ! Local Variables
        LOGICAL  :: AllocateEOS
        INTEGER  :: nThermo, nPairs, nQuads, nMicro
        INTEGER  :: iThermo, iCompo, iMicro, iRho, iT, iYp
        REAL(dp) :: NeutronMass, ProtonMass
        REAL(dp) :: TotalInternalEnergy
        REAL(dp), ALLOCATABLE, DIMENSION(:,:,:,:) :: ThermoTable, CompoTable, MicroTable
        REAL(dp), ALLOCATABLE, DIMENSION(:,:,:,:) :: ZbarTable, AbarTable, YhTable

        INTEGER, ALLOCATABLE, DIMENSION(:) :: ThermoIndices, CompoIndices, MicroIndices
        
        INTEGER :: iPress_in_Table, iEntr_in_Table, iEps_in_Table, iMub_in_Table, iMuq_in_Table, &
                   iXp_in_Table, iXn_in_Table, iXa_in_Table, &
                   iNeutSelfEnergy_in_Table, iProtSelfEnergy_in_Table
        
        LOGICAL :: fix_composition
        
        REAL(DP) :: NewAbarNumerator, NewZbarNumerator, NewHeavyDenominator, NewFracHeavy, &
            Light_Abar, Light_Zbar
        
        IF ( PRESENT(AllocateEOS_Optional) ) THEN
          AllocateEOS = AllocateEOS_Optional
        ELSE
          AllocateEOS = .TRUE.
        END IF

        NeutronMass = mn
        ProtonMass = mp
        
        CALL OpenFileHDF( FileName, .false., file_id )
        
        CALL OpenGroupHDF( "Parameters", .false., file_id, group_id )
        CALL hdf5_read_int4_attr0D(group_id, "pointsnb", nRho)
        CALL hdf5_read_int4_attr0D(group_id, "pointst", nTemp)
        CALL hdf5_read_int4_attr0D(group_id, "pointsyq", nYp)
        
        ALLOCATE( nbCompOSE(nRho) )
        IF (AllocateEOS) THEN
          ALLOCATE( RhoCompOSE(nRho) )
          ALLOCATE( TempCompOSE(nTemp) )
          ALLOCATE( YpCompOSE(nYp) )
        ENDIF
        
        datasize1d = nRho
        CALL ReadHDF( "nb", nbCompOSE(:), group_id, datasize1d )
        datasize1d = nTemp
        CALL ReadHDF( "t", TempCompOSE(:), group_id, datasize1d )		
        datasize1d = nYp
        CALL ReadHDF( "yq", YpCompOSE(:), group_id, datasize1d )
        CALL CloseGroupHDF( group_id )
        
        ! Read thermo quantities
        CALL OpenGroupHDF( "Thermo_qty", .false., file_id, group_id )
        datasize1d(1) = 1
        CALL hdf5_read_int4_attr0D(group_id, "pointsqty", nThermo)
        
        ALLOCATE( ThermoTable(nRho, nTemp, nYp, nThermo) )
        ALLOCATE( ThermoIndices(nThermo) )
        
        datasize4d = (/ nRho, nTemp, nYp, nThermo /)
        CALL ReadHDF( "thermo", ThermoTable, group_id, datasize4d )
        datasize1d = nThermo
        CALL ReadHDF( "index_thermo", ThermoIndices(:), group_id, datasize1d )
        CALL CloseGroupHDF( group_id )
        
        ! Read pairs composition data
        CALL OpenGroupHDF( "Composition_pairs", .false., file_id, group_id )
        CALL hdf5_read_int4_attr0D(group_id, "pointspairs", nPairs)
        
        ALLOCATE( CompoTable(nRho, nTemp, nYp, nPairs) )
        ALLOCATE( CompoIndices(nPairs) )
        
        datasize4d = (/ nRho, nTemp, nYp, nPairs /) 
        CALL ReadHDF( "yi", CompoTable, group_id, datasize4d )		
        datasize1d = nPairs
        CALL ReadHDF( "index_yi", CompoIndices(:), group_id, datasize1d )
        CALL CloseGroupHDF( group_id )
        
        ! Read quadruples (i.e. heavy nucleus) composition data
        CALL OpenGroupHDF( "Composition_quadrupels", .false., file_id, group_id )
        datasize1d(1) = 1
        CALL hdf5_read_int4_attr0D(group_id, "pointsav", nQuads)
        
        ALLOCATE( ZbarTable(nRho, nTemp, nYp, 1) )
        ALLOCATE( AbarTable(nRho, nTemp, nYp, 1) )
        ALLOCATE( YhTable(nRho, nTemp, nYp, 1) )
        
        datasize4d = (/ nRho, nTemp, nYp, 1 /) 
        CALL ReadHDF( "yav", YhTable, group_id, datasize4d )
        CALL ReadHDF( "aav", AbarTable, group_id, datasize4d )
        CALL ReadHDF( "zav", ZbarTable, group_id, datasize4d )
        CALL CloseGroupHDF( group_id )
        
        IF (nQuads .gt. 1) THEN
            WRITE(*,*) 'Cannot handle more than one quadruple, i.e. the average heavy nucleus'
        END IF
        
        ! Read microscopic quantities
        CALL OpenGroupHDF( "Micro_qty", .false., file_id, group_id )
        datasize1d(1) = 1
        CALL hdf5_read_int4_attr0D(group_id, "pointsmicro", nMicro)
        
        ALLOCATE( MicroTable(nRho, nTemp, nYp, nMicro) )
        ALLOCATE( MicroIndices(nMicro) )
        
        datasize4d = (/ nRho, nTemp, nYp, nMicro /)
        CALL ReadHDF( "micro", MicroTable, group_id, datasize4d )
        datasize1d = nMicro
        CALL ReadHDF( "index_micro", MicroIndices(:), group_id, datasize1d )
        CALL CloseGroupHDF( group_id )
        CALL CloseFileHDF( file_id )

        ! Convert to grams per cm^3
        RhoCompOSE(:) = nbCompOSE(:)*rmu/cm3fm3
        ! Convert to kelvin
        TempCompOSE(:) = TempCompOSE(:)*kmev_inv
        
        ! Find thermal quantities you need
        iPress_in_Table = 0
        iEntr_in_Table = 0
        iEps_in_Table = 0
        iMub_in_Table = 0
        iMuq_in_Table = 0
        
        ! Find the indices you want in the table
        DO iThermo=1,nThermo
            IF (ThermoIndices(iThermo) == iThermoPressure) THEN
                iPress_in_Table = iThermo
            ELSE IF (ThermoIndices(iThermo) == iThermoEntropy) THEN
                iEntr_in_Table = iThermo
            ELSE IF (ThermoIndices(iThermo) == iThermoEps) THEN
                iEps_in_Table = iThermo
            ELSE IF (ThermoIndices(iThermo) == iThermoMub) THEN
                iMub_in_Table = iThermo
            ELSE IF (ThermoIndices(iThermo) == iThermoMuq) THEN
                iMuq_in_Table = iThermo
            END IF
        END DO		
        
        ! Handle exceptions
        IF ( (iPress_in_Table == 0) .OR. (iEntr_in_Table == 0) .OR. (iEps_in_Table == 0) &
        .OR. (iMub_in_Table == 0) .OR. (iMuq_in_Table == 0) ) THEN
            WRITE(*,*) 'Not enough quantities in the CompOSE table', iPress_in_Table, iEntr_in_Table, &
            iEps_in_Table, iMub_in_Table, iMuq_in_Table
        ENDIF
        
        ! Find microscopic quantities you need
        iNeutSelfEnergy_in_Table = 0
        iProtSelfEnergy_in_Table = 0
        ! Find the indices you want in the table
        DO iMicro=1,nMicro
          WRITE(*,*) iMicro, MicroIndices(iMicro), iMicroNeutSelfEnergy
            IF (MicroIndices(iMicro) == iMicroNeutSelfEnergy) THEN
              iNeutSelfEnergy_in_Table = iMicro
            ELSE IF (MicroIndices(iMicro) == iMicroProtSelfEnergy) THEN
              iProtSelfEnergy_in_Table = iMicro
            END IF
        END DO		

        IF ( iNeutSelfEnergy_in_Table == 0 ) THEN
          WRITE(*,*) 'No Vector self energy for Neutrons in Table, will be set to zero!!!'
        END IF
        IF ( iProtSelfEnergy_in_Table == 0 ) THEN
          WRITE(*,*) 'No Proton self energy for Neutrons in Table, will be set to zero!!!'
        END IF      

        IF (AllocateEOS) THEN
          ALLOCATE( EOSCompOSE(nRho,nTemp,nYp,nVariablesCompOSE) )
        ENDIF

        ! Now copy into Compose array
        DO iRho=1,nRho
            DO iT=1,nTemp
                DO iYp=1,nYp
                    EOSCompOSE(iRho,iT,iYp,iPressCompOSE) = ThermoTable(iRho,iT,iYp,iPress_in_Table) * ergmev / cm3fm3
                    EOSCompOSE(iRho,iT,iYp,iEntropyCompOSE) = ThermoTable(iRho,iT,iYp,iEntr_in_Table)
                    ! This is how I used to do it
                    ! TotalInternalEnergy = (1.0_dp + ThermoTable(iRho,iT,iYp,iEps_in_Table)) * NeutronMass
                    ! EOSCompOSE(iRho,iT,iYp,iInternalEnergyDensityCompOSE) = &
                    !     TotalInternalEnergy*nbCompOSE(iRho)/cm3fm3 * ergmev / RhoCompOSE(iRho) - cvel**2.0d0
                    
                    ! This is how it should be done I think
                    !EOSCompOSE(iRho,iT,iYp,iInternalEnergyDensityCompOSE) = &
                    !    ThermoTable(iRho,iT,iYp,iEps_in_Table) * ergmev / rmu
                    EOSCompOSE(iRho,iT,iYp,iInternalEnergyDensityCompOSE) = &
                        (ThermoTable(iRho,iT,iYp,iEps_in_Table) + 1.0_dp) * &
                         NeutronMass * ergmev / rmu - cvel**2.0_dp
          
                    ! now handle the chemical potentials. These include rest mass already.
                    EOSCompOSE(iRho,iT,iYp,iNeutronChemPotCompOSE) = ThermoTable(iRho,iT,iYp,iMub_in_Table) + NeutronMass
                    EOSCompOSE(iRho,iT,iYp,iProtonChemPotCompOSE) = &
                        ThermoTable(iRho,iT,iYp,iMuq_in_Table) + ThermoTable(iRho,iT,iYp,iMub_in_Table) + ProtonMass
                    
                    ! Take care of microscopic quantities
                    IF (iNeutSelfEnergy_in_Table > 0 ) THEN
                        EOSCompOSE(iRho,iT,iYp,iNeutronSelfEnergyCompOSE) = &
                            MicroTable(iRho,iT,iYp,iNeutSelfEnergy_in_Table)
                    ELSE
                        EOSCompOSE(iRho,iT,iYp,iNeutronSelfEnergyCompOSE) = 0.0_dp
                    ENDIF
                    IF (iProtSelfEnergy_in_Table > 0 ) THEN
                      EOSCompOSE(iRho,iT,iYp,iProtonSelfEnergyCompOSE) = &
                          MicroTable(iRho,iT,iYp,iProtSelfEnergy_in_Table)
                    ELSE
                        EOSCompOSE(iRho,iT,iYp,iProtonSelfEnergyCompOSE) = 0.0_dp
                    ENDIF

                    ! Finally set other things to zero for now, CompOSE does not tell you what they are
                    EOSCompOSE(iRho,iT,iYp,iHeavyBindingEnergyCompOSE) = 0.0_dp
                    EOSCompOSE(iRho,iT,iYp,iThermalEnergyCompOSE)      = 0.0_dp
                END DO
            END DO
        END DO
        
        ! Now worry about the composition
        IF (nPairs > 3) THEN
            fix_composition = .true.
            WRITE(*,*) 'Handling more than 3 Species, be careful what you are writing in the composition!'
        ELSE
            fix_composition = .false.
        END IF
        
        !fix_composition = .false.
        WRITE(*,*)
        WRITE(*,*) 'Composition table has these particles: '
        DO iCompo=1,nPairs
            IF (CompoIndices(iCompo) == iCompoProton) THEN
                iXp_in_Table = iCompo
            ELSE IF (CompoIndices(iCompo) == iCompoNeutron) THEN
                iXn_in_Table = iCompo
            ELSE IF (CompoIndices(iCompo) == iCompoHe4) THEN
                iXa_in_Table = iCompo
            END IF
            WRITE(*,*) 'Particle ', CompoIndices(iCompo), 'present'
        END DO
                
        DO iRho=1,nRho
            DO iT=1,nTemp
                DO iYp=1,nYp
                    EOSCompOSE(iRho,iT,iYp,iProtonMassFractionCompOSE)  = &
                        CompoTable(iRho,iT,iYp,iXp_in_Table)
                    EOSCompOSE(iRho,iT,iYp,iNeutronMassFractionCompOSE) = &
                        CompoTable(iRho,iT,iYp,iXn_in_Table)
                    EOSCompOSE(iRho,iT,iYp,iAlphaMassFractionCompOSE)   = &
                        CompoTable(iRho,iT,iYp,iXa_in_Table) * 4.0d0
                
                    IF (fix_composition) THEN
                        NewAbarNumerator = YhTable(iRho,iT,iYp,1) * AbarTable(iRho,iT,iYp,1)
                        NewFracHeavy = YhTable(iRho,iT,iYp,1)
                        NewZbarNumerator = ZbarTable(iRho,iT,iYp,1) * YhTable(iRho,iT,iYp,1)
                        NewHeavyDenominator = AbarTable(iRho,iT,iYp,1) * YhTable(iRho,iT,iYp,1)
                        DO iCompo=4,nPairs
                            IF (iCompo == iCompoHe3) THEN
                                Light_Abar = 3
                                Light_Zbar = 2
                            ELSEIF (iCompo == iCompoH3) THEN
                                Light_Abar = 3
                                Light_Zbar = 1
                            ELSEIF (iCompo == iCompoH2) THEN
                                Light_Abar = 2
                                Light_Zbar = 1
                            ENDIF
                            
                            IF (CompoTable(iRho,iT,iYp,iCompo) > 0) THEN
                                NewAbarNumerator = NewAbarNumerator + CompoTable(iRho,iT,iYp,iCompo) * Light_Abar
                                NewZbarNumerator = NewZbarNumerator + CompoTable(iRho,iT,iYp,iCompo) * Light_Zbar
                                NewHeavyDenominator = NewHeavyDenominator + CompoTable(iRho,iT,iYp,iCompo)
                                NewFracHeavy = NewFracHeavy + CompoTable(iRho,iT,iYp,iCompo)
                            ENDIF
                        ENDDO
                        
                        IF (NewHeavyDenominator .eq. 0.0d0) THEN
                            EOSCompOSE(iRho,iT,iYp,iHeavyMassNumberCompOSE) = 0.0_dp
                            EOSCompOSE(iRho,iT,iYp,iHeavyChargeNumberCompOSE) = 0.0_dp
                            EOSCompOSE(iRho,iT,iYp,iHeavyMassFractionCompOSE) = 0.0_dp
                        ELSE
                            EOSCompOSE(iRho,iT,iYp,iHeavyMassNumberCompOSE) = &
                                NewAbarNumerator / NewHeavyDenominator
                            EOSCompOSE(iRho,iT,iYp,iHeavyChargeNumberCompOSE) = &
                                NewZbarNumerator / NewHeavyDenominator
                            EOSCompOSE(iRho,iT,iYp,iHeavyMassFractionCompOSE) =  &
                                NewFracHeavy * EOSCompOSE(iRho,iT,iYp,iHeavyMassNumberCompOSE)
                        ENDIF
                    ELSE
                        EOSCompOSE(iRho,iT,iYp,iHeavyMassFractionCompOSE) = YhTable(iRho,iT,iYp,1) * AbarTable(iRho,iT,iYp,1)
                        EOSCompOSE(iRho,iT,iYp,iHeavyMassNumberCompOSE) = AbarTable(iRho,iT,iYp,1)
                        EOSCompOSE(iRho,iT,iYp,iHeavyChargeNumberCompOSE) = ZbarTable(iRho,iT,iYp,1)
                    ENDIF
                    
                END DO
            END DO
        END DO
        
        ! Fix the composition???? How to handle light clusters? No idea.
        ! This is a VERY PRELIMINARY FIX
        DO iRho=1,nRho
            DO iT=1,nTemp
                DO iYp=1,nYp
                    EOSCompOSE(iRho,iT,iYp,iProtonMassFractionCompOSE) = CompoTable(iRho,iT,iYp,iXp_in_Table)
                    EOSCompOSE(iRho,iT,iYp,iNeutronMassFractionCompOSE) = CompoTable(iRho,iT,iYp,iXn_in_Table)
                    EOSCompOSE(iRho,iT,iYp,iAlphaMassFractionCompOSE) = CompoTable(iRho,iT,iYp,iXa_in_Table)
                    EOSCompOSE(iRho,iT,iYp,iHeavyMassFractionCompOSE) = YhTable(iRho,iT,iYp,1)
                    EOSCompOSE(iRho,iT,iYp,iHeavyMassNumberCompOSE) = AbarTable(iRho,iT,iYp,1)
                    EOSCompOSE(iRho,iT,iYp,iHeavyChargeNumberCompOSE) = ZbarTable(iRho,iT,iYp,1)
                END DO
            END DO
        END DO		
        
        DEALLOCATE( nbCompOSE )
        DEALLOCATE( ThermoTable )
        DEALLOCATE( CompoTable )
        DEALLOCATE( ThermoIndices )
        DEALLOCATE( CompoIndices )
        DEALLOCATE( ZbarTable )
        DEALLOCATE( AbarTable )
        DEALLOCATE( YhTable )
        
    END SUBROUTINE ReadCompOSEHDFTable
    
    
    !===============================================!
    !         ATTRIBUTE INPUT ROUTINES
        !===============================================!
        ! These are taken directly from compose
    !=======================================================================
    !> Read an integer4 attribute in a hdf5 file
    subroutine hdf5_read_int4_attr0D(id, name, attr)
        
        implicit none
        
        integer(kind=hid_t), intent(in) :: id              !< id of the file/group where the attribute will be read
        character(len=*), intent(in) :: name              !< name of the attribute
        integer(kind=4), intent(inout) :: attr             !< attribute value
        
        integer(kind=hsize_t), dimension(1) :: dim1D       !< dimension of the attribute: here = 1
        integer(kind=hid_t) :: attr_id                     !< attribute id
        integer(kind=hid_t) :: atype_id                    !< attribute type id
        
        dim1D(1) = 1
        
        call h5aopen_name_f(id, name, attr_id, h5err)
        if(h5err/=0) then
            print *,'Attribute ',name,'not found in the hdf5 file.'
            return
        end if                                  
        call h5aget_type_f(attr_id, atype_id, h5err)
        call h5aread_f(attr_id, atype_id, attr, dim1D, h5err)
        call h5aclose_f(attr_id, h5err)
        
    end subroutine hdf5_read_int4_attr0D                                                                                                                

END MODULE wlCompOSEInterface
        
