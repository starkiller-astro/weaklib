PROGRAM wlOpacityFieldsTest

  USE wlKindModule, ONLY: &
    dp
  USE wlOpacityTableModule, ONLY: &
    OpacityTableType, &
    AllocateOpacityTable, &
    DeAllocateOpacityTable, &
    DescribeOpacityTable
  USE wlOpacityFieldsModule, ONLY: &
    iNu_e, iNu_e_bar, iNu_x, iNu_x_bar
  USE wlGridModule, ONLY: &
    MakeLogGrid

  IMPLICIT NONE

  INTEGER :: iE, iEp, iD, iT, iY, il
  TYPE(OpacityTableType) :: OpacityTable

  CALL AllocateOpacityTable &
         ( OpacityTable, nOpacA = 4, nOpacB = 1, nMomB = 1, &
           nOpacB_NES = 0, nMomB_NES = 0, nOpacB_TP = 0, nMomB_TP = 0, &
           nOpacC = 1, nMomC = 1, &
           nPointsE = 10, nPointsEta = 0 )

  ! -- Energy Grid -- 

  ASSOCIATE( EnergyGrid => OpacityTable % EnergyGrid )

  EnergyGrid % Name &
    = 'Comoving Frame Neutrino Energy  '

  EnergyGrid % Unit &
    = 'MeV                             '

  EnergyGrid % MinValue = 3.0d-1
  EnergyGrid % MaxValue = 3.0d+2

  EnergyGrid % LogInterp = 1

  CALL MakeLogGrid &
         ( EnergyGrid % MinValue, EnergyGrid % MaxValue, &
           EnergyGrid % nPoints, EnergyGrid % Values )

  END ASSOCIATE ! EnergyGrid

  ! -- Absorptivity -- 

  ASSOCIATE( thermEmAb => OpacityTable % thermEmAb )

  thermEmAb % Names &
    = [ 'Electron Capture Absoptivity    ', &
        'Electron Capture Absoptivity    ', &
        'Electron Capture Absoptivity    ', &
        'Electron Capture Absoptivity    ' ]

  thermEmAb % Species &
    = [ 'Electron Neutrino               ', &
        'Electron Antineutrino           ', &
        'Mu/Tau Neutrino                 ', &
        'Mu/Tau Antineutrino             ' ]

  thermEmAb % Units &
    = [ 'Per Centimeter                  ', &
        'Per Centimeter                  ', &
        'Per Centimeter                  ', &
        'Per Centimeter                  ' ]

  ASSOCIATE( Chi_Nu_e => thermEmAb % Absorptivity(iNu_e) % Values )

  DO iY = 1, thermEmAb % nPoints(4)
    DO iT = 1, thermEmAb % nPoints(3)
      DO iD = 1, thermEmAb % nPoints(2)
        DO iE = 1, thermEmAb % nPoints(1)

          Chi_Nu_e(iE,iD,iT,iY) = 1.0_dp

        END DO
      END DO
    END DO
  END DO

  END ASSOCIATE ! Chi_Nu_e

  END ASSOCIATE ! thermEmAb

  ! -- Inelastic Scattering --

  ASSOCIATE( scatt_Iso => OpacityTable % scatt_Iso )

  scatt_Iso % Names &
    = [ 'Elastic Scattering on Nuclei  ' ]

  scatt_Iso % Species &
    = [ 'Electron Neutrino             ' ]

  scatt_Iso % Units &
    = [ 'Per Centimeter                ' ]


  ASSOCIATE( Sig_Nu_e => scatt_Iso % Kernel(iNu_e) % Values )

  DO il = 1, scatt_Iso % nMoments
    DO iY = 1, scatt_Iso % nPoints(4)
      DO iT = 1, scatt_Iso % nPoints(3)
        DO iD = 1, scatt_Iso % nPoints(2)
          DO iE = 1, scatt_Iso % nPoints(1)

            Sig_Nu_e(iE,iD,iT,iY,il) = 2.0_dp
    
          END DO
        END DO
      END DO
    END DO
  END DO

  END ASSOCIATE ! Sig_Nu_e

  END ASSOCIATE ! scatt_Iso

  ! -- Elastic Scattering --

  ASSOCIATE( scatt_nIso => OpacityTable % scatt_nIso )

  scatt_nIso % Names &
    = [ 'Inelastic Scattering on Elect.' ]

  scatt_nIso % Species &
    = [ 'Electron Neutrino             ' ]

  scatt_nIso % Units &
    = [ 'Per Centimeter Per MeV**3     ' ]

  ASSOCIATE( Sig_Nu_e => scatt_nIso % Kernel(iNu_e) % Values )

  DO il = 1, scatt_nIso % nMoments
    DO iY = 1, scatt_nIso % nPoints(4)
      DO iT = 1, scatt_nIso % nPoints(3)
        DO iD = 1, scatt_nIso % nPoints(2)
          DO iEp = 1, scatt_nIso % nPoints(1)
            DO iE = 1, scatt_nIso % nPoints(1)

              Sig_Nu_e(iE,iEp,iD,iT,iY,il) = 3.0_dp

            END DO
          END DO
        END DO
      END DO
    END DO
  END DO

  END ASSOCIATE ! Sig_Nu_e

  END ASSOCIATE ! scatt_nIso

  CALL DescribeOpacityTable( OpacityTable )

  CALL DeAllocateOpacityTable( OpacityTable )

END PROGRAM wlOpacityFieldsTest
