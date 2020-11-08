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
         ( OpacityTable, nOpac_EmAb = 4, nOpac_Iso = 1, nMom_Iso = 1, &
           nOpac_NES = 0, nMom_NES = 0, nOpac_Pair = 0, nMom_Pair = 0, &
           nOpac_Brem = 0, nMom_Brem = 0, nPointsE = 10, nPointsEta = 0 )

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

  ! -- Opacity -- 

  ASSOCIATE( EmAb => OpacityTable % EmAb )

  EmAb % Names &
    = [ 'Electron Neutrino               ', &
        'Electron Antineutrino           ', &
        'Mu/Tau Neutrino                 ', &
        'Mu/Tau Antineutrino             ' ]

  EmAb % Units &
    = [ 'Per Centimeter                  ', &
        'Per Centimeter                  ', &
        'Per Centimeter                  ', &
        'Per Centimeter                  ' ]

  ASSOCIATE( Chi_Nu_e => EmAb % Opacity(iNu_e) % Values )

  DO iY = 1, EmAb % nPoints(4)
    DO iT = 1, EmAb % nPoints(3)
      DO iD = 1, EmAb % nPoints(2)
        DO iE = 1, EmAb % nPoints(1)

          Chi_Nu_e(iE,iD,iT,iY) = 1.0_dp

        END DO
      END DO
    END DO
  END DO

  END ASSOCIATE ! Chi_Nu_e

  END ASSOCIATE ! EmAb

  ! -- Inelastic Scattering --

  ASSOCIATE( Scat_Iso => OpacityTable % Scat_Iso )

  Scat_Iso % Names &
    = [ 'Elastic Scattering on Nuclei  ' ]

!  Scat_Iso % Species &
!    = [ 'Electron Neutrino             ' ]

  Scat_Iso % Units &
    = [ 'Per Centimeter                ' ]


  ASSOCIATE( Sig_Nu_e => Scat_Iso % Kernel(iNu_e) % Values )

  DO il = 1, Scat_Iso % nMoments
    DO iY = 1, Scat_Iso % nPoints(4)
      DO iT = 1, Scat_Iso % nPoints(3)
        DO iD = 1, Scat_Iso % nPoints(2)
          DO iE = 1, Scat_Iso % nPoints(1)

            Sig_Nu_e(iE,iD,iT,iY,il) = 2.0_dp
    
          END DO
        END DO
      END DO
    END DO
  END DO

  END ASSOCIATE ! Sig_Nu_e

  END ASSOCIATE ! Scat_Iso


  CALL DescribeOpacityTable( OpacityTable )

  CALL DeAllocateOpacityTable( OpacityTable )

END PROGRAM wlOpacityFieldsTest
