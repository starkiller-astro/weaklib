PROGRAM wlOpacityFieldsTest

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

  TYPE(OpacityTableType) :: OpacityTable

  CALL AllocateOpacityTable &
         ( OpacityTable, nOpacA = 4, nOpacB = 1, nMomB = 1, &
           nOpacC = 1, nMomC = 1, nPointsE = 10 )

  ASSOCIATE( EnergyGrid => OpacityTable % EnergyGrid )

  EnergyGrid % Name &
    = 'Comoving Frame Neutrino Energy  '

  EnergyGrid % Unit &
    = 'MeV                             '

  EnergyGrid % MinValue = 3.0d-1
  EnergyGrid % MaxValue = 3.0d+2

  CALL MakeLogGrid &
         ( EnergyGrid % MinValue, EnergyGrid % MaxValue, &
           EnergyGrid % nPoints, EnergyGrid % Values )

  END ASSOCIATE ! EnergyGrid

  ASSOCIATE( ecap => OpacityTable % ecap )

  ecap % Names &
    = [ 'Electron Capture Absoptivity    ', &
        'Electron Capture Absoptivity    ', &
        'Electron Capture Absoptivity    ', &
        'Electron Capture Absoptivity    ' ]

  ecap % Species &
    = [ 'Electron Neutrino               ', &
        'Electron Antineutrino           ', &
        'Mu/Tau Neutrino                 ', &
        'Mu/Tau Antineutrino             ' ]

  ecap % Units &
    = [ 'Per Centimeter                  ', &
        'Per Centimeter                  ', &
        'Per Centimeter                  ', &
        'Per Centimeter                  ' ]

  END ASSOCIATE ! ecap


  ASSOCIATE( scatt_Iso => OpacityTable % scatt_Iso )

    scatt_Iso % Names &
      = [ 'Elastic Scattering on Nuclei  ' ]

    scatt_Iso % Species &
      = [ 'Electron Neutrino             ' ]

    scatt_Iso % Units &
      = [ 'Per Centimeter                ' ]

  END ASSOCIATE ! scatt_Iso


  ASSOCIATE( scatt_nIso => OpacityTable % scatt_nIso )

    scatt_nIso % Names &
      = [ 'Inelastic Scattering on Elect.' ]

    scatt_nIso % Species &
      = [ 'Electron Neutrino             ' ]

    scatt_nIso % Units &
      = [ 'Per Centimeter Per MeV**3     ' ]

  END ASSOCIATE ! scatt_nIso

  CALL DescribeOpacityTable( OpacityTable )

  CALL DeAllocateOpacityTable( OpacityTable )

END PROGRAM wlOpacityFieldsTest
