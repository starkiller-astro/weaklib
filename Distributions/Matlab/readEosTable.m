function [ D, T, Y, nD, nT, nY, P, S, E, Me, Mp, Mn, Xp, Xn, Xa, Xh,...
           Zh, Ah, Eh, Eth, Gm ] = readEosTable( fileName )

    % Reads HDF5 EOS Table

    disp( fprintf( 'INFO: Reading EOS from file: %s', fileName ) );

    % Independent Variables:
    D = h5read( fileName, '/ThermoState/Density' );
    T = h5read( fileName, '/ThermoState/Temperature' );
    Y = h5read( fileName, '/ThermoState/Electron Fraction' );
    
    Dims = h5read( fileName, '/ThermoState/Dimensions' );
    
    nD = Dims( 1 );
    nT = Dims( 2 );
    nY = Dims( 3 );

    minD = min( D );
    maxD = max( D );
    minT = min( T );
    maxT = max( T );
    minY = min( Y );
    maxY = max( Y );
    
    disp( fprintf( '  INFO: nD, nT, nY = %i, %i, %i', nD, nT, nY ) );
    
    disp( fprintf( '  INFO: Range in D: %d to %d', minD, maxD ) );
    disp( fprintf( '  INFO: Range in T: %d to %d', minT, maxT ) );
    disp( fprintf( '  INFO: Range in Y: %d to %d', minY, maxY ) );

    % Dependent Variables:
    
    OS = h5read( fileName, '/DependentVariables/Offsets' );
    
    % Pressure:
    iP = h5read( fileName, '/DependentVariables/iPressure' );
    P  = h5read( fileName, '/DependentVariables/Pressure' );
    P  = 10.^( P ) - OS(iP);
    
    % Entropy Per Baryon:
    iS = h5read( fileName, '/DependentVariables/iEntropyPerBaryon' );
    S  = h5read( fileName, '/DependentVariables/Entropy Per Baryon' );
    S  = 10.^( S ) - OS(iS);
    
    % Ineternal Energy Density:
    iE = h5read( fileName, '/DependentVariables/iInternalEnergyDensity' );
    E  = h5read( fileName, '/DependentVariables/Internal Energy Density' );
    E  = 10.^( E ) - OS(iE);
    
    % Electron Chemical Potential:
    iMe = h5read( fileName, '/DependentVariables/iElectronChemicalPotential' );
    Me  = h5read( fileName, '/DependentVariables/Electron Chemical Potential' );
    Me  = 10.^( Me ) - OS(iMe);
    
    % Proton Chemical Potential:
    iMp = h5read( fileName, '/DependentVariables/iProtonChemicalPotential' );
    Mp  = h5read( fileName, '/DependentVariables/Proton Chemical Potential' );
    Mp  = 10.^( Mp ) - OS(iMp);
    
    % Neutron Chemical Potential:
    iMn = h5read( fileName, '/DependentVariables/iNeutronChemicalPotential' );
    Mn  = h5read( fileName, '/DependentVariables/Neutron Chemical Potential' );
    Mn  = 10.^( Mn ) - OS(iMn);
    
    % Proton Mass Fraction:
    iXp = h5read( fileName, '/DependentVariables/iProtonMassFraction' );
    Xp  = h5read( fileName, '/DependentVariables/Proton Mass Fraction' );
    Xp  = 10.^( Xp ) - OS(iXp);
    
    % Neutron Mass Fraction:
    iXn = h5read( fileName, '/DependentVariables/iNeutronMassFraction' );
    Xn  = h5read( fileName, '/DependentVariables/Neutron Mass Fraction' );
    Xn  = 10.^( Xn ) - OS(iXn);
    
    % Alpha Mass Fraction:
    iXa = h5read( fileName, '/DependentVariables/iAlphaMassFraction' );
    Xa  = h5read( fileName, '/DependentVariables/Alpha Mass Fraction' );
    Xa  = 10.^( Xa ) - OS(iXa);
    
    % Heavy Mass Fraction:
    iXh = h5read( fileName, '/DependentVariables/iHeavyMassFraction' );
    Xh  = h5read( fileName, '/DependentVariables/Heavy Mass Fraction' );
    Xh  = 10.^( Xh ) - OS(iXh);
    
    % Heavy Charge Number:
    iZh = h5read( fileName, '/DependentVariables/iHeavyChargeNumber' );
    Zh  = h5read( fileName, '/DependentVariables/Heavy Charge Number' );
    Zh  = 10.^( Zh ) - OS(iZh);
    
    % Heavy Mass Number:
    iAh = h5read( fileName, '/DependentVariables/iHeavyMassNumber' );
    Ah  = h5read( fileName, '/DependentVariables/Heavy Mass Number' );
    Ah  = 10.^( Ah ) - OS(iAh);
    
    % Heavy Binding Energy:
    iEh = h5read( fileName, '/DependentVariables/iHeavyBindingEnergy' );
    Eh  = h5read( fileName, '/DependentVariables/Heavy Binding Energy' );
    Eh  = 10.^( Eh ) - OS(iEh);
    
    % Thermal Energy:
    iEth = h5read( fileName, '/DependentVariables/iThermalEnergy' );
    Eth  = h5read( fileName, '/DependentVariables/Thermal Energy' );
    Eth  = 10.^( Eth ) - OS(iEth);
    
    % Gamma1:
    iGm = h5read( fileName, '/DependentVariables/iGamma1' );
    Gm  = h5read( fileName, '/DependentVariables/Gamma1' );
    Gm  = 10.^( Gm ) - OS(iGm);
    
end

