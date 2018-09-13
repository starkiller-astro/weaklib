function [r, rho, T, Ye, E, opacity] = ...
ReadProfileOutputHDF( fileName, opacityname )

%%% option of opacity name is
%%%    [ 'EmAb_Electron',  'EmAb_ElecAnti',
%%%      'ES_Electron',    'ES_ElecAnti']

r = h5read( fileName, '/ProfileInfo/Radius' );
rho = h5read( fileName, '/ProfileInfo/Density' );
T = h5read( fileName, '/ProfileInfo/Temperature' );
Ye = h5read( fileName, '/ProfileInfo/Electron Fraction' );
E = h5read( fileName, '/ProfileInfo/Energy' );

opacity =  h5read( fileName, ['/Opacities/',opacityname] );

end