function [ Data2D ] = ReadData2D( FileName, M, N )

    fid = fopen( FileName, 'rb' );
    Data2D = fread( fid, M * N, 'double' );
    fclose( fid );
    
    Data2D = reshape( Data2D, [ M, N ] );

end