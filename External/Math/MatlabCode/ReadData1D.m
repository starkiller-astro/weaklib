function [ Data1D ] = ReadData1D( FileName, N )

    fid = fopen( FileName, 'rb' );
    Data1D = fread(fid, N, 'double' );
    fclose( fid );

end

