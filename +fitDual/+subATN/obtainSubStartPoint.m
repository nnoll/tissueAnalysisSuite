function [ q, theta, tri, i0 ] = obtainSubStartPoint( Struct, XRange, YRange, r0 )
    % FIT PRESSURE 

    if (nargin == 4)
        [ Q, ~, ~, ~, ~, i0 ] = fitDual.subATN.fitSubTensionGraph( Struct, XRange, YRange, r0  );
    else
        [ Q, ~, ~, ~, ~, i0 ] = fitDual.subATN.fitSubTensionGraph( Struct, XRange, YRange );
    end
    [ q, i0 ] = fitDual.subATN.rescaleSubQ( Q, i0, Struct, XRange, YRange  );
    [ q, theta, tri ] = fitDual.subATN.fitSubThetaModes( q, Struct, XRange, YRange );
    
end

