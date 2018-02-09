function [ q, theta, tri, i0 ] = obtainStartPoint( Struct, extCell )
    % FIT PRESSURE 

    [ Q, ~, ~, ~, ~, i0 ] = fitDual.ATN.fitTensionGraph( Struct, extCell );
    [ q, i0 ] = fitDual.ATN.rescaleQ( Q, i0, Struct, extCell );
    [ q, theta, tri ] = fitDual.ATN.fitThetaModes( q, Struct, extCell );
    
end

