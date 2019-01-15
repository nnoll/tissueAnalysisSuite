function [ sigmaI ] = globalStress( Struct, L, timePts )
%GLOBALSTRESS Summary of this function goes here
%   Detailed explanation goes here

    for t = timePts
        [ sigma, i0 ] = measure.stressTensor( Struct(t), L(:,:,t) );
        [ Xg, Yg, sigmaI{t} ] = measure.interpolateCellHM( sigma, i0, L(:,:,t) );
    end
    
end

