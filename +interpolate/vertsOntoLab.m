function [ phiLab ] = vertsOntoLab( smScale, rv, phi )
%VERTSONTOLAB Summary of this function goes here
%   Detailed explanation goes here

    [X,Y] = meshgrid(linspace(1,1738,218),linspace(1,2050,257));

    D = pdist2(rv,rv);
    L = exp(-D.^2 / (2*smScale^2) );
    L = bsxfun(@rdivide,L,sum(L,2));
    F = scatteredInterpolant(rv(:,1),rv(:,2),L*phi,'natural','none');
    phiLab = F(X,Y);

end

