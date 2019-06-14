function [ PN, ERes, r0 ] = returnGridDual( Struct, xG, yG )
%RETURNGRIDDUAL 
%Optional file header info 
%
% Syntax:  [rgb] = trackedVerts(Struct, vertexPairs, cellPairs, skel, T)
%
% Inputs:
%    Struct - Description
%    vertexPairs - Description
%    cellPairs - Description
%    skel - 
%    T - 
%
% Outputs:
%    rgb - Description
%
% Example: 
%    Line 1 of example
%    Line 2 of example
%    Line 3 of example
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

    Zc = 15;
    Z = zeros(length(yG),length(xG));
    
    for ii = 1:(length(yG))
        for jj = 1:(length(xG))
            [Z(ii,jj), nCC(ii,jj)] = measure.numberOfBulkCells(Struct,xG(:,jj),yG(:,ii));
        end
    end

    PN = cell(length(yG),length(xG));
    for ii = 1:(length(yG))
        for jj = 1:(length(xG))
            if ( Z(ii,jj) >= Zc && nCC(ii,jj) == 1 )
                [ tmp, ERes(ii,jj) ] = fitDual.returnSubDual( Struct, 1, xG(:,jj), yG(:,ii) );
                PN{ii,jj} = tmp{1};
                [ ~, ~, ~, ~, ~, r0{ii,jj} ] = fitDual.returnSubGraph( Struct, xG(:,jj), yG(:,ii) );
                [ ~, dV, ~, i0 ] = fitDual.subATN.computeSubDiffOperators( Struct, xG(:,jj), yG(:,ii) );
                D = sqrt( sum( (dV*[Struct.Vdat(i0).vertxcoord;Struct.Vdat(i0).vertycoord]').^2, 2) );
                
%                 [ T, Ta ] = fitDual.compareTension( Struct, tmp{1} );
%                 ERes(ii,jj) / mean(D).^2
%                 corr(T,Ta)
%                 figure(1)
%                 plot.skel(Struct,'k',0)
%                 hold on
%                 tmp{1}.plotPrimal();
%                 figure(2)
%                 scatter(T,Ta)
%                 pause
                
                ERes(ii,jj) = ERes(ii,jj) / (mean(D).^2);
            else
                ERes(ii,jj) = 100;
            end
        end
    end

%     [ scale, goodTiles ] = fitDual.rescaleGridDual(PN,ERes);
   
end

