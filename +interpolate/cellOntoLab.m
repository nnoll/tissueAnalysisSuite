function [ sigmaLab ] = cellOntoLab( Struct, smScale, ERes, xG, yG )
    %CELLONTOLAB Summary of this function goes here
    %   Detailed explanation goes here

    [X,Y] = meshgrid(linspace(1,1738,218),linspace(1,2050,257));
%     [X,Y] = meshgrid(linspace(1,1024,218),linspace(1,1024,257));
    
    for t = 1:length(Struct)
        sigmaLab{t} = zeros(size(X));
        
        [ ~, ~, ~, ~, iCells ] = fitDual.ATN.computeDiffOperators( Struct(t), 1 );
        if (nargin  > 2)
            goodCells = measure.goodCells(Struct(t),ERes,xG,yG);
            iCells = iCells(ismember(iCells,goodCells));
        end
        
        sigma = zeros(length(iCells),3);
        Rc = zeros(length(iCells),2);

        for c = 1:length(iCells)
            sigma(c,1) = Struct(t).Cdat(iCells(c)).stress(1,1);
            sigma(c,2) = Struct(t).Cdat(iCells(c)).stress(1,2);
            sigma(c,3) = Struct(t).Cdat(iCells(c)).stress(2,2);

            Rc(c,:) = Struct(t).Cdat(iCells(c)).centroid.coord;
        end
        
        D = pdist2(Rc,Rc);
        L = exp(-D.^2 / (2*smScale^2) );
        L = bsxfun(@rdivide,L,sum(L,2));
    
        for ii = 1:3
            F = scatteredInterpolant(Rc(:,1),Rc(:,2),L*sigma(:,ii),'natural','none');
            sigmaLab{t}(:,:,ii) = F(X,Y);
        end
        
        [ sigmaLab{t} ] = generate.symmetrizedStress( sigmaLab{t} );
    end
end

