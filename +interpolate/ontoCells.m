function [ PhiC ] = ontoCells( x, y, goodTiles, PhiG, Struct )
    % ONTO CELLS 
    
    mask = zeros(size(x));
    mask(goodTiles) = 1;
    [x,y] = meshgrid(x,y);
    
    x = x(mask == 1);
    y = y(mask == 1);
    Phi = zeros(length(x),size(PhiG,3));
    for ii = 1:size(Phi,2)
       tmp = PhiG(:,:,ii);
       Phi(:,ii) = tmp(mask==1); 
    end
    
    [ ~, ~, ~, ~, iCells ] = fitDual.ATN.computeDiffOperators( Struct, 1 );  
    Rc = zeros(length(iCells),2);
    for c = 1:length(iCells)
       Rc(c,:) = Struct.Cdat(iCells(c)).centroid.coord; 
    end
    
    PhiC = zeros(size(Rc,1),size(Phi,2));
    for ii = 1:size(Phi,2)
        F = scatteredInterpolant(x,y,Phi(:,ii),'natural');
        PhiC(:,ii) = F(Rc(:,1),Rc(:,2));
    end
end

