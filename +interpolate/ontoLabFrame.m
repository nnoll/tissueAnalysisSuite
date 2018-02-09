function [ X, Y, Phi ] = ontoLabFrame( x, y, phi, goodTiles, ERes )
    % ONTO LAB FRAME 

    mask = zeros(size(phi,1),size(phi,2));
    mask(goodTiles) = 1;
    mask(ERes > 10) = 0;
    mask = mask(3:(size(mask,1)-2),3:(size(mask,2)-2));

    goodIdx = mask==1;
    
    x = x(3:end-2);
    y = y(3:end-2);
    
    [yGI,xGI] = find(goodIdx);

    rGI = 1.0*[x(xGI);y(yGI)]';
    D = pdist2(rGI,rGI);
    L = exp(-D.^2 / (2*125^2) );
    L = bsxfun(@rdivide,L,sum(L,2));
    
    [xg,yg] = meshgrid(x,y);
    [X,Y] = meshgrid(linspace(1,1738,218),linspace(1,2050,257));
%     [X,Y] = meshgrid(linspace(1,1024,218),linspace(1,1024,257));
    
    for ii = 1:size(phi,3)
        z = phi(3:(size(phi,1)-2),3:(size(phi,2)-2),ii);
        dat = L*z(goodIdx);
        
        Zg = griddata(xg(goodIdx),yg(goodIdx),dat,X,Y,'cubic'); 
        Phi(:,:,ii) = Zg;
    end

end

