function [ Xg, Yg, Zg ] = interpolateHM( x, y, z, goodTiles, ERes )
    % INTERPOLATE HM 
    
    mask = zeros(size(z));
    mask(goodTiles) = 1;
    mask(ERes > 10) = 0;
    
    z = z(3:(size(z,1)-2),3:(size(z,2)-2));
    mask = mask(3:(size(mask,1)-2),3:(size(mask,2)-2));
    x = x(3:end-2);
    y = y(3:end-2);
    
    goodIdx = mask==1;
    [xg,yg] = meshgrid(x,y);
    [Yg,Xg] = meshgrid(linspace(min(y),max(y),218),linspace(min(x),max(x),257));
    Zg = griddata(yg(goodIdx),xg(goodIdx),z(goodIdx),Yg,Xg,'cubic'); 
    Zg = medfilt2(Zg,[30,30],'symmetric');
    
%     % Symmetrize across L/R axis.
%     if (nargin == 4 || mode == 0)
%         Zg(1:128,:) = .5*(Zg(1:128,:) + Zg(257:-1:130,:));     
%         Zg(257:-1:130,:) = Zg(1:128,:);
%     else
%         Zg(1:128,:) = -.5*(Zg(1:128,:) - Zg(257:-1:130,:)); 
%         Zg(257:-1:130,:) = -Zg(1:128,:);
%     end
    

    
end

