function [ ] = plot3DScalar( field, time_series )
%PLOT3DSCALAR 
    
    X = double(time_series(1).eGrids{1});
    Y = double(time_series(1).eGrids{2});
    Z = double(time_series(1).eGrids{3});
    
    % Refine grid
    [x,y] = meshgrid(1:size(X,2),1:size(X,1));
    [xr,yr] = meshgrid(1:.25:size(X,2),1:.25:size(X,1));
    
    X = interp2(x,y,X,xr,yr);
    Y = interp2(x,y,Y,xr,yr);
    Z = interp2(x,y,Z,xr,yr);

%     field = field(:,size(field,2):-1:1,:);
    
    if (size(field,1) ~= size(Z,1) && size(field,2) ~= size(Z,2))
        [Xi,Yi] = meshgrid(1:size(field,2),1:size(field,1));
        [Xg,Yg] = meshgrid(linspace(1,size(field,2),size(Z,2)),linspace(1,size(field,1),size(Z,1)));
        
        if (length(size(field))==2)
            phi = interp2(Xi,Yi,field,Xg,Yg);
        else
            for ii = 1:3
                phi(:,:,ii) = interp2(Xi,Yi,field(:,:,ii),Xg,Yg);
            end
        end
    else
        phi = field;
    end
    
    surf(Z,Y,-X,phi,'EdgeColor','none')
    shading flat
    axis equal
    view([0,0])
    
end

