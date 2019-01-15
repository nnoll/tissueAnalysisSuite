function [ xI, yI, vI ] = refineGrid( x, y, v )
    % RE INTERPOLATE 

    xI = linspace(min(x),max(x),218);
    yI = linspace(min(y),max(y),257);
    [Xq,Yq] = meshgrid(xI,yI);

    [ x, y ] = meshgrid(x,y);
    vI = cell(size(v));
    for t = 1:length(v)
        vI{t} = zeros(length(yI),length(xI),2);
        vI{t}(:,:,1) = interp2(x,y,squeeze(v{t}(1,:,:)),Xq,Yq,'spline');
        
        vI{t}(:,:,2) = interp2(x,y,squeeze(v{t}(2,:,:)),Xq,Yq,'spline');
    end
    
    
end

