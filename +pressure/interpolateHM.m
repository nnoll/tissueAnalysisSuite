function [ Xg, Yg, Zg ] = interpolateHM( x, y, z )
    % INTERPOLATE HM 
    
    g = fspecial('gaussian',2,2);
    z = imfilter(z,g);
    
    [yg,xg] = ndgrid(y,x);
    [Yg,Xg] = ndgrid(linspace(min(y),max(y),100),linspace(min(x),max(x),100));
    ZF = griddedInterpolant(yg,xg,z,'linear');
    Zg = ZF(Yg,Xg);

end

