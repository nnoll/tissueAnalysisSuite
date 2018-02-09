function [ Xg, Yg, SigmaI ] = interpolateTensor( x, y, sigma, Eres )
    % INTERPOLATE HM 
    
    [~, ~, SigmaI(:,:,1)] = measure.interpolateHM(x,y,sigma(:,:,1),Eres,0);
    [~, ~, SigmaI(:,:,2)] = measure.interpolateHM(x,y,sigma(:,:,2),Eres,1);
    [Xg, Yg, SigmaI(:,:,3)] = measure.interpolateHM(x,y,sigma(:,:,3),Eres,0);
    
end

