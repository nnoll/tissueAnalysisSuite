function [ vIGrad ] = reinterpolateTensor( xG, yG, vGrad, x, y )
    % REINTERPOLATE TENSOR 
    
    [xG,yG] = meshgrid(xG,yG);
    [x,y] = meshgrid(x,y);
    vIGrad = cell(size(vGrad));
    for t = 1:length(vGrad)
        vIGrad{t} = zeros(2,2,size(x,1),size(x,2));
        for ii = 1:2
            for jj = 1:2
                vIGrad{t}(ii,jj,:,:) = interp2(xG,yG,squeeze(vGrad{t}(ii,jj,:,:)),x,y,'spline');
            end
        end
    end

end

