function [ sM, mov ] = tensorMagnitude( sigma )
%TENSORMAGNITUDE 

    for t = 1:length(sigma)
       if (~isempty(sigma{t}))
          sM(:,:,t) = sqrt( sigma{t}(:,:,1).^2 + 2*sigma{t}(:,:,2).^2 + sigma{t}(:,:,3).^2 );
       end
    end
    
    clim = [inf,-inf];
    for t = 1:size(sM,3)
       pcolor(sM(:,:,t));
       shading flat
       mov(t) = getframe();
    end

end

