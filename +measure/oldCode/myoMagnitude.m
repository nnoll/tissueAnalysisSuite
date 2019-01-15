function [ mM, mov ] = myoMagnitude( time_series )
    % MYOM AGNITUDE 
    
    for t = 1:length(time_series)
       if (~isempty(time_series(t).myosin_2D))
          mM(:,:,t) = sqrt( (time_series(t).myosin_2D{1,1}+.5*time_series(t).myosin_trace).^2 + ...
                      2*(time_series(t).myosin_2D{1,2}).^2 + ...
                      (time_series(t).myosin_2D{2,2}+.5*time_series(t).myosin_trace).^2);
       end
    end
    
    for t = 1:size(mM,3)
       pcolor(mM(:,:,t));
       shading flat
       mov(t) = getframe();
    end


end

