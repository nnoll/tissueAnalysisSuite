function [ c ] = correlateVelocity( this, vP, timePts )
    %CORRELATE VELOCITY 

    if (nargin == 2)
        timePts = 1:length(vP);
    end
    
    for t = timePts
        c(t) = median( dot(this.v{t},vP{t},2) ./ sqrt(dot(this.v{t},this.v{t},2).*dot(vP{t},vP{t},2)) );
    end
end

