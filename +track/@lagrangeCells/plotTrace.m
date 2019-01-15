function [ trace ] = plotTrace( this, id, L )
    % PLOT TRACE 
    
    trace = zeros(size(L,1),size(L,2),3,length(this.label{id}));
    for t = this.t0(id):(length(this.label{id})+this.t0(id)-1)
        tmp = zeros(size(L,1),size(L,2));
        ii = t - this.t0(id) + 1;
        
        trace(:,:,1,t) = L(:,:,t) == 0;
        trace(:,:,2,t) = L(:,:,t) == 0;
        trace(:,:,3,t) = L(:,:,t) == 0;
        tmp(L(:,:,t) == this.label{id}(ii)) = 1;
        trace(:,:,1,t) = trace(:,:,1,t) + tmp;
    end


end

