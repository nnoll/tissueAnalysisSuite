function [ IDs, labels ] = returnIDMap( this, t )
    %RETURN ID MAP 

    IDs = find((this.t0 <= t) & (t <= (this.t0 + cellfun(@length,this.label) - 1)));
    labels = zeros(size(IDs));
    for ii = 1:length(IDs)
        jj = t - this.t0(IDs(ii)) + 1;
        labels(ii) = this.label{IDs(ii)}(jj);
    end
    
end

