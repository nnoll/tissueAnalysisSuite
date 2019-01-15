function [ ] = storeChem( this, bStruct, ind )
%STORECHEM 

    if (nargin == 2)
        ind = 1;
    end
    
    for b = 1:length(this.label)
       for ii = 1:length(this.label{b})
           t = this.times{b}(ii);
           this.myo{b}(ii) = bStruct(t).Bdat(this.label{b}(ii)).chem(ind);
       end
        
    end

end

