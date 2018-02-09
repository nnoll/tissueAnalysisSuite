function [ ] = storeLength( this, bStruct )
    %STORE LENGTH 


    for b = 1:length(this.label)
       for ii = 1:length(this.label{b})
           t = this.times{b}(ii);
           rB = [bStruct(t).Vdat(bStruct(t).Bdat(this.label{b}(ii)).verts(1)).vertxcoord,bStruct(t).Vdat(bStruct(t).Bdat(this.label{b}(ii)).verts(1)).vertycoord] - ...
                [bStruct(t).Vdat(bStruct(t).Bdat(this.label{b}(ii)).verts(2)).vertxcoord,bStruct(t).Vdat(bStruct(t).Bdat(this.label{b}(ii)).verts(2)).vertycoord] ;
           this.r{b}(ii) = sqrt(sum(rB.^2));
       end
        
    end
end

