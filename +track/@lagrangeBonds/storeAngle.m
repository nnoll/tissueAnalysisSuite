function [ ] = storeAngle( this, bStruct )
%STOREANGLE Summary of this function goes here
%   Detailed explanation goes here

    for b = 1:length(this.label)
       for ii = 1:length(this.label{b})
           t = this.times{b}(ii);
           rB = [bStruct(t).Vdat(bStruct(t).Bdat(this.label{b}(ii)).verts(1)).vertxcoord,bStruct(t).Vdat(bStruct(t).Bdat(this.label{b}(ii)).verts(1)).vertycoord] - ...
                [bStruct(t).Vdat(bStruct(t).Bdat(this.label{b}(ii)).verts(2)).vertxcoord,bStruct(t).Vdat(bStruct(t).Bdat(this.label{b}(ii)).verts(2)).vertycoord] ;
           rB = rB/sqrt(sum(rB.^2));
           theta = atan2(rB(2),rB(1));
           if (theta < 0)
                theta = theta + pi;
           end
           this.theta{b}(ii) = theta;
       end
        
    end
end

