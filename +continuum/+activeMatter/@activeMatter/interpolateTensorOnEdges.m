function [ T ] = interpolateTensorOnEdges( this, Tv )
%INTERPOLATETENSORONEDGES

    Te = zeros(3,3,size(this.d0,1));
    for ii = 1:3
        for jj = 1:3
            Te(ii,jj,:) = .5*abs(this.d0)*squeeze(Tv(ii,jj,:));
        end
    end
    
    rB = this.d0 * this.Mesh.X;
    rB = bsxfun(@rdivide,rB,this.Lp);
    rB = reshape(rB',[3,1,size(rB,1)]);
    
    T = squeeze(mtimesx(permute(rB,[2,1,3]),mtimesx(Te,rB)));

end

