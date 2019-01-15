function [ sigma ] = computeStressTensor( this )
    % COMPUTE STRESS TENSOR 
 
    [ F ] = this.computeForces();
    [ r ] = this.computePrimalVerts();
    sigma = zeros(size(this.d0,2),3);
    
    for v = 1:size(this.tri,1)
       for ii = 1:3
           % xx component
           sigma(this.tri(v,ii),1) = sigma(this.tri(v,ii),1) + F(v,ii,1)*r(v,1);
           % xy component
           sigma(this.tri(v,ii),2) = sigma(this.tri(v,ii),2) + .5*(F(v,ii,1)*r(v,2) + F(v,ii,2)*r(v,1));
           % yy component
           sigma(this.tri(v,ii),3) = sigma(this.tri(v,ii),3) + F(v,ii,2)*r(v,2);
       end
    end
    
    [ A1 ] = this.computeEdgeAreas();
    A0 = (abs(this.d0'==1) * A1(:,1)) + (abs(this.d0'==-1) * A1(:,2));
%     A0
%     sigma = bsxfun(@rdivide,sigma,A0);
    
    % Set sigma of boundary to zero.
    bndryEdges = sum(abs(this.d1),1) == 1;
    bndryCells = sum(abs(this.d0(bndryEdges,:)),1) >=1;
    sigma(bndryCells,:) = 0;
    
end

