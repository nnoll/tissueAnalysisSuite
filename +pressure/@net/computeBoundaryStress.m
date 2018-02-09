function [ sigmaB, Apatch ] = computeBoundaryStress( this, r )
    % COMPUTE BOUNDARY STRESS 

    if (nargin == 1)
        [ r ] = this.computePrimalVerts();
    end
    r = bsxfun(@minus,r,mean(r,1));
    
    sigmaB = zeros(3,1);
    Rot = [0,1;-1,0];
    
    % Find boundary edges
    bndryE = find(sum(abs(this.d1),1) == 1);
    bndryVs = [];
    for bE = bndryE
        qAB = this.d0(bE,:)*this.q;
        qAB = qAB * Rot';
        TAB = sqrt(sum(qAB.^2));
        bndryV = (this.d1(:,bE)~=0);
        bndryVs = [bndryVs,find(bndryV)];
        S = this.d1(bndryV,bE);
        qAB = -S*qAB;
        TAB = sqrt(sum(qAB.^2));
%         this.plotPrimal()
%         hold on
%         scatter(r(bndryV,1),r(bndryV,2),'g','filled')
%         quiver(r(bndryV,1),r(bndryV,2),qAB(1),qAB(2),'LineWidth',2,'Color','r')
        sigmaB(1) = sigmaB(1) + qAB(1)*r(bndryV,1);
        sigmaB(2) = sigmaB(2) + .5*(qAB(1)*r(bndryV,2)+qAB(2)*r(bndryV,1));
        sigmaB(3) = sigmaB(3) + qAB(2)*r(bndryV,2);
%         sigmaB(1) = sigmaB(1) + qAB(1)*qAB(1)/TAB;
%         sigmaB(2) = sigmaB(2) + qAB(1)*qAB(2)/TAB;
%         sigmaB(3) = sigmaB(3) + qAB(2)*qAB(2)/TAB;
        
    end
    
    [~,Apatch] = convhull(r);
    sigmaB = sigmaB / Apatch;
  
%     [ F ] = this.computeForces();
%     [ r ] = this.computePrimalVerts();
%     sigma = zeros(size(this.d0,2),3);
%     
%     for v = 1:size(this.tri,1)
%        for ii = 1:3
%            % xx component
%            sigma(this.tri(v,ii),1) = sigma(this.tri(v,ii),1) + F(v,ii,1)*r(v,1);
%            % xy component
%            sigma(this.tri(v,ii),2) = sigma(this.tri(v,ii),2) + .5*(F(v,ii,1)*r(v,2) + F(v,ii,2)*r(v,1));
%            % yy component
%            sigma(this.tri(v,ii),3) = sigma(this.tri(v,ii),3) + F(v,ii,2)*r(v,2);
%        end
%     end
%     
%     bndryEdges = sum(abs(this.d1),1) == 1;
%     bndryCells = sum(abs(this.d0(bndryEdges,:)),1) >=1;
%     sigma(bndryCells,:) = [];
%     sigmaB = squeeze(sum(sigma,1));
%     [~,Apatch] = convhull(r);
%     sigmaB = sigmaB / Apatch;
    
end

