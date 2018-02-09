function [ F ] = computeForces( this )
    % COMPUTE FORCES 
    
    F = zeros(size(this.tri,1),3,2); % 2 component force vector for every vertex/edge combo.
    rot = [0,1;-1,0];
    
    % Compute edge centroids from our dual parameters. 
    [ rho, R, edgeDualVerts, edgeVerts ] = this.computeEdgeArcs( 1 );
    
    [ r ] = this.computePrimalVerts();
    [ T ] = this.returnTension();
    
%     Tscale = this.d0*this.q;
%     Tscale = mean(sqrt(sum(Tscale.^2,2)));
    
    for e = 1:size(rho,1)
        if (~isinf(R(e)) && ~isnan(R(e)))
            F1 = r(edgeDualVerts(e,1),:) - rho(e,:);
            F2 = rho(e,:) - r(edgeDualVerts(e,2),:);
            F1 = F1*rot;
            F2 = F2*rot;
            
            F1 = T(e) * F1 ;
            F2 = T(e) * F2 ;
            ind1 = ~ismember(this.tri(edgeDualVerts(e,1),:) ,edgeVerts(e,:));
            ind2 = ~ismember(this.tri(edgeDualVerts(e,2),:) ,edgeVerts(e,:));
            F(edgeDualVerts(e,1),ind1,:) = F1;
            F(edgeDualVerts(e,2),ind2,:) = F2;
        else
            F1 = this.q(edgeVerts(e,2),:) - this.q(edgeVerts(e,1),:);
            F2 = this.q(edgeVerts(e,1),:) - this.q(edgeVerts(e,2),:);
            F1 = F1*rot;
            F2 = F2*rot;
            S = this.d0(e,edgeVerts(e,1));
            
            ind1 = ~ismember(this.tri(edgeDualVerts(e,1),:),edgeVerts(e,:));
            F(edgeDualVerts(e,1),ind1,:) = -S*this.d1(edgeDualVerts(e,1),e)*F1;
            if (edgeDualVerts(e,2) > 0)
                ind2 = ~ismember(this.tri(edgeDualVerts(e,2),:) ,edgeVerts(e,:));
                F(edgeDualVerts(e,2),ind2,:) = S*this.d1(edgeDualVerts(e,2),e)*F2;
            end
        end
    end

end

