function [ rho, R, edgeVerts, edgeCells, badVerts ] = computeEdgeArcs( this, mode )
    % COMPUTE EDGE ARCS
    
    edgeCells = zeros(size(this.d0,1),2);
    edgeVerts = zeros(size(this.d0,1),2);

    badVerts = [];
    for b = 1:size(this.d0,1)
        edgeCells(b,1) = find(this.d0(b,:)==1);
        edgeCells(b,2) = find(this.d0(b,:)==-1);
        if (sum(abs(this.d1(:,b))) == 2)
            edgeVerts(b,1) = find(this.d1(:,b)==1);
            edgeVerts(b,2) = find(this.d1(:,b)==-1);
        else
            edgeVerts(b,1) = find(this.d1(:,b));
            badVerts = [badVerts,b];
        end
    end
    
    
    rho = bsxfun(@times,this.p(edgeCells(:,1)),this.q(edgeCells(:,1),:)) ...
        - bsxfun(@times,this.p(edgeCells(:,2)),this.q(edgeCells(:,2),:));
    rho = bsxfun(@rdivide,rho,this.p(edgeCells(:,1))-this.p(edgeCells(:,2)));

%     R = sum(rho.^2,2) + bsxfun(@rdivide,this.p(edgeCells(:,2)).*(sum(this.q(edgeCells(:,2),:).^2,2)-this.theta(edgeCells(:,2))) - ...
%         this.p(edgeCells(:,1)).*(sum(this.q(edgeCells(:,1),:).^2,2)-this.theta(edgeCells(:,1))),this.p(edgeCells(:,1))-this.p(edgeCells(:,2)));
    R = this.p(edgeCells(:,1)).*this.p(edgeCells(:,2)) .* sum( (this.d0 * this.q).^2,2) - (this.d0 * this.p) .* (this.d0 * (this.theta));
    R = R ./ (this.d0*this.p).^2;
    R = real(sqrt(R));

    if (nargin == 1 || mode == 0)
        edgeVerts(badVerts,:) = [];
        R(badVerts) = [];
        rho(badVerts,:) = [];
    end
    
end

