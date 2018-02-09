function [ PNnew ] = rescaleDual( this, lambda, R )

    % RE-NORMALIZE THETA
    theta = lambda*(this.theta + (1-lambda)*sum(this.q.^2,2));
    q = lambda*this.q;
    q = bsxfun(@plus,q,R);
    PNnew = pressure.net(q,theta,this.p,this.tri,this.cellLabels);
    
end

