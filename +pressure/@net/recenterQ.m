function [ ] = recenterQ( this )
    % RE-NORMALIZE THETA

    R = (2*this.q) \ this.theta;
    this.q = bsxfun(@plus,this.q,R');
    this.theta = this.theta - 2*this.q(:,1)*R(1) - 2*this.q(:,2)*R(2) + sum(R.^2);
%     this.theta = this.theta - mean(this.theta);

end

