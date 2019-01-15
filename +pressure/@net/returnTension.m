function [ T ] = returnTension( this )

    T = this.d0 * this.q;
    T = (sum(T.^2,2));
    T = bsxfun(@times,T,(abs( ( (this.d0==1)*this.p ).* ( (this.d0==-1)*this.p ) )));
    T = T - bsxfun(@times,this.d0*this.p,this.d0*(this.theta));
    T = sqrt(T);
%     T = T/mean(T);

end

