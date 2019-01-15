function [ S ] = removeTrace( this, S )

    [ Ldcp ] = this.decomposeTensor( );
    Phi = Ldcp * [S;zeros(6,1)];
    V = length(Phi)/3;
    Phi = reshape(Phi,V,3);
    n = bsxfun(@rdivide,this.K,sqrt(sum(this.K.^2,2)));
    Phi = Phi - bsxfun(@times,dot(Phi,n,2),n);
    S = Ldcp \ Phi(:);
    S = S(1:end-6);
    
end


