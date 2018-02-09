function [ Ldcp ] = decomposeTensor( this )
    %DECOMPOSE TENSOR 
    
    D = sparse(this.d0);

    rB = D*this.Mesh.X;
    rB = bsxfun(@rdivide,rB,this.Lp);
    rB = bsxfun(@times,rB,this.Ld);

    DivX = bsxfun(@times,D,rB(:,1))';
    DivY = bsxfun(@times,D,rB(:,2))';
    DivZ = bsxfun(@times,D,rB(:,3))';

    DivX = bsxfun(@ldivide,this.Ad,DivX);
    DivY = bsxfun(@ldivide,this.Ad,DivY);
    DivZ = bsxfun(@ldivide,this.Ad,DivZ);

    Ldcp = [DivX;DivY;DivZ];
    Ldcp = sparse(Ldcp);
        
end

