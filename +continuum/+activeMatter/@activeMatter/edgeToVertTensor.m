function [ L ] = edgeToVertTensor( this )
    
    rB = this.d0 * this.Mesh.X;
    rB = bsxfun(@rdivide,rB,this.Lp);
    
    Dx = bsxfun(@times,this.d0,rB(:,1));
    Dy = bsxfun(@times,this.d0,rB(:,2));
    Dz = bsxfun(@times,this.d0,rB(:,3));

    D1 = bsxfun(@times,Dx',this.e1(:,1)) + bsxfun(@times,Dy',this.e1(:,2)) + bsxfun(@times,Dz',this.e1(:,3));
    D2 = bsxfun(@times,Dx',this.e2(:,1)) + bsxfun(@times,Dy',this.e2(:,2)) + bsxfun(@times,Dz',this.e2(:,3));

    D1 = bsxfun(@times,D1,sqrt(.5*(this.Lp.*this.Ld)'));
    D2 = bsxfun(@times,D2,sqrt(.5*(this.Lp.*this.Ld)'));

    D1 = bsxfun(@ldivide,sqrt(this.Ad),D1);
    D2 = bsxfun(@ldivide,sqrt(this.Ad),D2);
    
    L = sparse([D1.*D1;D1.*D2;D2.*D2]);
    
end

