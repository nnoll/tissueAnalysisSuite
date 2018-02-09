function [ S ] = traceToEdgeTensor( this, Strace )
%TRACETOEDGETENSOR 

    tr = abs(this.d0');
    tr = bsxfun(@ldivide,this.Ad,tr);
    tr = bsxfun(@times,tr',.5*this.Lp.*this.Ld)';
    
    Div = this.returnDiv();
    Div = bsxfun(@ldivide,[this.Ad;this.Ad],Div);
    
    rB = this.d0*this.Mesh.X;
    rB = bsxfun(@times,rB,this.Ld./this.Lp);

    Dx = bsxfun(@times,this.d0,rB(:,1));
    Dy = bsxfun(@times,this.d0,rB(:,2));
    Dz = bsxfun(@times,this.d0,rB(:,3));

    D1 = bsxfun(@times,Dx,this.e1(:,1)') + bsxfun(@times,Dy,this.e1(:,2)') + bsxfun(@times,Dz,this.e1(:,3)');
    D2 = bsxfun(@times,Dx,this.e2(:,1)') + bsxfun(@times,Dy,this.e2(:,2)') + bsxfun(@times,Dz,this.e2(:,3)');
    
    D = [D1,D2];
    vDiv = .5*abs(this.d0)'*D;
    vDiv = bsxfun(@ldivide,this.Ad,vDiv);
    
    vGrad = -vDiv';
    L = [tr;Div];
    
    Phi = [Strace;vGrad*Strace];
    S = L \ Phi;
    
end

