function [ L, Lbulk, Div, r1, r2, r3 ] = computeOperators( this )
    % COMPUTE DISPLACEMENTS 
    
    Div = this.returnDiv();
    Grad = this.returnGrad();
    Div = bsxfun(@ldivide,[this.Ad;this.Ad],Div);
    L = Div*Grad;

    % While L has full rank, the rotation modes have eigenvalues at around
    % 1e-6 which dominates the inverse. Thus we should gap these modes
    % EXPLICITLY!

    r1 = [this.Mesh.X(:,2),-this.Mesh.X(:,1),zeros(size(this.Mesh.X(:,3)))];
    r2 = [zeros(size(this.Mesh.X(:,2))),-this.Mesh.X(:,3),this.Mesh.X(:,2)];
    r3 = [-this.Mesh.X(:,3),zeros(size(this.Mesh.X(:,2))),this.Mesh.X(:,1)];

    r1 = r1/sum(this.Ap);
    r2 = r2/sum(this.Ap);
    r3 = r3/sum(this.Ap);
    
    r1 = [dot(r1,this.e1,2)',dot(r1,this.e2,2)'];
    r2 = [dot(r2,this.e1,2)',dot(r2,this.e2,2)'];
    r3 = [dot(r3,this.e1,2)',dot(r3,this.e2,2)'];

    [ gradS ] = this.returnGradScalar();
    vDiv = -gradS';
    Lbulk = gradS*vDiv;
    
%     rB = this.d0*this.Mesh.X;
%     rB = bsxfun(@times,rB,this.Ld./this.Lp);
% 
%     Dx = bsxfun(@times,this.d0,rB(:,1));
%     Dy = bsxfun(@times,this.d0,rB(:,2));
%     Dz = bsxfun(@times,this.d0,rB(:,3));
% 
%     D1 = bsxfun(@times,Dx,this.e1(:,1)') + bsxfun(@times,Dy,this.e1(:,2)') + bsxfun(@times,Dz,this.e1(:,3)');
%     D2 = bsxfun(@times,Dx,this.e2(:,1)') + bsxfun(@times,Dy,this.e2(:,2)') + bsxfun(@times,Dz,this.e2(:,3)');
% 
%     D = [D1,D2];
%     
%     vDiv = .5*abs(this.d0)'*D;
%     vDiv = bsxfun(@ldivide,this.Ad,vDiv);
%     vGrad = -vDiv';
% 
%     Lbulk = vGrad*vDiv;      
%     
%     xB = (this.d1')*this.Mesh.circumcenters;
% 
%     Dx = bsxfun(@times,this.d0,xB(:,1));
%     Dy = bsxfun(@times,this.d0,xB(:,2));
%     Dz = bsxfun(@times,this.d0,xB(:,3));
% 
%     D1 = bsxfun(@times,Dx,this.e1(:,1)') + bsxfun(@times,Dy,this.e1(:,2)') + bsxfun(@times,Dz,this.e1(:,3)');
%     D2 = bsxfun(@times,Dx,this.e2(:,1)') + bsxfun(@times,Dy,this.e2(:,2)') + bsxfun(@times,Dz,this.e2(:,3)');
% 
%     D = [D1,D2];
%     
%     vDiv = .5*abs(this.d0)'*D;
%     vDiv = bsxfun(@ldivide,this.Ad,vDiv);
%     vGrad = -vDiv';
% 
%     LpBulk = vGrad*vDiv;      
end

