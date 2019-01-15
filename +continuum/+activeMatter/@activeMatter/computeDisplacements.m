function [ u, U, T, Div ] = computeDisplacements( this, S, B )
    % COMPUTE DISPLACEMENTS 
    
    if (nargin == 2)
        B = 0;
    end
    
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

    % Compute the u's given our constraint.     
    if (B ~= 0) % Apply bulk modulus.

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

        LbulkM = B*vGrad*vDiv;
        L = LbulkM + L;
        L = [L;r1;r2;r3];
%             Lflat = [Lflat;r1;r2;r3];

%             Strace = this.computeTrace(S);
        Fact = Div*S; % + B*vGrad*Strace;

        u = sparse(L) \ [Fact;0;0;0];
        U = Grad*u;
        Utrace = this.computeTrace(U);
        bulkTrace = (Utrace); % - Strace);

        T = (U - S) + B*this.traceToEdgeTensor(bulkTrace);

    else
        L = sparse([L;r1;r2;r3]);
        Fact = Div*S;
        u = L \ [Fact;0;0;0];
        U = Grad*u;
        T = U - S;
    end

    u = bsxfun(@times,u(1:(length(u)/2)),this.e1) + bsxfun(@times,u((length(u)/2) + (1:(length(u)/2))),this.e2);
        
    
end

