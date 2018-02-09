function [ mesh ] = smoothMesh( this, nIter )
    % SMOOTH MESH 

    omega = .01;
    mesh = this;
    for ii = 1:nIter
        rB = mesh.d0*mesh.Mesh.X;
%         rB = bsxfun(@rdivide,rB,mesh.Lp);
        
        F = -mesh.d0'*rB;
        u = bsxfun(@times,dot(F,mesh.e1,2),mesh.e1) + bsxfun(@times,dot(F,mesh.e2,2),mesh.e2);
        u = omega*u;
        
        mesh = mesh.updatePrimal(u);
    end
end

