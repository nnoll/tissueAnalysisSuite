function [ vP ] = predictVelocity( this, B, sigma, time_series )
    % PREDICT VELOCITY 
   
   % Compute Divergence of Tensor Sigma
   Z = time_series(1).charts{1};
   Phi = time_series(1).charts{2};
   f = this.computeDiv(sigma,Z,Phi);
    
   V = size(this.emb.Mesh.X,1);
   P1_DoFmap = uint32(this.emb.Mesh.Triangulation);
   FEM = inferV([],this.emb.Mesh.X,uint32(this.emb.Mesh.Triangulation),[],[],P1_DoFmap,f);
   
   delV = FEM(1).MAT;
   gradDiv = FEM(2).MAT;
   curv = zeros(size(gradDiv));
   
   kappa = 1.1*(this.emb.kappa.*this.emb.Ad);
   curv((1:V) + 3*V*(0:(V-1))) = kappa;
   curv(((V+1):(2*V)) + 3*V*(V:(2*V-1))) = kappa;
   curv(((2*V+1):(3*V)) + 3*V*(2*V:(3*V-1))) = kappa;
   
   L = delV + (B*gradDiv) - sparse(curv);
   f = FEM(3).MAT;

   % Project onto tangential plane.
   E1 = [diag(this.emb.e1(:,1)),diag(this.emb.e1(:,2)),diag(this.emb.e1(:,3))];
   E2 = [diag(this.emb.e2(:,1)),diag(this.emb.e2(:,2)),diag(this.emb.e2(:,3))];
   
   N = sparse([E1;E2]);
   L = N*L*N';
   f = N*f;
   
   v = L \ f;
   vP = bsxfun(@times,v(1:V),this.emb.e1) + bsxfun(@times,v(V+(1:V)),this.emb.e2);
   
end

