function [ L, Div, r1, r2, r3, tx, ty, tz, N ] = compute3D_displacements( this )
    % COMPUTE DISPLACEMENTS 
    
    Div = this.decomposeTensor;
    
    D = this.d0;
    rB = D*this.Mesh.X;
    rB = bsxfun(@rdivide,rB,this.Lp.^2);

    GradX = bsxfun(@times,D,rB(:,1));
    GradY = bsxfun(@times,D,rB(:,2));
    GradZ = bsxfun(@times,D,rB(:,3));
    
    Grad = [GradX,GradY,GradZ];
    
    L = Div*Grad;
    n = bsxfun(@rdivide,this.K,sqrt(sum(this.K.^2,2)));
    
    N = [diag(n(:,1)),diag(n(:,2)),diag(n(:,3))];
    
    r1 = [this.Mesh.X(:,2),-this.Mesh.X(:,1),zeros(size(this.Mesh.X(:,3)))];
    r2 = [zeros(size(this.Mesh.X(:,2))),-this.Mesh.X(:,3),this.Mesh.X(:,2)];
    r3 = [-this.Mesh.X(:,3),zeros(size(this.Mesh.X(:,2))),this.Mesh.X(:,1)];

    r1 = r1/sum(this.Ap);
    r2 = r2/sum(this.Ap);
    r3 = r3/sum(this.Ap);
        
    tx = [ones(1,size(this.Mesh.X,1)),zeros(1,size(this.Mesh.X,1)),zeros(1,size(this.Mesh.X,1))];
    ty = [zeros(1,size(this.Mesh.X,1)),ones(1,size(this.Mesh.X,1)),zeros(1,size(this.Mesh.X,1))];
    tz = [zeros(1,size(this.Mesh.X,1)),zeros(1,size(this.Mesh.X,1)),ones(1,size(this.Mesh.X,1))];

    tx = tx/sum(this.Ap);
    ty = ty/sum(this.Ap);
    tz = tz/sum(this.Ap);
    
end

