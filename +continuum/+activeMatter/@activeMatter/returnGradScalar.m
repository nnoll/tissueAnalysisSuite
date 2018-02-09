function [ gradS ] = returnGradScalar( this )
    %RETURN GRAD SCALAR
    gradS = zeros(3*size(this.Mesh.X,1),size(this.Mesh.X,1));
    
    rHat = this.d0*this.Mesh.X;
    rHat = bsxfun(@rdivide,rHat,sqrt(sum(rHat.^2,2)));
    NV = size(this.Mesh.X,1);
    n = bsxfun(@rdivide,this.K,sqrt(sum(this.K.^2,2)));
    for v = 1:NV
        conEdges = find(this.d0(:,v));
        for e = conEdges'
            rHatPrime = rHat(e,:) - bsxfun(@times,dot(rHat(e,:),n(v,:),2),n(v,:));
            rHatPrime = bsxfun(@rdivide,rHatPrime,sqrt(sum(rHatPrime.^2,2)));
            vE = find(this.d0(e,:));
            vE = vE(vE~=v);
            
            % X components
            gradS(v,v) = gradS(v,v) + this.d0(e,v) * rHatPrime(:,1) * this.Ld(e) / (2*this.Ad(v));
            gradS(v,vE) = gradS(v,vE) + this.d0(e,vE) * rHatPrime(:,1) * this.Ld(e) / (2*this.Ad(v));

            % Y components
            gradS(v+NV,v) = gradS(v+NV,v) + this.d0(e,v) * rHatPrime(:,2) * this.Ld(e) / (2*this.Ad(v));
            gradS(v+NV,vE) = gradS(v+NV,vE) + this.d0(e,vE) * rHatPrime(:,2) * this.Ld(e) / (2*this.Ad(v));
            
            % Z components
            gradS(v+2*NV,v) = gradS(v+2*NV,v) + this.d0(e,v) * rHatPrime(:,3) * this.Ld(e) / (2*this.Ad(v));
            gradS(v+2*NV,vE) = gradS(v+2*NV,vE) + this.d0(e,vE) * rHatPrime(:,3) * this.Ld(e) / (2*this.Ad(v));
        end
    end
    
    E1 = [diag(this.e1(:,1)),diag(this.e1(:,2)),diag(this.e1(:,3))];
    E2 = [diag(this.e2(:,1)),diag(this.e2(:,2)),diag(this.e2(:,3))];

    N = sparse([E1;E2]);
    gradS = N*gradS;

end

