function [ Lp, Ld, Ap, Ad, K ] = computeMetricQuantities( this )
    %COMPUTE DUAL GEO 
    
    Mesh = this.Mesh;
    d0 = this.d0;
    d1 = this.d1;
    
    X = Mesh.Points;
    Tri = Mesh.ConnectivityList;
    Xd = Mesh.circumcenter;
    
    % Compute primal geometric quanties.
    Rb = d0*X;
    Lp = sqrt(sum(Rb.^2,2));
    
    R1 = X(Tri(:,2),:) - X(Tri(:,1),:);
    R2 = X(Tri(:,3),:) - X(Tri(:,1),:);

    Ap = .5 * sqrt(sum(cross(R1,R2).^2,2));
        
    % Compute dual bond vectors and their lengths.
    Xb = d1' * Xd;
    Ld = sqrt(sum(Xb.^2,2));
    
    % Compute cell areas.
    Ad = .25*abs(d0')*(Ld.*Lp);
    
    K = .5*bsxfun(@rdivide,(d0'*bsxfun(@times,d0,Ld./Lp)) * X,Ad);
   
    
end

