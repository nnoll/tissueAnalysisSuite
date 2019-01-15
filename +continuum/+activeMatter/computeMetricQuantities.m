function [ Lp, Ld, Ap, Ad, K, kappa ] = computeMetricQuantities( Mesh, d0, d1 )
    %COMPUTE DUAL GEO 
    
    X = Mesh.X;
    Tri = Mesh.Triangulation;
    Xd = Mesh.circumcenters;
    
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
    
    K = .5*bsxfun(@rdivide,(d0'*bsxfun(@times,d0,Ld./Lp)) * Mesh.X,Ad);
    
    % Compute Gauss curvature.
    kappa = zeros(size(Ad));
    Angles = Mesh.Angles;
    for t = 1:size(Tri,1)
        kappa(Tri(t,:)) = kappa(Tri(t,:)) + Angles(t,:)';
    end
    kappa = (2*pi - kappa)./Ad;
    
end

