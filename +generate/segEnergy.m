function [ E, dE ] = segEnergy( Q, d0, rB, mask )
    % FORCE BALANCE ENERGY 
    
    q = Q(:,1:2);
    theta = Q(:,3);
    
    % Distance from each edge assuming uniform presure everywhere
    rho = .5 * abs(d0) * q;
    qVec = d0*q; 
    dT = d0*theta;
    QL = sqrt(sum(qVec.^2,2));
    qVec = bsxfun(@rdivide,qVec,QL);
    
    deltaT = bsxfun(@minus,rB,rho);
    
    dEdge = permute(bsxfun(@times,deltaT(:,1,:),qVec(:,1)) + bsxfun(@times,deltaT(:,2,:),qVec(:,2)),[1,3,2]);
    dEdge = bsxfun(@plus,dEdge,(.5*dT)./QL);
    
    E = .5*mean( sum(mask.*(dEdge.^2),2) ./ sum(mask,2) ); 

    if (nargout > 1)
        
        NB = size(d0,1);
        gradPF = bsxfun(@rdivide,mask.*dEdge,sum(mask,2));
        IP = bsxfun(@times,deltaT(:,1,:),qVec(:,1)) + bsxfun(@times,deltaT(:,2,:),qVec(:,2));
        
        orthProj = zeros(size(deltaT));
        orthProj(:,1,:) = bsxfun(@rdivide,deltaT(:,1,:) - bsxfun(@times,IP,qVec(:,1)),QL);
        orthProj(:,2,:) = bsxfun(@rdivide,deltaT(:,2,:) - bsxfun(@times,IP,qVec(:,2)),QL);
        
        % Qx gradient
        dE_x = -.5*abs(d0') * sum(bsxfun(@times,gradPF,qVec(:,1)),2) + ...
                d0' * sum( gradPF .* permute(orthProj(:,1,:),[1,3,2]), 2) - ...
                d0' * sum( bsxfun(@times,gradPF, (dT.*qVec(:,1))./QL.^2), 2);
            
        % Qy gradient
        dE_y = -.5*abs(d0') * sum(bsxfun(@times,gradPF,qVec(:,2)),2) + ...
                d0' * sum( gradPF .* permute(orthProj(:,2,:),[1,3,2]), 2) - ...
                d0' * sum( bsxfun(@times,gradPF, (dT.*qVec(:,2))./QL.^2), 2);
            
        % Theta gradient
        dE_th = d0' * sum( bsxfun(@rdivide,.5*gradPF,QL), 2);
        
        dE = [dE_x,dE_y,dE_th]/NB;
        
    end
    
end

