function [ E, dE ] = cSegEnergy( Q, d0, bCells, rB, mask )
    % FORCE BALANCE ENERGY 
    
    q = Q(:,1:2);
    theta = Q(:,3);
    p = Q(:,4);
    
    % Distance from each edge assuming finite presure everywhere
    rho = bsxfun(@rdivide,d0*bsxfun(@times,p,q),d0*p);
    R = sqrt(  p(bCells(:,1)).*p(bCells(:,2)).*sum((d0*q).^2,2) - ((d0*p) .* (d0*theta)) ) ./ (d0*p);
    R = abs(R);

    delta = bsxfun(@minus,rB,rho);
    dMag = sqrt(permute( sum(delta.^2,2), [1,3,2]));
    dEdge = bsxfun(@minus,dMag,R).^2;
    
    E = .5*mean( sum(mask.*dEdge,2) ./ sum(mask,2) ); 

    if (nargout == 2)
        
        % Compute pre-factor of each edge.
        gradPF = bsxfun(@minus,dMag,R);
        gradPF = gradPF .* mask;
        gradPF = bsxfun(@rdivide,gradPF,sum(mask,2));
        
        % Compute both needed linear operators
        Lbare = d0';
        Lp = bsxfun(@times,p,d0');
        Lsym = bsxfun(@rdivide,abs(d0)',p);
        Lx = bsxfun(@times,q(:,1),d0');
        Ly = bsxfun(@times,q(:,2),d0');
        
        % Normalized radial vector from each data point.
        radDir = permute(bsxfun(@minus,rB,rho),[1,3,2]);
        radDir = bsxfun(@rdivide,radDir,sqrt(sum(radDir.^2,3)));
        radDir(:,:,1) = radDir(:,:,1) .* mask;
        radDir(:,:,2) = radDir(:,:,2) .* mask;
        
        % Compute needed bond geometric quantities.
        dP = d0*p;
        dT = d0*theta;
        dQ = d0*q;
        QL = sum(dQ.^2,2);
        NB = size(d0,1);
        
        % grad Qx
        dE_x = Lp*sum((gradPF .* bsxfun(@rdivide,-radDir(:,:,1),dP)),2) - ...
               Lbare*sum(bsxfun(@times,gradPF,(p(bCells(:,1)).*p(bCells(:,2)).*dQ(:,1))./(dP.^2.*R)),2);
        
        % grad Qy
        dE_y = Lp*sum((gradPF .* bsxfun(@rdivide,-radDir(:,:,2),dP)),2) - ...
               Lbare*sum(bsxfun(@times,gradPF,(p(bCells(:,1)).*p(bCells(:,2)).*dQ(:,2))./(dP.^2.*R)),2);
           
        % grad Theta
        dE_t = Lbare*sum(bsxfun(@times,gradPF,1./(2*R.*dP)),2);
        
        % grad p
        dE_p = Lbare*sum( gradPF .* ( bsxfun(@times,radDir(:,:,1),rho(:,1)./dP) + bsxfun(@times,radDir(:,:,2),rho(:,2)./dP) ),2) + ...
               Lbare*sum( bsxfun(@times,gradPF, QL.*p(bCells(:,1)).*p(bCells(:,2))./(R.*dP.^3)) ,2) - ...
               Lbare*sum( bsxfun(@times,gradPF, dT./(2*R.*dP.^2)) ,2) - ...
               Lsym*sum( bsxfun(@times,gradPF, (QL.*p(bCells(:,1)).*p(bCells(:,2)))./(2*R.*dP.^2)) ,2) - ...
               Lx*sum( bsxfun(@times,gradPF.*radDir(:,:,1),1./dP) ,2) - ...
               Ly*sum( bsxfun(@times,gradPF.*radDir(:,:,2),1./dP) ,2); %abs(Lp)*sum( bsxfun(@times,gradPF, QL./(2*R.*dP.^2)) ,2) - ...

        
        dE = [dE_x,dE_y,dE_t,dE_p]/NB;
    end
end

