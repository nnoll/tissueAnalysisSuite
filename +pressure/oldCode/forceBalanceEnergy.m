function [ E, dE ] = forceBalanceEnergy( x, conBonds, rV1, rV2, rB, kappa )
    % FORCE BALANCE ENERGY 
    
    x21 = x(conBonds(:,2),:) - x(conBonds(:,1),:);
    x13 = x(conBonds(:,1),:) - x(conBonds(:,3),:);
    x32 = x(conBonds(:,3),:) - x(conBonds(:,2),:);
    
    L21 = sqrt(sum(x21.^2,2));
    L13 = sqrt(sum(x13.^2,2));
    L32 = sqrt(sum(x32.^2,2));
    
    S0 = (x21(:,1).*x32(:,2) - x21(:,2).*x32(:,1));
    S = sign(S0);
    S0 = S.*S0;
    
    Ev = S0 .* ((1./(L21.*L13) + 1./(L13.*L32) + 1./(L32.*L21)));
    Ev = Ev / 3;
    
    D = bsxfun(@minus,rB,x);
    D = sqrt(sum(D.^2,2));
    R = .5*( sqrt(sum( (rV1-x).^2, 2) ) + sqrt(sum( (rV2-x).^2, 2) ) );
    
    d = bsxfun(@minus,D,R);
    d(rB(:,1,:)==0) = 0;
    nPix = sum(rB(:,1,:)~=0,3);

    d = bsxfun(@rdivide,sum(d.^2,3),nPix);
    E = mean(Ev) + .5*kappa*mean(d);

    if (nargout == 2)
        
        x1 = x(conBonds(:,1),:);
        x2 = x(conBonds(:,2),:);
        x3 = x(conBonds(:,3),:);

        EvG = zeros(size(conBonds,1),2,3);
        Rot = [0,-1;1,0];
        
        t21 = bsxfun(@times,S,x21);
        t32 = bsxfun(@times,S,x32);
        t13 = bsxfun(@times,S,x13);
        
        % Contribution from first term of energy.
        
        EvG(:,:,1) = bsxfun(@times,Ev./S0,t32*Rot');
        EvG(:,:,1) = EvG(:,:,1) - bsxfun(@times,(S0./(3*L21.^3)) .* (1./L13 + 1./L32),x1-x2);
        EvG(:,:,1) = EvG(:,:,1) - bsxfun(@times,(S0./(3*L13.^3)) .* (1./L32 + 1./L21),x1-x3);

        EvG(:,:,2) = bsxfun(@times,Ev./S0,t13*Rot');
        EvG(:,:,2) = EvG(:,:,2) - bsxfun(@times,(S0./(3*L21.^3)) .* (1./L13 + 1./L32),x2-x1);
        EvG(:,:,2) = EvG(:,:,2) - bsxfun(@times,(S0./(3*L32.^3)) .* (1./L13 + 1./L21),x2-x3);
        
        EvG(:,:,3) = bsxfun(@times,Ev./S0,t21*Rot');
        EvG(:,:,3) = EvG(:,:,3) - bsxfun(@times,(S0./(3*L32.^3)) .* (1./L13 + 1./L21),x3-x2);
        EvG(:,:,3) = EvG(:,:,3) - bsxfun(@times,(S0./(3*L13.^3)) .* (1./L32 + 1./L21),x3-x1);
        
        % Sum up over vertices.
        dEg = zeros(size(conBonds,1),2*size(x,1));
        for ii = 1:3
            for jj = 1:2
                dEg( (1:size(conBonds,1)) + size(conBonds,1)*( size(x,1)*(jj-1) + conBonds(:,ii) - 1 )' ) = EvG(:,jj,ii);
            end
        end
        dE = sum(dEg,1)'/size(conBonds,1);
        
        % Add in contribution from second term.
        dG = bsxfun(@minus,D,R);
        dG = (rB(:,1,:)~=0) .* dG;
        nPix = sum(rB(:,1,:)~=0,3);
        denom = D;
        denom(denom==0) = 1;
        
        dD = bsxfun(@rdivide,bsxfun(@minus,x,rB),denom);
        
        dR1 = (x - rV1);
        dR1 = bsxfun(@rdivide,dR1,sqrt(sum(dR1.^2,2)));
        dR2 = (x - rV2);
        dR2 = bsxfun(@rdivide,dR2,sqrt(sum(dR2.^2,2)));
        dR = .5*(dR1 + dR2);
        
        dVec = bsxfun(@minus,dD,dR);
        dE1 = sum(bsxfun(@times,dG,dVec(:,1,:)),3)./nPix;
        dE2 = sum(bsxfun(@times,dG,dVec(:,2,:)),3)./nPix;
        dE = dE + (kappa*[dE1;dE2])/size(x,1);
        
    end

end

