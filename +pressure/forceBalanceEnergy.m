function [ E ] = forceBalanceEnergy( y, nB, x0, L0, conBonds, rB, kappa )
    % FORCE BALANCE ENERGY 
    
    x = x0 + bsxfun(@times,y,nB);
 
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
    R = sqrt(y.^2 + L0.^2);

    d = bsxfun(@minus,D,R);
    d(rB(:,1,:)==0) = 0;
  
    nPix = sum(rB(:,1,:)~=0,3);
    d = bsxfun(@rdivide,sum(d.^2,3),nPix);

    E = kappa*mean(Ev) + .5*mean(d);

end

