function [ E ] = segEnergy( x, bVerts, rB, mask )
    % FORCE BALANCE ENERGY 
    
    nB = x(bVerts(:,1),:) - x(bVerts(:,2),:);
    D = sqrt(sum(nB.^2,2));

    nB = bsxfun(@rdivide,nB,D);
    nB = nB*[0,-1;1,0];
    
    x0 = .5*(x(bVerts(:,1),:) + x(bVerts(:,2),:));
    
    delta = bsxfun(@minus,rB,x0);
    
    IP = permute(bsxfun(@times,delta(:,1,:),nB(:,1)) + bsxfun(@times,delta(:,2,:),nB(:,2)),[1,3,2]);
    IP(mask==0) = 0;

    dMag = permute(sum(delta.^2,2),[1,3,2]);
    
    A = 2*sum(IP.^2,2);
    B = sum(  bsxfun(@minus,dMag,(D/2).^2) .* IP ,2 );
    
    y = B./A;
    R = (y.^2 + (D/2).^2);
   
    rho = x0+bsxfun(@times,y,nB);
    delta = bsxfun(@minus,rB,rho);

    dMag = (permute(sum(delta.^2,2),[1,3,2]));
    dEdge = bsxfun(@minus,dMag,R).^2;
    E = mean( sum(mask.*dEdge,2) ./ sum(mask,2) ); 

end

