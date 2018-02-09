function [ q, theta, p, tri, iCells, bulkV, ERes ] = returnDual( Struct, extCell, Q )
    
    [ d0, bCells, rBX, rBY, mask, iCells ] = fitDual.AFN.getCurvedGraph( Struct, extCell );
    [ q, p, theta, tri, ~, ~, bulkV ] = fitDual.AFN.obtainStartPoint( Struct, 1 );
    q0 = [q,theta,p];

    Aeq = [[zeros(1,3*size(q,1)),ones(1,size(q,1))/size(q,1)]; ...
           [zeros(1,2*size(q,1)),ones(1,size(q,1))/size(q,1),zeros(1,size(q,1))]]; 
       
    beq = [mean(p);mean(q0(:,3))];
    d0 = sparse(d0);
    
    optimset = optimoptions('fmincon','Display','none','MaxFunEvals',1e6,'MaxIter',2e3,'TolFun',1e-6,...
               'Algorithm','interior-point','GradObj','on','GradConstr','on','DerivativeCheck','off',...
               'Hessian','on','HessFcn',@(x,lambda) hessian(x,lambda,bCells,d0,rBX,rBY,mask));
           
    lb = -inf*ones(size(q0));
    lb(:,4) = .001;
    ub = inf*ones(size(q0));
%     ub(:,4) = 1000;

    [qF,ERes] = fmincon(@(q) fitDual.AFN.energy(q, d0, bCells, rBX, rBY, mask), ...
                q0,[],[],Aeq,beq,lb,ub,@(q) ensurePositive(q,bCells,d0),optimset);
    
    q = qF(:,1:2);
    theta = qF(:,3);
    p = qF(:,4);
    
end

function [ H ] = hessian(x,lambda,bCells,d0,RBx,RBy,mask)

    Q = reshape(x,length(x)/4,4);
    NB = size(d0,1);
    NC = size(d0,2);
    
    q = Q(:,1:2);
    theta = Q(:,3);
    p = Q(:,4);
    
    % Compute terms.
%     dP = d0*p;
%     dT = d0*theta;
%     dQ = d0*q;
%     QL = sqrt(sum(dQ.^2,2));
    dP = p(bCells(:,1)) - p(bCells(:,2));
    dT = theta(bCells(:,1)) - theta(bCells(:,2));
    dQ = q(bCells(:,1),:) - q(bCells(:,2),:);
    QL = sum(dQ.^2,2);
    
    rho = bsxfun(@rdivide,d0*bsxfun(@times,p,q),dP);
%     rho = bsxfun(@rdivide,bsxfun(@times,p(bCells(:,1)),q(bCells(:,1),:)) - bsxfun(@times,p(bCells(:,2)),q(bCells(:,2),:)),dP);
    Rsq = (p(bCells(:,1)).*p(bCells(:,2)).*QL - (dP .* dT))  ./ dP.^2;
    indZ = Rsq <= 0;
    R = sqrt( Rsq );
    R(indZ) = 0;
    indZ = ~indZ;
    
    deltaX = -bsxfun(@minus,RBx,rho(:,1));
    deltaY = -bsxfun(@minus,RBy,rho(:,2));
    dMag = sqrt(deltaX.^2 + deltaY.^2);
    deltaX = deltaX ./ dMag;
    deltaY = deltaY ./ dMag;

    d = bsxfun(@minus,dMag,R);

    % Compute gradients.
    [ dRhoX, dRhoY, dR ] = fitDual.AFN.returnBondGrads( q, theta, p, bCells );
    dR(~indZ,:) = 0;
    
    %% Build objective hessian.
    % 1a + 2c
    rRatio = bsxfun(@rdivide,R,dMag);
    dNormXX = sum( rRatio .*( (deltaX .* deltaX) ) .* mask, 2);
    dNormXY = sum( rRatio .*( (deltaX .* deltaY) ) .* mask, 2);
%     dNormYX = sum( rRatio .*( (deltaY .* deltaX) ) .* mask, 2);
    dNormYY = sum( rRatio .*( (deltaY .* deltaY) ) .* mask, 2);

    dRhoX = reshape(dRhoX',[1,8,NB]);
    dRhoY = reshape(dRhoY',[1,8,NB]);
    dR = reshape(dR',[1,8,NB]);

    dRhoXT = permute(dRhoX,[2,1,3]);
    dRhoYT = permute(dRhoY,[2,1,3]);
    dRT = permute(dR,[2,1,3]);
    
    Hobj = bsxfun(@times,dRhoXT,permute(bsxfun(@times,dNormXX,permute(dRhoX,[3,1,2])),[2,3,1])) + ...
           bsxfun(@times,dRhoXT,permute(bsxfun(@times,dNormXY,permute(dRhoY,[3,1,2])),[2,3,1])) + ...
           bsxfun(@times,dRhoYT,permute(bsxfun(@times,dNormXY,permute(dRhoX,[3,1,2])),[2,3,1])) + ...
           bsxfun(@times,dRhoYT,permute(bsxfun(@times,dNormYY,permute(dRhoY,[3,1,2])),[2,3,1]));

    % 1b    
    dNormX = indZ .* sum( deltaX .* mask, 2);
    dNormY = indZ .* sum( deltaY .* mask, 2);
     
    Hobj = Hobj - bsxfun(@times,dRT,permute(bsxfun(@times,dNormX,permute(dRhoX,[3,1,2])),[2,3,1])) - ...
           bsxfun(@times,dRT,permute(bsxfun(@times,dNormY,permute(dRhoY,[3,1,2])),[2,3,1])) - ... 
           bsxfun(@times,dRhoXT,permute(bsxfun(@times,dNormX,permute(dR,[3,1,2])),[2,3,1])) - ...
           bsxfun(@times,dRhoYT,permute(bsxfun(@times,dNormY,permute(dR,[3,1,2])),[2,3,1]));
       
    % 1c
    nPix = indZ .* sum(mask,2);
    Hobj = Hobj + bsxfun(@times,dRT,permute(bsxfun(@times,nPix,permute(dR,[3,1,2])),[2,3,1]));
       
    % 2b
    avgRatio = sum( (d./dMag) .* mask, 2);
    Hobj = Hobj + bsxfun(@times,dRhoXT,permute(bsxfun(@times,avgRatio,permute(dRhoX,[3,1,2])),[2,3,1])) + ...
                  bsxfun(@times,dRhoYT,permute(bsxfun(@times,avgRatio,permute(dRhoY,[3,1,2])),[2,3,1]));
    Hobj = permute(Hobj,[3,1,2])/NB;
    
    % Return bond hessians. 2a + 2d + cons
    dNormX = sum( deltaX .* d .* mask, 2 ) / NB;
    dNormY = sum( deltaY .* d .* mask, 2 ) / NB;
    dAvg = indZ .* sum( d .* mask, 2 ) / NB;
    
    [ bondHess ] = fitDual.AFN.returnBondHessian( q, theta, p, bCells, dNormX, dNormY, dAvg, lambda.ineqnonlin ); 
    H = Hobj + bondHess;
    H = fitDual.AFN.deriv.loadHessian(H,bCells,NC,8);
    H = H + H';
    
end

function [ c, ceq, dc, dceq ] = ensurePositive(Q,bCells,d0) 
    
    q = Q(:,1:2);
    theta = Q(:,3);
    p = Q(:,4);

    ceq = [];
    c = ((d0*p) .* (d0*theta)) - p(bCells(:,1)).*p(bCells(:,2)).*sum((d0*q).^2,2);
    
    if (nargout > 2)
        dP = d0*p;
        dT = d0*theta;
        dQ = d0*q;
        QL = sum(dQ.^2,2);
        
        % X gradient
        gX = -2*(p(bCells(:,1)).*p(bCells(:,2)).*dQ(:,1));
        gX = bsxfun(@times,gX,d0);
        
        % Y gradient
        gY = -2*(p(bCells(:,1)).*p(bCells(:,2)).*dQ(:,2));
        gY = bsxfun(@times,gY,d0);

        % T gradient
        gTh = bsxfun(@times,dP,d0);
        
        % P gradient
        gP = bsxfun(@times,dT,d0);
        gP = gP - bsxfun(@rdivide,bsxfun(@times,(QL.*p(bCells(:,1)).*p(bCells(:,2))),abs(d0))',p)';
        
        dc = [gX,gY,gTh,gP]';
        dceq = [];
    end
end
