function [ q, p, theta, tri, bulk0, ext0, bulkV, ERes ] = obtainStartPoint( Struct, extCell, q, p, theta )
    % FIT PRESSURE 

    % Obtain initial conditions from ATN fit.
    if (nargin == 2)
        [ q, p, tri, bulk0, ext0, bulkV, ERes ] = fitDual.AFN.returnReducedDual( Struct, extCell );
    else
        [ tri, bulk0, ext0 ] = fitDual.returnGraph( Struct, extCell );
    end
    
    if (std(p) == 0)
        p = p + .01*randn(size(p));
    end
    
    if (nargin <= 4)
        [ x0, q, p ] = fitDual.AFN.fitThetaModes( q, p, Struct, extCell );
    else
        x0 = theta;
    end
    
    pTri = fitDual.AFN.returnPTri(q,p,tri);
    
    [ ~, ~, ~, ~, ~, r0 ] = fitDual.returnGraph( Struct, extCell );
    [ d0, bCells, rBX, rBY, mask ] = fitDual.AFN.getCurvedGraph( Struct, extCell );

    energyFunc = @(x) fitDual.AFN.thetaEnergy(x, q, p, d0, bCells, rBX, rBY, mask);

%     nonlinFunc = @(x) fitDual.AFN.thetaCon( x, q, p, pTri );
       
%     max(nonlinFunc(zeros(size(x0))))
%     max(nonlinFunc(x0))
%     pause
    
    if (energyFunc(zeros(size(x0))) < energyFunc(x0))
        x0 = zeros(size(x0));
    end

    Aeq = ones(1,length(x0));
    beq = 0;
    
    dP = p(bCells(:,1)) - p(bCells(:,2));
    A = bsxfun(@times,dP,d0);
    dQ = q(bCells(:,1),:) - q(bCells(:,2),:);
    QL = sum(dQ.^2,2);
    b = p(bCells(:,1)).*p(bCells(:,2)).*QL;
    optimset = optimoptions('fmincon','Display','none','Algorithm','interior-point','TolFun',1e-6, ... 
                            'MaxFunEvals',5e6,'MaxIter',1e3,'GradObj','on','GradConstr','on','Hessian','lbfgs',...
                            'HessFcn',@(x,lambda) hessian(x,lambda,q,p,tri,pTri,r0,bCells,d0,rBX,rBY,mask),...
                            'DerivativeCheck','off');
    [theta,ERes] = fmincon(energyFunc,x0,A,b,Aeq,beq,[],[],[],optimset);
    
end

function [ H ] = hessian(theta,lambda,q,p,tri,pTri,rV,bCells,d0,RBx,RBy,mask)

    %% Compute hessian of objective.
    NB = size(d0,1);
    NC = size(d0,2);
%     NV = size(pTri,1);

    % Compute terms.
%     dP = d0*p;
%     dT = d0*theta;
%     dQ = d0*q;
%     QL = sqrt(sum(dQ.^2,2));
    dP = p(bCells(:,1))-p(bCells(:,2));
    dT = theta(bCells(:,1))-theta(bCells(:,2));
    dQ = q(bCells(:,1),:)-q(bCells(:,2),:);
    QL = sum(dQ.^2,2);
    
    rho = bsxfun(@rdivide,d0*bsxfun(@times,p,q),dP);
    Rsq = (p(bCells(:,1)).*p(bCells(:,2)).*QL - (dP .* dT))  ./ dP.^2;
    indZ = Rsq <= 0;
    R = sqrt( Rsq );
    R(indZ) = 0;
%     indZ = ~indZ;
    
    deltaX = -bsxfun(@minus,RBx,rho(:,1));
    deltaY = -bsxfun(@minus,RBy,rho(:,2));
    dMag = sqrt(deltaX.^2 + deltaY.^2);

    d = bsxfun(@minus,dMag,R);

    % Compute gradients.
    [ dR ] = fitDual.AFN.returnBondGradsTheta( q, theta, p, bCells );
    dR(indZ,:) = 0;
%     dR = bsxfun(@times,indZ,dR);
    
    %% Build objective hessian.
    nPix = sum(mask,2);
    Hobj = bsxfun(@times,reshape(dR',[2,1,NB]),reshape(bsxfun(@times,nPix,dR)',[1,2,NB]));
    Hobj = permute(Hobj,[3,1,2]);
    
    % Return bond hessians. 2a + 2d
    dAvg = sum( d .* mask, 2 );
    [ hR ] = fitDual.AFN.returnBondHessianTheta( q, theta, p, bCells );
    hR(indZ,:,:) = 0;
    Hobj = real(Hobj + bsxfun(@times,hR,dAvg))/NB;
    Hobj = fitDual.AFN.deriv.loadHessian(Hobj, bCells, NC, 2);
    H = Hobj + Hobj';
    
    %% Compute Hessian of each constraint.
%     row = zeros(9*NV,1);
%     col = zeros(9*NV,1);
%     value = zeros(9*NV,1);
%     
%     [ r, indB, r0, lambdaV ] = fitDual.AFN.returnVertexPositions( q, theta, p, tri, pTri, rV );
%     NV = size(r,1);
%     
%     [ dR, dR0, r0, etaStar ] = fitDual.AFN.returnThetaGrad( r, q, p, pTri, indB, r0, lambdaV );
% %     [ dR0 ] = fitDual.AFN.returnThetaGrad( r, q, p, pTri );
% %     [ r0, ~, etaStar ] = fitDual.AFN.returnR0( q, p, theta, pTri );
%     
%     r0Mag = sqrt(dot(r0,r0,2));
%     etaMag = sqrt(sum(etaStar.^2,2));
%     
%     n = 1;
%     for ii = 1:3
%         for jj = 1:3
%             II = ((n-1)*NV) + (1:NV);
%             row(II) = pTri(:,ii);
%             col(II) = pTri(:,jj);
%             IPi = bsxfun(@times,r0(:,1),dR(:,1,ii)) + bsxfun(@times,r0(:,2),dR(:,2,ii));
%             IPj = bsxfun(@times,r0(:,1),dR(:,1,jj)) + bsxfun(@times,r0(:,2),dR(:,2,jj));
%             value(II) = lambda.ineqnonlin.*((dR0(:,1,ii) .* dR0(:,1,jj)) + (dR0(:,2,ii) .* dR0(:,2,jj)) - (IPi.*IPj ./ r0Mag)) .* (etaMag./r0Mag);
%             n = n + 1;
%         end
%     end
%     
%     Hcon = sparse(row,col,value,NC,NC);
%     Hcon = Hcon + Hcon';
%     H = Hobj + Hcon;
    
end

