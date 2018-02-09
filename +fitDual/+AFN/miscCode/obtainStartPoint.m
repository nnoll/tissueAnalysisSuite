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
    
    if (nargin <= 2)
        [ x0, q, p ] = fitDual.AFN.fitThetaModes( q, p, Struct, extCell );
    else
        x0 = theta;
    end
    
    pTri = fitDual.AFN.returnPTri(q,p,tri);
    
    [ ~, ~, ~, ~, ~, r0 ] = fitDual.returnGraph( Struct, extCell );
    dim = size(Struct.labelMat);
    [ d0, ~, ~, ~, iCells ] = fitDual.ATN.computeDiffOperators( Struct, extCell );  
    [ involvedBonds ] = generate.bondMap(Struct);
    involvedBonds = involvedBonds{1};

    z = 0;
    bCells = zeros(length(involvedBonds),2);
    badBonds = [];
    for b = 1:length(involvedBonds)
        if (involvedBonds(b) > 0)
%             c1 = find(iCells == Struct.Bdat(involvedBonds(b)).cells(1));
%             c2 = find(iCells == Struct.Bdat(involvedBonds(b)).cells(2));
            c1 = find(d0(b,:)==1);
            c2 = find(d0(b,:)==-1);
            if (length(Struct.Bdat(involvedBonds(b)).pix) > z)
                z = length(Struct.Bdat(involvedBonds(b)).pix);
            end
            bCells(b,:) = [c1,c2];
        else
            badBonds = [badBonds,b];
        end
    end
    
    involvedBonds(badBonds) = [];
    bCells(badBonds,:) = [];
    d0(badBonds,:) = [];

    rBX = zeros(length(involvedBonds),z+2);
    rBY = zeros(length(involvedBonds),z+2);
    for ii = 1:length(involvedBonds)
        ind = Struct.Bdat(involvedBonds(ii)).pix;
        [y,x] = ind2sub(dim,ind);
        x = double(x);
        y = double(y);
        v1 = Struct.Bdat(involvedBonds(ii)).verts(1);
        v2 = Struct.Bdat(involvedBonds(ii)).verts(2);
        rBX(ii,1:(length(x)+2)) = [double(x);Struct.Vdat(v1).vertxcoord;Struct.Vdat(v2).vertxcoord];
        rBY(ii,1:(length(x)+2)) = [double(y);Struct.Vdat(v1).vertycoord;Struct.Vdat(v2).vertycoord];
    end

    mask = (rBX ~= 0);
%     energyFunc = @(x) fitDual.AFN.thetaEnergy(x, q, p, r0, tri, pTri);
    energyFunc = @(x) fitDual.AFN.thetaEnergy(x, q, p, d0, bCells, rBX, rBY, mask);

%     if (energyFunc(x0) > energyFunc(zeros(size(x0))))
%         x0 = zeros(size(x0));
%     end

    nonlinFunc = @(x) fitDual.AFN.thetaCon( x, q, p, pTri );
       
    Aeq = ones(1,length(x0));
    beq = 0;
    
    optimset = optimoptions('fmincon','Display','iter','Algorithm','interior-point','TolFun',1e-6, ... 
                            'MaxFunEvals',5e6,'MaxIter',3e3,'GradObj','on','GradConstr','on','Hessian','on',...
                            'HessFcn',@(x,lambda) hessian(x,lambda,q,p,tri,pTri,r0,bCells,d0,rBX,rBY,mask),...
                            'DerivativeCheck','off');
%      optimset = optimoptions('fmincon','Display','iter','Algorithm','sqp','TolFun',1e-5, ... 
%                         'MaxFunEvals',5e6,'MaxIter',3e3,'GradObj','on','GradConstr','on','Hessian','off',...
%                         'HessFcn',@(x,lambda) hessian(x,lambda,q,p,tri,pTri,r0,bCells,d0,rBX,rBY,mask),...
%                         'DerivativeCheck','off');     
                
    profile on
    [theta,ERes] = fmincon(energyFunc,x0,[],[],Aeq,beq,[],[],nonlinFunc,optimset);
    profile viewer
    pause
end

function [ H ] = hessian(theta,lambda,q,p,tri,pTri,rV,bCells,d0,RBx,RBy,mask)

    %% Compute hessian of objective.
    NB = size(d0,1);
    NC = size(d0,2);
    
    % Compute terms.
    dP = d0*p;
    dT = d0*theta;
    dQ = d0*q;
    QL = sqrt(sum(dQ.^2,2));
    
    rho = bsxfun(@rdivide,d0*bsxfun(@times,p,q),dP);
    Rsq = (p(bCells(:,1)).*p(bCells(:,2)).*QL.^2 - (dP .* dT))  ./ dP.^2;
    indZ = Rsq < 0;
    R = sqrt( Rsq );
    R(indZ) = 0;
    indZ = ~indZ;
    
    deltaX = -bsxfun(@minus,RBx,rho(:,1));
    deltaY = -bsxfun(@minus,RBy,rho(:,2));
    dMag = sqrt(deltaX.^2 + deltaY.^2);

    d = bsxfun(@minus,dMag,R);

    % Compute gradients.
    [ dR ] = fitDual.AFN.returnBondGradsTheta( q, theta, p, bCells );
    dR = bsxfun(@times,indZ,dR);
    
    %% Build objective hessian.
    NB = size(d0,1);
    nPix = sum(mask,2);
    Hobj =  bsxfun(@times,reshape(dR',[2,1,NB]),reshape(bsxfun(@times,nPix,dR)',[1,2,NB]));
    
    % Return bond hessians. 2a + 2d
    dAvg = indZ .* sum( d .* mask, 2 )/NB;
    [bondHess] = fitDual.AFN.returnBondHessianTheta( q, theta, p, bCells, dAvg );
    Hobj = real(Hobj/NB + bondHess + bondHess');
    
    %% Compute Hessian of each constraint.
    NV = size(pTri,1);
    NC = size(q,1);
    row = zeros(9*NV,1);
    col = zeros(9*NV,1);
    value = zeros(9*NV,1);
    
    [ r, indB, r0, lambdaV ] = fitDual.AFN.returnVertexPositions( q, theta, p, tri, pTri, rV );
    NV = size(r,1);
    
    [ dR ] = fitDual.AFN.returnThetaGrad( r, q, p, pTri, indB, r0, lambdaV );
    [ dR0 ] = fitDual.AFN.returnThetaGrad( r, q, p, pTri );
    [ r0, ~, etaStar ] = fitDual.AFN.returnR0( q, p, theta, pTri );
    r0Mag = sqrt(dot(r0,r0,2));
    etaMag = sqrt(sum(etaStar.^2,2));
    
    n = 1;
    for ii = 1:3
        for jj = 1:3
            II = ((n-1)*NV) + (1:NV);
            row(II) = pTri(:,ii);
            col(II) = pTri(:,jj);
            IPi = bsxfun(@times,r0(:,1),dR(:,1,ii)) + bsxfun(@times,r0(:,2),dR(:,2,ii));
            IPj = bsxfun(@times,r0(:,1),dR(:,1,jj)) + bsxfun(@times,r0(:,2),dR(:,2,jj));
            value(II) = lambda.ineqnonlin.*((dR0(:,1,ii) .* dR0(:,1,jj)) + (dR0(:,2,ii) .* dR0(:,2,jj)) - (IPi.*IPj ./ r0Mag)) .* (etaMag./r0Mag);
            n = n + 1;
        end
    end
    Hcon = sparse(row,col,value,NC,NC);

    H = Hobj + Hcon;
end

function [ H ] = oldHessian(x,lambda,q,p,tri,pTri,rV)

    % Compute vertex positions and gradient of vertex positions.
    [ r, indB, r0, lambdaV ] = fitDual.AFN.returnVertexPositions( q, x, p, tri, pTri, rV );
    delta = r - rV;
    NV = size(r,1);
    
    [ dR ] = fitDual.AFN.returnThetaGrad( r, q, p, pTri, indB, r0, lambdaV );
    [ dR0 ] = fitDual.AFN.returnThetaGrad( r, q, p, pTri );
    
    [ r0, S0, etaStar ] = fitDual.AFN.returnR0( q, p, x, pTri );
    r0Mag = sqrt(dot(r0,r0,2));
    etaMag = sqrt(sum(etaStar.^2,2));
    
    %% Compute hessian of objective.
    
    % Embed gradient into larger matrix.
    row = zeros(3*NV,1);
    col = zeros(3*NV,1);
    valueX = zeros(3*NV,1);
    valueY = zeros(3*NV,1);

    n = 1;
    for ii = 1:3
        II = ((n-1)*NV) + (1:NV);
        row(II) = 1:NV;
        col(II) = pTri(:,ii);
        valueX(II) = dR(:,1,ii);
        valueY(II) = dR(:,2,ii);
        n = n + 1;
    end
    
    dRx = sparse(row,col,valueX,NV,size(q,1));
    dRy = sparse(row,col,valueY,NV,size(q,1));

    % First component of Hessian is the outerproduct of these derivatives.
    Hobj = (dRx'*dRx) + (dRy'*dRy);
    
    % Second component of Hessian is computed from the second derivative of
    % vertex position w.r.t. theta.
    NC = size(q,1);
%     ind = (1:NV)';
    row = zeros(9*NV,1);
    col = zeros(9*NV,1);
    value = zeros(9*NV,1);
    
    denom = 2*S0 - (lambdaV .* etaMag.^2) - dot(r0,etaStar,2);
    n = 1;
    for ii = 1:3
        for jj = 1:3
            II = ((n-1)*NV) + (1:NV);
%             row(II) = ind;
%             col(II) = pTri(:,ii) + NC*(pTri(:,jj)-1);
            row(II) = pTri(:,ii);
            col(II) = pTri(:,jj);
            normalVal = ((dR(:,1,ii) .* dR(:,1,jj)) + (dR(:,2,ii) .* dR(:,2,jj))) .* dot(etaStar,delta,2) ./ denom;
%             normalVal = ((dR(:,1,ii) .* dR(:,1,jj)) + (dR(:,2,ii) .* dR(:,2,jj))) .* dot(etaStar,delta,2) ./ etaMag;  
            normalVal(indB) = 0;
            value(II) = normalVal;
            n = n + 1;
        end
    end
%     H2dv = sparse(row,col,value,NV,NC^2);
    H2dv = sparse(row,col,value,NC,NC);
%     Hobj = Hobj + reshape(sum(H2dv,1)',NC,NC);
    Hobj = Hobj + H2dv;
    Hobj = Hobj / size(r,1);
%     
    %% Compute Hessian of each constraint.
    row = zeros(9*NV,1);
    col = zeros(9*NV,1);
    value = zeros(9*NV,1);

    n = 1;
    for ii = 1:3
        for jj = 1:3
            II = ((n-1)*NV) + (1:NV);
%             row(II) = ind;
%             col(II) = pTri(:,ii) + NC*(pTri(:,jj)-1);
            row(II) = pTri(:,ii);
            col(II) = pTri(:,jj);
            IPi = bsxfun(@times,r0(:,1),dR(:,1,ii)) + bsxfun(@times,r0(:,2),dR(:,2,ii));
            IPj = bsxfun(@times,r0(:,1),dR(:,1,jj)) + bsxfun(@times,r0(:,2),dR(:,2,jj));
            value(II) = lambda.ineqnonlin.*((dR0(:,1,ii) .* dR0(:,1,jj)) + (dR0(:,2,ii) .* dR0(:,2,jj)) - (IPi.*IPj ./ r0Mag)) .* (etaMag./r0Mag);
            n = n + 1;
        end
    end
%     Hcon = sparse(row,col,value,NV,NC^2);
%     Hcon = sum(Hcon,1)';
    Hcon = sparse(row,col,value,NC,NC);
    % Add hessians together.
%     H = Hobj + reshape(Hcon,NC,NC);
    H = Hobj + Hcon;
end




