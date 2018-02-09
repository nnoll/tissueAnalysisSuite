function [ q, p, tri, bulk0, ext0, bulkV, ERes ] = returnReducedDual( Struct, extCell, q, p )
    % RETURN REDUCED DUAL
    
    [ tri, bulk0, ext0, bulkV, x0, ~ ] = fitDual.returnGraph( Struct, extCell );

    [ ~, d0, t1, t2, ~, ~, ~, bCells, r1, r2, S ] = MI.returnBonds( Struct );

    t1 = bsxfun(@times,t1,S);
    t2 = bsxfun(@times,t2,S);
    
    if (nargin <= 2)
%         [ q, ~, tri, ~, ERes ] = fitDual.ATN.returnDual( Struct, extCell );
%         if (ERes > 1e2)
            q = x0(:,1:2);
%         end
        x0(:,1:2) = q;
        x0(:,3) = fitDual.AFN.seedPressure( q, bCells, t1, t2, r1, r2 );
        p = x0(:,3);
    else
        x0(:,1:2) = q;
        x0(:,3) = p;
    end

    B = d0*bsxfun(@times,q,p);
    deltaP = d0*p;

    B1 = B - bsxfun(@times,deltaP,r1);
    B2 = B - bsxfun(@times,deltaP,r2);

    L1 = sqrt(sum(B1.^2,2));    
    L2 = sqrt(sum(B2.^2,2));
    
    Lscale = .5*(mean(L1) + mean(L2));
    
    energyFunc = @(x) fitDual.AFN.reducedEnergy( x, sparse(d0), t1, t2, r1, r2, bCells );
    nonlinFunc = @(x) mycon( x, sparse(d0), Lscale, r1, r2, bCells );
    
    lb = -inf*ones(size(x0));
    lb(:,3) = .001;
    ub = inf*ones(size(x0));
    ub(:,3) = 1000;

    Aeq = [zeros(1,size(d0,2)),zeros(1,size(d0,2)),ones(1,size(d0,2))];
    beq = size(d0,2)*mean(p);
    
    optimset = optimoptions('fmincon','Display','none','Algorithm','interior-point','MaxFunEvals',2e6,...
                            'MaxIter',2e3,'TolFun',1e-5,'GradObj','on','GradConstr','on','Hessian','lbfgs',...
                            'HessFcn',@(x,lambda) hessian(x,lambda,sparse(d0),r1,r2,t1,t2,bCells),...
                            'DerivativeCheck','off'); 
                  
    [x,ERes] = fmincon(energyFunc,x0,[],[],Aeq,beq,lb,ub,nonlinFunc,optimset);

    q = x(:,1:2);
    p = x(:,3);

end

function [ H ] = hessian(x,lambda,d0,r1,r2,t1,t2,bCells)

    NB = size(d0,1);
    NC = size(d0,2);
    x = reshape(x,length(x)/3,3);
    q = x(:,1:2);
    p = x(:,3);
    
    B = d0*bsxfun(@times,q,p);
    deltaP = d0*p;

    B1 = B - bsxfun(@times,deltaP,r1);
    B2 = B - bsxfun(@times,deltaP,r2);

    L1 = sqrt(sum(B1.^2,2));
    B1 = bsxfun(@rdivide,B1,L1);
    
    L2 = sqrt(sum(B2.^2,2));
    B2 = bsxfun(@rdivide,B2,L2);

    [ drX1, drY1, drX2, drY2 ] = fitDual.AFN.reducedLocalGrads( q, p, bCells, r1 ,r2 );
    
    % Compute gradients and take an outer product
    IP1 = dot(B1,t1,2);
    IP2 = dot(B2,t2,2);

    gradE1 = bsxfun(@times,t1(:,1),drX1) + bsxfun(@times,t1(:,2),drY1);
    gradE2 = bsxfun(@times,t2(:,1),drX2) + bsxfun(@times,t2(:,2),drY2);

    Hobj = bsxfun(@times,reshape(gradE1',[6,1,NB]),reshape(gradE1',[1,6,NB])) + ...
           bsxfun(@times,reshape(gradE2',[6,1,NB]),reshape(gradE2',[1,6,NB]));
    
    Hobj = permute(Hobj,[3,1,2]);
    
    %% Add in components related to second derivatives of r
    T1X = IP1.*t1(:,1);
    T1Y = IP1.*t1(:,2);
    T2X = IP2.*t2(:,1);
    T2Y = IP2.*t2(:,2);
    
    [ H2nd ] = fitDual.AFN.returnReducedHessian( q, p, bCells, r1, r2, T1X, T1Y, T2X, T2Y );
    Hobj = (Hobj + H2nd)/NB;
    
    % Compute hessian of nonlinear constraint.
    [ dRX1, dRY1, dRX2, dRY2 ] = fitDual.AFN.reducedLocalGrads( q, p, bCells, r1 ,r2, 0 );
    
    Hcon = bsxfun(@times,reshape(drX1',[6,1,NB]),reshape(dRX1',[1,6,NB])) + ...
           bsxfun(@times,reshape(drY1',[6,1,NB]),reshape(dRY1',[1,6,NB])) + ...
           bsxfun(@times,reshape(drX2',[6,1,NB]),reshape(dRX2',[1,6,NB])) + ...
           bsxfun(@times,reshape(drY2',[6,1,NB]),reshape(dRY2',[1,6,NB]));
    Hcon = lambda.eqnonlin * permute(Hcon,[3,1,2]);
    
    [ HCon2nd ] = fitDual.AFN.returnReducedHessian( q, p, bCells, r1, r2 , B1(:,1), B1(:,2), B2(:,1), B2(:,2), 0 );
    Hcon = (Hcon + lambda.eqnonlin*HCon2nd) / (2*NB);

    H = Hobj + Hcon;
    [ H ] = fitDual.AFN.deriv.loadHessian( H, bCells, NC, 6 );
    H = H + H';
    
end

function [ clin, ceq, dlin, deq ] = mycon( x, d0, Lscale, r1, r2, bCells )

    NB = size(d0,1);
    NC = size(d0,2);
    q = x(:,1:2);
    p = x(:,3);
    clin = [];
    
%     dQ = d0*q;
%     ceq = mean(sqrt(sum(dQ.^2,2))) - Lscale;
    B = d0*bsxfun(@times,q,p);
    deltaP = d0*p;

    B1 = B - bsxfun(@times,deltaP,r1);
    B2 = B - bsxfun(@times,deltaP,r2);

    L1 = sqrt(sum(B1.^2,2));
    B1 = bsxfun(@rdivide,B1,L1);
    
    L2 = sqrt(sum(B2.^2,2));
    B2 = bsxfun(@rdivide,B2,L2);
    
    ceq = .5*(mean(L1) + mean(L2)) - Lscale;
    if (nargout > 2)
        dlin = [];
%         deq  = [d0'*bsxfun(@rdivide,dQ,sqrt(sum(dQ.^2,2))),zeros(size(p))];
%         deq = deq(:) / size(d0,1);
        [ drX1, drY1, drX2, drY2 ] = fitDual.AFN.reducedLocalGrads( q, p, bCells, r1 ,r2, 0 );
%         II = 1:size(d0,1);
%         deq = .5*mean( (sparse(II,II,B1(:,1))*drX1) + (sparse(II,II,B1(:,2))*drY1),1)' + ...
%               .5*mean( (sparse(II,II,B2(:,1))*drX2) + (sparse(II,II,B2(:,2))*drY2),1)';
        deq = .5*(bsxfun(@times,B1(:,1),drX1) + bsxfun(@times,B1(:,2),drY1)) + ...
              .5*(bsxfun(@times,B2(:,1),drX2) + bsxfun(@times,B2(:,2),drY2));
        rows = [bCells(:,1);bCells(:,1)+NC;bCells(:,1)+2*NC;bCells(:,2);bCells(:,2)+NC;bCells(:,2)+2*NC];
        deq = accumarray(rows,deq(:)/NB,[3*NC,1]);
    end
    
end

