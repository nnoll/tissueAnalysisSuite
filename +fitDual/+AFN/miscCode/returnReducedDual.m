function [ q, p, tri, bulk0, ext0, bulkV, ERes ] = returnReducedDual( Struct, extCell, q, p )
    % RETURN REDUCED DUAL
    
    [ tri, bulk0, ext0, bulkV, x0, rV ] = fitDual.returnGraph( Struct, extCell );

    [ ~, d0, t1, t2, ~, ~, rhoA, bCells, r1, r2, S ] = MI.returnBonds( Struct );

    t1 = bsxfun(@times,t1,S);
    t2 = bsxfun(@times,t2,S);
    
    if (nargin <= 2)
        [ q, ~, tri, ~, ERes ] = fitDual.ATN.returnDual( Struct, extCell );
        if (ERes > 1e3)
            q = x0(:,1:2);
        end
        Qcom = mean(q,1);
        Lscale = mean(sqrt(sum((d0*q(:,1:2)).^2,2)));
        x0(:,1:2) = q;
        x0(:,3) = fitDual.AFN.seedPressure( q, bCells, t1, t2, r1, r2 );
        p = x0(:,3);
    else
        Qcom = mean(q,1);
        Lscale = mean(sqrt(sum((d0*q(:,1:2)).^2,2)));
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
    lb(:,3) = .01;
    ub = inf*ones(size(x0));
    ub(:,3) = 100;
    
%     Aeq = [[ones(1,size(d0,2)),zeros(1,size(d0,2)),zeros(1,size(d0,2))]; ...
%           [zeros(1,size(d0,2)),ones(1,size(d0,2)),zeros(1,size(d0,2))]; ...
%           [zeros(1,size(d0,2)),zeros(1,size(d0,2)),ones(1,size(d0,2))]];
%     
%     beq = zeros(3,1);
%     
%     beq(1) = size(d0,2)*Qcom(1);
%     beq(2) = size(d0,2)*Qcom(2);
%     beq(3) = size(d0,2)*mean(p);

    Aeq = [zeros(1,size(d0,2)),zeros(1,size(d0,2)),ones(1,size(d0,2))];
    beq = size(d0,2)*mean(p);
    
    optimset = optimoptions('fmincon','Display','iter','Algorithm','interior-point','MaxFunEvals',2e6,...
                            'MaxIter',2e3,'TolFun',1e-5,'GradObj','on','GradConstr','on','Hessian','on',...
                            'HessFcn',@(x,lambda) hessian(x,lambda,sparse(d0),r1,r2,t1,t2,bCells),...
                            'DerivativeCheck','on');
%     optimset = optimoptions('fmincon','Display','iter','Algorithm','sqp','MaxFunEvals',2e6,...
%                             'MaxIter',2e3,'TolFun',1e-5,'GradObj','on','GradConstr','on','Hessian','on',...
%                             'HessFcn',@(x,lambda) hessian(x,lambda,sparse(d0),r1,r2,t1,t2,bCells),...
%                             'DerivativeCheck','off');
    profile on                   
    [x,ERes] = fmincon(energyFunc,x0,[],[],Aeq,beq,lb,ub,nonlinFunc,optimset);
    profile viewer
    q = x(:,1:2);
    p = x(:,3);

%     rho = bsxfun(@rdivide,d0*(bsxfun(@times,p,q)) , (d0*p) );
%     scatter(rho(:,1),rhoA(:,1))
%     pause
%     
%     scatter(rho(:,2),rhoA(:,2))
%     pause
end

function [ H ] = hessian(x,lambda,d0,r1,r2,t1,t2,bCells)

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

%     II = 1:length(L1);
% 
%     T1X = sparse(II,II,t1(:,1));
%     T1Y = sparse(II,II,t1(:,2));
%     T2X = sparse(II,II,t2(:,1));
%     T2Y = sparse(II,II,t2(:,2));

%     gradE1 = (T1X * drX1) + (T1Y * drY1);
%     gradE2 = (T2X * drX2) + (T2Y * drY2);
    gradE1 = bsxfun(@times,t1(:,1),drX1) + bsxfun(@times,t1(:,2),drY1);
    gradE2 = bsxfun(@times,t2(:,1),drX2) + bsxfun(@times,t2(:,2),drY2);

%     Hobj = (gradE1' * gradE1) + (gradE2' * gradE2);
    Hobj = bsxfun(@times,reshape(gradE1',[size(gradE1,2),1,size(gradE1,1)]),reshape(gradE1',[1,size(gradE1,2),size(gradE1,1)])) + ...
           bsxfun(@times,reshape(gradE2',[size(gradE2,2),1,size(gradE2,1)]),reshape(gradE2',[1,size(gradE2,2),size(gradE2,1)]));
%     size(Hobj)
%     pause
    
    %% Add in components related to second derivatives of r
    T1X = IP1.*t1(:,1);
    T1Y = IP1.*t1(:,2);
    T2X = IP2.*t2(:,1);
    T2Y = IP2.*t2(:,2);
    
    [ H2nd ] = fitDual.AFN.returnReducedHessian( q, p, bCells, r1, r2 , T1X, T1Y, T2X, T2Y );
    Hobj = Hobj + H2nd;
    
    % Compute hessian of nonlinear constraint.
    [ dRX1, dRY1, dRX2, dRY2 ] = fitDual.AFN.reducedLocalGrads( q, p, bCells, r1 ,r2, 0 );
    
    Hcon = bsxfun(@times,reshape(drX1',[6,1,size(drX1,1)]))*bsxfun(@times,reshape(dRX1',[1,6,size(dRX1,1)])) + ...
           bsxfun(@times,reshape(drY1',[6,1,size(drY1,1)]))*bsxfun(@times,reshape(dRY1',[1,6,size(dRY1,1)])) + ...
           bsxfun(@times,reshape(drX2',[6,1,size(drX2,1)]))*bsxfun(@times,reshape(dRX2',[1,6,size(dRX2,1)])) + ...
           bsxfun(@times,reshape(drY2',[6,1,size(drY2,1)]))*bsxfun(@times,reshape(dRY2',[1,6,size(dRY2,1)]));
    Hcon = lambda.eqnonlin * Hcon;
    
    [ HCon2nd ] = fitDual.AFN.returnReducedHessian( q, p, bCells, r1, r2 , B1(:,1), B1(:,2), B2(:,1), B2(:,2), 0 );
    HCon2nd = (HCon2nd + HCon2nd');
    Hcon = Hcon + HCon2nd) / (2*size(d0,1));

    H = Hobj/size(d0,1) + Hcon;

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

