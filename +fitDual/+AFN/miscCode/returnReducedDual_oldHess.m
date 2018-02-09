function [ q, p, tri, bulk0, ext0, ERes ] = returnReducedDual( Struct, extCell, q, p )
    % RETURN REDUCED DUAL
    
    [ tri, bulk0, ext0, ~, x0, rV ] = fitDual.returnGraph( Struct, extCell );

    [ ~, d0, t1, t2, ~, ~, ~, bCells, r1, r2, S ] = MI.returnBonds( Struct );

    t1 = bsxfun(@times,t1,S);
    t2 = bsxfun(@times,t2,S);
    
    if (nargin <= 2)
%         [ q, p, theta ] = MI.seedDual( Struct );
        [ q, ~, tri, ~, ~ ] = fitDual.ATN.returnDual( Struct, extCell );
%         q = x0(:,1:2);
        Qcom = mean(q,1);
        Lscale = mean(sqrt(sum((d0*q(:,1:2)).^2,2)));
        x0(:,1:2) = q;
        x0(:,3) = fitDual.AFN.seedPressure( q, bCells, t1, t2, r1, r2 );
    else
        Qcom = mean(q,1);
        Lscale = mean(sqrt(sum((d0*q(:,1:2)).^2,2)));
        x0(:,1:2) = q;
        x0(:,3) = p;
    end

%     d0_1 = (d0 == 1)';
%     d0_2 = -(d0 == -1)';

    energyFunc = @(x) fitDual.AFN.reducedEnergy( x, sparse(d0), t1, t2, r1, r2, bCells );
    nonlinFunc = @(x) mycon( x, sparse(d0), Lscale );
    
    lb = -inf*ones(size(x0));
    lb(:,3) = .01;
    ub = inf*ones(size(x0));
    ub(:,3) = 100;
    
    Aeq = [[ones(1,size(d0,2)),zeros(1,size(d0,2)),zeros(1,size(d0,2))]; ...
          [zeros(1,size(d0,2)),ones(1,size(d0,2)),zeros(1,size(d0,2))]; ...
          [zeros(1,size(d0,2)),zeros(1,size(d0,2)),ones(1,size(d0,2))]];
    beq = zeros(3,1);
    
    beq(1) = size(d0,2)*Qcom(1);
    beq(2) = size(d0,2)*Qcom(2);
    
    if (nargin == 2)
        beq(end) = size(d0,2);
    else
        beq(end) = size(d0,2)*mean(p);
    end

    optimset = optimoptions('fmincon','Display','iter','Algorithm','interior-point','MaxFunEvals',2e6,...
                            'MaxIter',1e3,'TolFun',1e-5,'GradObj','on','GradConstr','on','Hessian','on',...
                            'HessFcn',@(x,lambda) hessian(x,lambda,sparse(d0),r1,r2,t1,t2,bCells));
   
    [x,ERes] = fmincon(energyFunc,x0,[],[],Aeq,beq,lb,ub,nonlinFunc,optimset);
        
    q = x(:,1:2);
    p = x(:,3);

end

function [ H ] = hessian(x,lambda,d0,r1,r2,t1,t2,bCells)

    C = size(d0,2);
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

    [ drX1, drY1, drX2, drY2 ] = fitDual.AFN.reducedLocalGrads( q, p, bCells, d0, r1 ,r2 );
    
    % Compute gradients and take an outer product
    IP1 = dot(B1,t1,2);
    IP2 = dot(B2,t2,2);

    II = 1:length(L1);

    sp1 = (sparse(II,II,B1(:,1)) * drX1) + (sparse(II,II,B1(:,2)) * drY1);
    sp2 = (sparse(II,II,B2(:,1)) * drX2) + (sparse(II,II,B2(:,2)) * drY2);

    dBX1 = (sparse(II,II,1./L1) * drX1) - (sparse(II,II,B1(:,1)./L1) * sp1);
    dBY1 = (sparse(II,II,1./L1) * drY1) - (sparse(II,II,B1(:,2)./L1) * sp1);

    dBX2 = (sparse(II,II,1./L2) * drX2) - (sparse(II,II,B2(:,1)./L2) * sp2);
    dBY2 = (sparse(II,II,1./L2) * drY2) - (sparse(II,II,B2(:,2)./L2) * sp2);

    gradE1 = (sparse(II,II,t1(:,1)) * dBX1) + (sparse(II,II,t1(:,2)) * dBY1); 
    gradE2 = (sparse(II,II,t2(:,1)) * dBX2) + (sparse(II,II,t2(:,2)) * dBY2); 
        
    Hobj = (gradE1' * gradE1) + (gradE2' * gradE2);

    %% Add in components related to second derivatives of r
    
    % First term
    Oxx = sparse(II,II, IP1 .* (3*IP1.*B1(:,1).*B1(:,1) - t1(:,1).*B1(:,1) - t1(:,1).*B1(:,1) - IP1)./L1.^2);
    Hobj = Hobj + (drX1'*Oxx*drX1);
    
    Oxy = sparse(II,II, IP1 .* (3*IP1.*B1(:,1).*B1(:,2) - t1(:,1).*B1(:,2) - t1(:,2).*B1(:,1))./L1.^2);
    Hobj = Hobj + (drX1'*Oxy*drY1);
    
    Oyx = sparse(II,II, IP1 .* (3*IP1.*B1(:,2).*B1(:,1) - t1(:,2).*B1(:,1) - t1(:,1).*B1(:,2))./L1.^2);
    Hobj = Hobj + (drY1'*Oyx*drX1);
    
    Oyy = sparse(II,II, IP1 .* (3*IP1.*B1(:,2).*B1(:,2) - t1(:,2).*B1(:,2) - t1(:,2).*B1(:,2) - IP1)./L1.^2);
    Hobj = Hobj + (drY1'*Oyy*drY1);
    
    % Second term
    Oxx = sparse(II,II, IP2 .* (3*IP2.*B2(:,1).*B2(:,1) - t2(:,1).*B2(:,1) - t2(:,1).*B2(:,1) - IP2)./L2.^2);
    Hobj = Hobj + (drX2'*Oxx*drX2);
    
    Oxy = sparse(II,II, IP2 .* (3*IP2.*B2(:,1).*B2(:,2) - t2(:,1).*B2(:,2) - t2(:,2).*B2(:,1))./L2.^2);
    Hobj = Hobj + (drX2'*Oxy*drY2);
    
    Oyx = sparse(II,II, IP2 .* (3*IP2.*B2(:,2).*B2(:,1) - t2(:,2).*B2(:,1) - t2(:,1).*B2(:,2))./L2.^2);
    Hobj = Hobj + (drY2'*Oyx*drX2);
    
    Oyy = sparse(II,II, IP2 .* (3*IP2.*B2(:,2).*B2(:,2) - t2(:,2).*B2(:,2) - t2(:,2).*B2(:,2) - IP2)./L2.^2);
    Hobj = Hobj + (drY2'*Oyy*drY2);
    
    %% Take care of last second derivative term.
    H1x = (t1(:,1) - IP1.*B1(:,1))./L1;
    H1y = (t1(:,2) - IP1.*B1(:,2))./L1;
    H2x = (t2(:,1) - IP2.*B2(:,1))./L2;
    H2y = (t2(:,2) - IP2.*B2(:,2))./L2;

    H1x = d0' * (IP1.*H1x);
    H1y = d0' * (IP1.*H1y);
    H2y = d0' * (IP2.*H2y);
    H2x = d0' * (IP2.*H2x);

    H1 = sparse([1:C,(C+1):2*C,(2*C+1):3*C,(2*C+1):3*C], [(2*C+1):3*C,(2*C+1):3*C,1:C,(C+1):2*C], [H1x',H1y',H1x',H1y'],3*C,3*C);
    H2 = sparse([1:C,(C+1):2*C,(2*C+1):3*C,(2*C+1):3*C], [(2*C+1):3*C,(2*C+1):3*C,1:C,(C+1):2*C], [H2x',H2y',H2x',H2y'],3*C,3*C);    
    Hobj = Hobj + H1 + H2;

    % Compute hessian of nonlinear constraint.
    dQ = d0*q;
    QBL = sqrt(sum(dQ.^2,2));
    II = 1:size(dQ,1);
    Hcon = sparse([ [d0' * sparse(II,II,lambda.eqnonlin*(1-dQ(:,1).*dQ(:,1))./QBL) * d0, ...
              d0' * sparse(II,II,-lambda.eqnonlin*dQ(:,1).*dQ(:,2)./QBL) * d0,sparse(C,C)];...
             [d0' * sparse(II,II,-lambda.eqnonlin*dQ(:,2).*dQ(:,1)./QBL) * d0,...
              d0' * sparse(II,II,lambda.eqnonlin*(1-dQ(:,2).*dQ(:,2))./QBL) * d0,sparse(C,C)];...
              sparse(C,3*C)]);
    
    H = Hobj/size(d0,1) + Hcon/size(d0,1);

end

function [ clin, ceq, dlin, deq ] = mycon( x, d0, Lscale )

    q = x(:,1:2);
    p = x(:,3);
%     clin = std(p) - 5;
    clin = [];
    
    dQ = d0*q;
    ceq = mean(sqrt(sum(dQ.^2,2))) - Lscale;

    if (nargout > 2)
%         dlin = [zeros(size(q)),2/length(p) * (p - sum(p))];
        dlin = [];
        deq  = [d0'*bsxfun(@rdivide,dQ,sqrt(sum(dQ.^2,2))),zeros(size(p))];
        dlin = dlin(:);
        deq = deq(:) / size(d0,1);
    end
    
end

