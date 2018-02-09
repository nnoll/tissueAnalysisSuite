function [ q, theta, p, tri, i0, ERes ] = returnDual( Struct, Text, Q )
    % RETURN DUAL 

    % Store bulk vertices.
    [ tri, bulk0, ext0, ~, ~, r0 ] = fitDual.returnGraph( Struct );
    optimset = optimoptions('fmincon','Display','none','Algorithm','active-set','MaxFunEvals',1e6,'MaxIter',2e3);

    if (nargin <= 2)
        if (nargin == 1)
            [ q, p, theta, tri, bulk0, ext0, ERes ] = fitDual.AFN.obtainStartPoint( Struct );
        else
            [ q, p, theta, tri, bulk0, ext0, ERes ] = fitDual.AFN.obtainStartPoint( Struct, Text );
        end
        i0 = [bulk0,ext0];
        x0(:,1:2) = q;
        x0(:,3) = theta;
        x0(:,4) = p;
        x0 = x0(:);
    else
        tri = fitDual.orderTri(Q(:,1:2),tri);
        i0 = [bulk0,ext0];
        x0 = Q(:);
    end
    
    [ d0, ~, ~, ~, ~ ] = fitDual.AFN.returnBonds( Struct, i0 );

    lb = -inf*ones(size(x0));
    lb((3/4*length(lb)+1):end) = .001;
    ub = inf*ones(size(x0));
    ub((3/4*length(lb)+1):end) = 1000;
    
    C = length(x0)/4;
    Aeq = [zeros(1,C),zeros(1,C),zeros(1,C),ones(1,C)];
    ExtCon = zeros(length(ext0),C);
    ExtCon( (1:length(ext0)) + length(ext0)*(length(bulk0):(length(i0)-1))) = 1;
    Aeq = [Aeq;[zeros(length(ext0),C),zeros(length(ext0),C),zeros(length(ext0),C),ExtCon]];
    beq = [C;ones(length(ext0),1)];
    
    d0_ext = full(d0(:,((length(bulk0)+1):size(d0,2))));
    badBonds = sum(abs(d0_ext),2) < 2;
    d0_ext(badBonds,:) = [];
    d0_ext = sparse(d0_ext);
    
    if (nargin == 1)
        Text = sum((d0_ext*x0((length(bulk0)+1):size(d0,2),1:2)).^2,2);
    end
    
    energyFunc = @(x) fitDual.AFN.energy(x,r0,tri);
    nonlinFunc = @(x) mycon( x, d0_ext, length(bulk0), Text );

    [x,ERes] = fmincon(energyFunc,x0,[],[],Aeq,beq,lb,ub,nonlinFunc,optimset);
       
    x = reshape(x,length(x)/4,4);
    q = x(:,1:2);
    theta = x(:,3);
    p = x(:,4);
        
end


function [ clin, ceq ] = mycon( x, d0, nBulk0, Tsq )

    x = reshape(x,length(x)/4,4);
    qE = x((nBulk0+1:size(x,1)),1:2);
    p = x(:,3);
    clin = std(p) - .5;
    ceq = sum( (d0*qE).^2, 2) - Tsq;
%     ceq = 0;
    
end

