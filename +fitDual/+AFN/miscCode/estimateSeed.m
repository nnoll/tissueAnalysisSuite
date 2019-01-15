function [ q, p ] = estimateSeed( Struct, extCell )
    % ESTIMATE SEED 

    [ ~, bulk0, ext0, ~, ~ ] = fitDual.returnGraph( Struct, extCell );
    i0 = [bulk0,ext0];
%     [ x0 ] = fitDual.ATN.fitTensionGraph( Struct, extCell );
%     x0 = [x0,ones(size(x0,1),1)];
    
    [ d0, t1, t2, r1, r2 ] = fitDual.AFN.returnBonds( Struct, i0 );
    
    L1_q = [bsxfun(@times,t1(:,1),d0),bsxfun(@times,t1(:,2),d0)];
    L1_p = bsxfun(@times,dot(t1,r1,2),d0);
    L1 = [L1_q,L1_p];
    
    L2_q = [bsxfun(@times,t2(:,1),d0),bsxfun(@times,t2(:,2),d0)];
    L2_p = bsxfun(@times,dot(t2,r2,2),d0);
    L2 = [L2_q,L2_p];
    
    L = full([L1;L2]);
%     Aeq = [[ones(1,size(d0,2)),zeros(1,size(d0,2)),zeros(1,size(d0,2))]; ...
%           [zeros(1,size(d0,2)),ones(1,size(d0,2)),zeros(1,size(d0,2))]; ...
%           [zeros(1,size(d0,2)),zeros(1,size(d0,2)),ones(1,size(d0,2))]];
%     beq = zeros(3,1);
%     beq(end) = size(d0,2);
    
    L = [L; [ones(1,size(d0,2)),zeros(1,size(d0,2)),zeros(1,size(d0,2))]; ...
            [zeros(1,size(d0,2)),ones(1,size(d0,2)),zeros(1,size(d0,2))]; ...
            [zeros(1,size(d0,2)),zeros(1,size(d0,2)),ones(1,size(d0,2))]];
    b = zeros(size(L,1),1);
    b(end) = size(d0,2);
    
    L = sparse(L);
    x = L \ b;
    plot(x)
    pause
%     optimset = optimoptions('fmincon','Display','iter','Algorithm','active-set','MaxFunEvals',2e6,'MaxIter',1e3,'GradObj','off','GradConstr','off');
    
%     F = @(x) mean((L*x).^2);
  
%     nonlinFunc = @(x) mycon( x, d0 );
    
%     R = mean(x0(:,1:2),1);
%     x0(:,1:2) = bsxfun(@minus,x0(:,1:2),R);
%     x0(:,1:2) = bsxfun(@rdivide,x0(:,1:2), mean(sqrt( sum((d0*x0(:,1:2)).^2, 2) )));
%     x0(:,3) = 1;
%     x0 = x0(:);
    
%     [x,ERes] = fmincon(F,x0,[],[],Aeq,beq,[],[],nonlinFunc,optimset);

    x = reshape(x,length(x)/3,3);
    q = x(:,1:2);
    p = x(:,3);
    q = bsxfun(@rdivide,q,p);
    
    % Rescale q so that we match onto cell centroid scale.
    [ d0, rho ] = fitDual.AFN.returnBondCentroid( Struct, i0 );
    
    rhoEst = bsxfun(@rdivide,d0*q,d0*p);
    lambda = mean(dot(rho,rho,2))/mean(dot(rho,rhoEst,2));
    
    q = lambda*q;
    p = lambda*p;
    p = p - mean(p) + 1;
    
    q = bsxfun(@plus,q,R);
    
end

function [ clin, ceq ] = mycon( x, d0 ) 
    q = reshape(x,length(x)/3,3);
    dQ = d0*q(:,1:2);
    ceq = mean(sqrt(sum(dQ.^2,2))) - 1;
    clin = 0;
end

