function [ q, p, theta ] = seedDual( Struct )
    % SEED DUAL 

    [ T, pI, dC, ~, ~, iCells, S1 ] = MI.invertMech( Struct, 1, 3 );
    [ ~, ~, t1, t2, ~, ~, rho, bCells, r1, r2, S2 ] = MI.returnBonds( Struct );
    
    S = S1.*S2;
    
    t1R = bsxfun(@times,S,t1*[0,-1;1,0]);
    t2R = bsxfun(@times,S,t2*[0,1;-1,0]);
    
    t1 = bsxfun(@times,S,t1);
    t2 = bsxfun(@times,S,-t2);
    
    b1R = -T;
    b2R = -T;
    b1 = zeros(size(T));
    b2 = zeros(size(T));
    
    L1R = [bsxfun(@times,t1R(:,1),dC),bsxfun(@times,t1R(:,2),dC),bsxfun(@times,dot(r1,t1R,2),dC)];
    L2R = [bsxfun(@times,t2R(:,1),dC),bsxfun(@times,t2R(:,2),dC),bsxfun(@times,dot(r2,t2R,2),dC)];

    L1 = [bsxfun(@times,t1(:,1),dC),bsxfun(@times,t1(:,2),dC),bsxfun(@times,dot(r1,t1,2),dC)];
    L2 = [bsxfun(@times,t2(:,1),dC),bsxfun(@times,t2(:,2),dC),bsxfun(@times,dot(r2,t2,2),dC)];

    L = sparse([L1R;L2R;L1;L2;[ones(1,size(dC,2))/size(dC,2),zeros(1,size(dC,2)),zeros(1,size(dC,2))]; ...
                              [zeros(1,size(dC,2)),ones(1,size(dC,2))/size(dC,2),zeros(1,size(dC,2))];
                              [zeros(1,size(dC,2)),zeros(1,size(dC,2)),ones(1,size(dC,2))/size(dC,2)] ]);
    b = [b1R;b2R;b1;b2;0;0;1];
    
    q = L \ b;

    q = reshape(q,length(q)/3,3);
    p = q(:,3);

    q = q(:,1:2);
    q = bsxfun(@rdivide,q,p);
    
    dQ = dC * q;
    dQ = sum(dQ.^2,2);
    dQ = dQ .* (p(bCells(:,1)).*p(bCells(:,2)));
    dP = dC * p;
    
    Tresid = T.^2 - dQ;
    Ltheta = bsxfun(@times,dP,dC);
    Ltheta = [Ltheta;ones(1,size(Ltheta,2))/size(Ltheta,2)];
    Tresid = [Tresid;0];
    theta = sparse(Ltheta) \ Tresid;    
    
    theta = theta./p;
    
end

