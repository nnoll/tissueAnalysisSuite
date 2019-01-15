function [ q, p, theta ] = seedDual( Struct )
    % SEED DUAL 

    [ T, p, dC, ~, ~, iCells, S1 ] = MI.invertMech( Struct, 1, 3 );
    [ ~, ~, t1, t2, ~, ~, rho, bCells, r1, r2, S2 ] = MI.returnBonds( Struct );
    
    S = S1.*S2;

%     for c = 1:length(iCells)
%         q2(c,:) = p(c) * Struct.Cdat(iCells(c)).q;
%     end
 
    dP = dC * p;
    
    r1 = bsxfun(@times,dP,r1);
    r2 = bsxfun(@times,dP,r2);
    
%     t1R = bsxfun(@times,S,t1*[0,-1;1,0]);
%     t2R = bsxfun(@times,S,t2*[0,1;-1,0]);
    
    t1 = bsxfun(@times,S,t1);
    t2 = bsxfun(@times,S,-t2);
    
%     b1R = dot(r1,t1,2) - T;
%     b2R = dot(r2,t2,2) - T;
    b1 = dot(r1,t1,2);
    b2 = dot(r2,t2,2);
    
%     L1R = [bsxfun(@times,t1R(:,1),dC),bsxfun(@times,t1R(:,2),dC)];
%     L2R = [bsxfun(@times,t2R(:,1),dC),bsxfun(@times,t2R(:,2),dC)];

    L1 = [bsxfun(@times,t1(:,1),dC),bsxfun(@times,t1(:,2),dC)];
    L2 = [bsxfun(@times,t2(:,1),dC),bsxfun(@times,t2(:,2),dC)];

%     L = sparse([L1R;L2R;L1;L2;[ones(1,size(dC,2)),zeros(1,size(dC,2))]; [zeros(1,size(dC,2)),ones(1,size(dC,2))] ]);
    L = sparse([L1;L2;[ones(1,size(dC,2)),zeros(1,size(dC,2))]; [zeros(1,size(dC,2)),ones(1,size(dC,2))] ]);
%     b = [b1R;b2R;b1;b2;0;0];
    b = [b1;b2;0;0];
%     plot(L1*q2(:)-b1)
%     pause
    
    q = L \ b;

    q = reshape(q,length(q)/2,2);
    q = bsxfun(@rdivide,q,p);
    
    dQ = dC * q;
    dQ = sum(dQ.^2,2);
    dQ = dQ .* (p(bCells(:,1)).*p(bCells(:,2)));
    
    Tresid = T.^2 - dQ;
    Ltheta = bsxfun(@times,dP,dC);
    Ltheta = [Ltheta;ones(1,size(Ltheta,2))/size(Ltheta,2)];
    Tresid = [Tresid;0];
    theta = sparse(Ltheta) \ Tresid;    
    
    theta = theta./p;
    
end

