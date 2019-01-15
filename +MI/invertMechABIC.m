function [ T, P, Imax ] = invertMechABIC( Struct, extCell, mu )
    % INVERT MECH 

    Rot = [0,-1;1,0];
    [ dC, dV, ~, iVerts, iCells ] = fitDual.ATN.computeDiffOperators( Struct, extCell );  
 
    %Log vertex positions.
    rv = zeros(length(iVerts),2);
    for ii = 1:length(iVerts)
        rv(ii,1) = double(Struct.Vdat(iVerts(ii)).vertxcoord);
        rv(ii,2) = double(Struct.Vdat(iVerts(ii)).vertycoord);
    end
    
    %Obtain bond vectors.
    rb = dV * rv;
    D = sqrt(sum(rb.^2,2));
    nb = rb*Rot;
    rb = bsxfun(@rdivide, rb, D);
    
    adj = dV'*dV;
    Z = diag(adj);
    bulkVerts = Z > 1;
    
    % Build MI matrix
    mTension = [bsxfun(@times,dV(:,bulkVerts),rb(:,1))';bsxfun(@times,dV(:,bulkVerts),rb(:,2))'];
        
    rC =  zeros(length(iCells),2);
    for ii = 1:length(iCells)
        rC(ii,:) = Struct.Cdat(iCells(ii)).centroid.coord;
    end
        
    nC = dC * rC;
    S = sign(dot(nC,nb,2));

    dC = bsxfun(@times,S,dC);
    mPressure = [bsxfun(@times,abs(dV(:,bulkVerts)),nb(:,1))'*dC;bsxfun(@times,abs(dV(:,bulkVerts)),nb(:,2))'*dC];
    
    A = [mTension,-mPressure];
    n = size(A,1); 
    m = size(A,2);
    B = zeros(m,m);
    B( (1:size(mTension,2)) + m*(0:(size(mTension,2)-1))) = 1;

    b = zeros(size(A,1),1);
    g = ones(m,1);
    
    detBtB = eig(sparse(B'*B));
    m0 = sum(detBtB == 0);
    
    if (nargin == 2)
        infoFunc = @(x) MI.maxABIC( x, A, b, B, g, m0, n, m );
        optimset = optimoptions('fmincon','Display','iter','Algorithm','active-set','MaxIter',500,'MaxFunEvals',500);
        x0 = 1;
        lb = 1e-8;
        ub = inf;
        [x,Imax] = fmincon(infoFunc,x0,[],[],[],[],lb,ub,[],optimset);
    else
        x = mu;
    end
    
    S = [ [A,b]; [sqrt(x)*B,sqrt(x)*g] ];
    [~,R] = qr(S);
    Rdiag = R(1:(m+1),(size(R,2)-m):size(R,2));
   
    h = Rdiag(1:m,m+1);
    H = Rdiag(1:m,1:m);
    Mech = pinv(H) * h;

    T = Mech(1:size(mTension,2));
    P = Mech((size(mTension,2)+1):end);
            
end

