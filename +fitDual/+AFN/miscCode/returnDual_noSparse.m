function [ q, theta, p, tri, iCells, bulkV, ERes ] = returnDual( Struct, extCell, Q )

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
%     mask = bsxfun(@rdivide,mask,sum(mask,2));
    
    [ q, p, theta, tri, ~, ~, bulkV ] = fitDual.AFN.obtainStartPoint( Struct, 1 );
    q0 = [q,theta,p];

    Aeq = [[zeros(1,3*size(q,1)),ones(1,size(q,1))/size(q,1)]; ...
           [zeros(1,2*size(q,1)),ones(1,size(q,1))/size(q,1),zeros(1,size(q,1))]]; %...
%            [ones(1,size(q,1))/size(q,1),zeros(1,3*size(q,1))] ; ...
%            [zeros(1,size(q,1)),ones(1,size(q,1))/size(q,1),zeros(1,2*size(q,1))]];
       
    beq = [mean(p);mean(q0(:,3))];%mean(q0(:,1));mean(q0(:,2))];
    d0 = sparse(d0);
    
    optimset = optimoptions('fmincon','Display','iter','MaxFunEvals',1e6,'MaxIter',2e3,'TolFun',1e-6,...
               'Algorithm','interior-point','GradObj','on','GradConstr','on','DerivativeCheck','off',...
               'Hessian','on','HessFcn',@(x,lambda) hessian(x,lambda,bCells,d0,rBX,rBY,mask),'MaxSQPIter',40*size(d0,2));
    
    profile on
    [qF,ERes] = fmincon(@(q) fitDual.AFN.energy(q, d0, bCells, rBX, rBY, mask), ...
                q0,[],[],Aeq,beq,[],[],@(q) ensurePositive(q,bCells,d0),optimset);
    profile viewer
    q = qF(:,1:2);
    theta = qF(:,3);
    p = qF(:,4);
    
end

function [ H ] = hessian(x,lambda,bCells,d0,RBx,RBy,mask)

    Q = reshape(x,length(x)/4,4);
    NB = size(d0,1);
    NC = size(d0,2);
    
    q = Q(:,1:2);
    theta = Q(:,3);
    p = Q(:,4);
    
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
    deltaX = deltaX ./ dMag;
    deltaY = deltaY ./ dMag;

    d = bsxfun(@minus,dMag,R);

    % Compute gradients.
    [ dRhoX, dRhoY, dR ] = fitDual.AFN.returnBondGrads( q, theta, p, bCells );
    dR = bsxfun(@times,indZ,dR);
    
    %% Build objective hessian.
    % 1a + 2c
    rRatio = bsxfun(@rdivide,R,dMag);
    dNormXX = sum( rRatio .*( (deltaX .* deltaX) ) .* mask, 2);
    dNormXY = sum( rRatio .*( (deltaX .* deltaY) ) .* mask, 2);
    dNormYX = sum( rRatio .*( (deltaY .* deltaX) ) .* mask, 2);
    dNormYY = sum( rRatio .*( (deltaY .* deltaY) ) .* mask, 2);

    dRhoX = reshape(dRhoX',[1,8,NB]);
    dRhoY = reshape(dRhoY',[1,8,NB]);
    dR = reshape(dR',[1,8,NB]);

    dRhoXT = permute(dRhoX,[2,1,3]);
    dRhoYT = permute(dRhoY,[2,1,3]);
    dRT = permute(dR,[2,1,3]);
    
    Hobj = bsxfun(@times,dRhoXT,reshape(bsxfun(@times,dNormXX,dRhoX)',[1,8,NB])) + ...
           bsxfun(@times,dRhoXT,reshape(bsxfun(@times,dNormXY,dRhoY)',[1,8,NB])) + ...
           bsxfun(@times,dRhoYT,reshape(bsxfun(@times,dNormYX,dRhoX)',[1,8,NB])) + ...
           bsxfun(@times,dRhoYT,reshape(bsxfun(@times,dNormYY,dRhoY)',[1,8,NB]));

    % 1b    
    dNormX = indZ .* sum( deltaX .* mask, 2);
    dNormY = indZ .* sum( deltaY .* mask, 2);
     
    Hobj = Hobj - bsxfun(@times,dRT,reshape(bsxfun(@times,dNormX,dRhoX)',[1,8,NB])) - ...
           bsxfun(@times,dRT,reshape(bsxfun(@times,dNormY,dRhoY)',[1,8,NB])) - ... 
           bsxfun(@times,dRhoXT,reshape(bsxfun(@times,dNormX,dR)',[1,8,NB])) - ...
           bsxfun(@times,dRhoYT,reshape(bsxfun(@times,dNormY,dR)',[1,8,NB]));
       
    % 1c
    nPix = indZ .* sum(mask,2);
    Hobj = Hobj + bsxfun(@times,dRT,reshape(bsxfun(@times,nPix,dR)',[1,8,NB]));
       
    % 2b
    avgRatio = sum( (d./dMag) .* mask, 2);
    Hobj = Hobj + bsxfun(@times,dRhoXT,reshape(bsxfun(@times,avgRatio,dRhoX)',[1,8,NB])) + ...
           bsxfun(@times,dRhoYT,reshape(bsxfun(@times,avgRatio,dRhoY)',[1,8,NB]));
    Hobj = permute(Hobj,[3,1,2])/NB;
    
    % Return bond hessians. 2a + 2d + cons
    dNormX = sum( deltaX .* d .* mask, 2 ) / NB;
    dNormY = sum( deltaY .* d .* mask, 2 ) / NB;
    dAvg = indZ .* sum( d .* mask, 2 ) / NB;
    
    [ bondHess ] = fitDual.AFN.returnBondHessian( q, theta, p, bCells, dNormX, dNormY, dAvg, lambda.ineqnonlin ); 
    H = Hobj + bondHess;
    H = fitDual.AFN.deriv.loadHessian(H,bCells,NC,8);
    H = H + H';
    
end

function [ c, ceq, dc, dceq ] = ensurePositive(Q,bCells,d0) 
    
    q = Q(:,1:2);
    theta = Q(:,3);
    p = Q(:,4);

    ceq = [];
    c = ((d0*p) .* (d0*theta)) - p(bCells(:,1)).*p(bCells(:,2)).*sum((d0*q).^2,2);
    
    if (nargout > 2)
        dP = d0*p;
        dT = d0*theta;
        dQ = d0*q;
        QL = sum(dQ.^2,2);
        
        % X gradient
        gX = -2*(p(bCells(:,1)).*p(bCells(:,2)).*dQ(:,1));
        gX = bsxfun(@times,gX,d0);
        
        % Y gradient
        gY = -2*(p(bCells(:,1)).*p(bCells(:,2)).*dQ(:,2));
        gY = bsxfun(@times,gY,d0);

        % T gradient
        gTh = bsxfun(@times,dP,d0);
        
        % P gradient
        gP = bsxfun(@times,dT,d0);
        gP = gP - bsxfun(@rdivide,bsxfun(@times,(QL.*p(bCells(:,1)).*p(bCells(:,2))),abs(d0))',p)';
        
        dc = [gX,gY,gTh,gP]';
        dceq = [];
    end
end
