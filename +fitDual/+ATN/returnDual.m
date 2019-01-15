function [ q, theta, tri, cells, bulkV, ERes ] = returnDual( Struct, extCell, Q )
    % RETURN DUAL 

    % Store bulk vertices.
    [ tri, bulk0, ext0, bulkV, xC, r0 ] = fitDual.returnGraph( Struct, extCell );
    [ ~, dV ] = fitDual.ATN.computeDiffOperators( Struct, 1 );  

    dV = dV(:,1:length(bulkV));
    badBonds = sum(abs(dV),2) == 1;
    dV(badBonds,:) = [];
    
%     rExt = [Struct.Vdat(extV).vertxcoord;Struct.Vdat(extV).vertycoord]';
    if (nargin <= 2)
        [ q, theta ] = fitDual.ATN.obtainStartPoint( Struct, extCell );
        x0 = zeros(size(xC,1),3);
%         q = xC(:,1:2);
        x0(:,1:2) = q;
        x0(:,3) = theta;
        x0 = x0(:);
    else
        tri = fitDual.orderTri(Q,tri);
        q = Q(:,1:2);
        theta = Q(:,3);
        x0 = Q(:);
    end
    
    [ tri ] = fitDual.orderTri( q, tri );
    
    [ extBonds, rBext ] = fitDual.ATN.findExtBonds( Struct, ext0 );
    [ A, b ] = fitDual.ATN.returnBoundaryConditions( x0, bulk0, ext0, extBonds, rBext );
    
    energyFunc = @(x) fitDual.ATN.energy(x,r0,tri,A,b);
    
    [ q2, theta2, tri ] = fitDual.ATN.fitThetaModes( xC(:,1:2), Struct, 1 );
    x1 = [q(:);zeros(size(theta))];
    x2 = [q2(:);theta2];
    x3 = [q2(:);zeros(size(theta2))];

    E = [energyFunc(x0),energyFunc(x1),energyFunc(x2),energyFunc(x3)];
    [~,ind] = min(E);
    
    if (ind == 2)
        x0 = x1;
    elseif (ind == 3)
        x0 = x2;
    elseif (ind == 4)
        x0 = x3;
    end

%     optimset = optimoptions('fmincon','Display','iter','Algorithm','sqp','MaxFunEvals',5e6,...
%                'MaxIter',2e3,'TolFun',1e-5,'GradObj','on','GradConstr','on','DerivativeCheck','off');
    optimset = optimoptions('fminunc','Display','none','Algorithm','Quasi-Newton','MaxFunEvals',5e6,...
               'MaxIter',1e3,'TolFun',1e-5,'GradObj','on');
%            
%     [x,ERes] = fmincon(energyFunc,x0,[],[],[],[],[],[],@(x) enforcePositivity(x,tri,dV), optimset);
    [x,ERes] = fminunc(energyFunc,x0,optimset);
    
    x = reshape(x,length(x)/3,3);
    q = x(:,1:2);
    theta = x(:,3);
    
    cells = [bulk0,ext0];
    
end


function [ c, ceq, dc, dceq ] = enforcePositivity(x,tri,d0) 

    x = reshape(x,length(x)/3,3);
    q = x(:,1:2);
    theta = x(:,3);
    [ tri ] = fitDual.orderTri(q,tri);
    
    if (nargout <= 2)
        [ rB, rTh ] = fitDual.ATN.returnPartitionedVertexPositions( q, theta, tri );
    else
        [ rB, rTh, dRb, dRTh ] = fitDual.ATN.returnPartitionedVertexPositions( q, theta, tri );
    end

    ceq = [];
    c = .5*(sum((d0*rTh).^2,2) -  sum((d0*rB).^2,2));
    
    if (nargout > 2)
        % Calculate bond vectors
        bTh = d0*rTh;
        bB = d0*rB;
        V = size(rTh,1);
        C = size(q,1);

%         rBX_grad = zeros(V,3*C);
%         rBY_grad = zeros(V,3*C);
%         rThX_grad = zeros(V,3*C);
%         rThY_grad = zeros(V,3*C);
        
        row = zeros(9*V,1);
        col = zeros(9*V,1);
        valBX = zeros(9*V,1);
        valBY = zeros(9*V,1);
        valThX = zeros(9*V,1);
        valThY = zeros(9*V,1);

        n = 1;
        for ii = 1:3 % Iterate over cells
           for jj = 1:3 % Iterate over dual dof
%                rBX_grad( (1:V)' + V*((jj-1)*C + tri(:,ii)-1) ) = dRb(:,1,ii,jj);
%                rBY_grad( (1:V)' + V*((jj-1)*C + tri(:,ii)-1) ) = dRb(:,2,ii,jj);
%                rThX_grad( (1:V)' + V*((jj-1)*C + tri(:,ii)-1) ) = dRTh(:,1,ii,jj);
%                rThY_grad( (1:V)' + V*((jj-1)*C + tri(:,ii)-1) ) = dRTh(:,2,ii,jj);
                II = (n-1)*V + 1:V;
                row(II) = 1:V;
                col(II) = (jj-1)*C + tri(:,ii);
                valBX(II) = dRb(:,1,ii,jj);
                valBY(II) = dRb(:,2,ii,jj);
                valThX(II) = dRTh(:,1,ii,jj);
                valThY(II) = dRTh(:,2,ii,jj);
           end
        end
        
        rBX_grad = sparse(row,col,valBX,V,3*C);
        rBY_grad = sparse(row,col,valBY,V,3*C);
        rThX_grad = sparse(row,col,valThX,V,3*C);
        rThY_grad =sparse(row,col,valThY,V,3*C);

%         rBX_grad = sparse(rBX_grad);
%         rBY_grad = sparse(rBY_grad);
%         rThX_grad = sparse(rThX_grad);
%         rThY_grad = sparse(rThY_grad);

        dceq = [];
        dc = (bsxfun(@times,bTh(:,1),d0) * rThX_grad + bsxfun(@times,bTh(:,2),d0) * rThY_grad) - ...
             (bsxfun(@times,bB(:,1),d0) * rBX_grad + bsxfun(@times,bB(:,2),d0) * rBY_grad);
        dc = dc';
        
    end
    
end

