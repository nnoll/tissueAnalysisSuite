function [ q, theta, tri, cells, ERes ] = returnSubDual( Struct, XRange, YRange )
    % RETURN DUAL 

    % Store bulk vertices.
    [ tri, bulk0, ext0, bulkVerts, rC, r0 ] = fitDual.returnSubGraph( Struct, XRange, YRange );
    [ extVerts ] = fitDual.returnExtVerts( Struct, bulkVerts );

    if (isfield(Struct,'embed'))
        [ r0, ERes ] = generate.embedSubGraph( Struct, bulkVerts, extVerts, Struct.embed  );
    end
    optimset = optimoptions('fminunc','Display','none','Algorithm','quasi-newton','MaxFunEvals',5e6,'MaxIter',5e2,'TolFun',1e-3,'GradObj','on');

    if (nargin <= 3)
        if (isfield(Struct,'embed'))
            [ q, theta ] = fitDual.subATN.obtainSubStartPoint( Struct, XRange, YRange, r0 );
            r0 = r0(1:length(bulkVerts),:);
        else
            [ q, theta ] = fitDual.subATN.obtainSubStartPoint( Struct, XRange, YRange );
        end
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
    if (energyFunc(x0) > energyFunc([q(:);zeros(size(theta))]))
        x0 = [q(:);zeros(size(theta))];
    end

    if (energyFunc(x0) < 5e3)
        [x,ERes] = fminunc(energyFunc,x0,optimset);
        x = reshape(x,length(x)/3,3);
        q = x(:,1:2);
        theta = x(:,3);
    else
        ERes = 10;
    end
    
    cells = [bulk0,ext0];
    
end

