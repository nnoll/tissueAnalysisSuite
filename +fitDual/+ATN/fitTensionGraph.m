function [ Q, T, dC, dV, bulkVerts, bulkCells ] = fitTensionGraph( Struct, extCell )
    %FIT_TENSIONGRAPH Takes in a data structure at a single time point and fits
    %the lattice to the nearest tension net with average tension set to
    %one. Will do so via least squares. Set c.o.m of tension graph to be at
    %zero.

    %%Inputs
    %1. Struct - Data structure.
    %2. L - Watershed image.

    %%Outputs
    %1. T - Tension Graph
    %2. cAdj - cell adjacency matrix
    %3. tension - vector of found tensions.

    [ dC, dV, ~, bulkVerts, bulkCells ] = fitDual.ATN.computeDiffOperators( Struct, extCell );  
    
    nCells = size(dC,2);
    nBonds = size(dC,1);
    
    %Log vertex positions.
    rv = zeros(length(bulkVerts),2);
    for ii = 1:length(bulkVerts)
        rv(ii,1) = double(Struct.Vdat(bulkVerts(ii)).vertxcoord);
        rv(ii,2) = double(Struct.Vdat(bulkVerts(ii)).vertycoord);
    end
    
    %Calculate bonds
    rb = dV * rv;
    rb = bsxfun(@rdivide, rb, sqrt(sum(rb.^2,2)));
    
%     % Calculate length constraint matrix.
%     I = [dC'*dC,zeros(nCells,nCells);zeros(nCells,nCells),dC'*dC];
%     I = sparse(I);
    
    %% Generate initial guess for tension net by just taking cell centroid
    
    [ x0 ] = fitDual.ATN.seedTN( Struct, extCell );
%     x0 = zeros(nCells,2);
%     for ii = 1:length(bulkCells)
%        x0(ii,:) = Struct.Cdat(bulkCells(ii)).centroid.coord; 
%     end
    
    x0 = bsxfun(@minus,x0,mean(x0,1));
    x0 = x0 / mean( sqrt( sum( (dC*x0).^2, 2) ) );
%     x02 = bsxfun(@minus,x02,mean(x02,1));
%     x02 = x02 / mean( sqrt( sum( (dC*x02).^2, 2) ) );
    x0 = x0(:);
%     x02 = x02(:);
%     residuals(x0,sparse(dC),rb)
%     residuals(x02,sparse(dC),rb)
%     pause
    % Rescale lengths to be consistent with imposed constraint.
%     av_len = (x0'*I*x0)/nBonds;
%     x0 = x0/sqrt(av_len);

    Aeq = zeros(2,length(x0));
    Aeq(1,1:(length(x0)/2)) = 1;
    Aeq(2,(length(x0)/2+1):length(x0)) = 1;
    beq = [0;0];
    
    %% Perform fit
        % Initialize optimization parameters
%     if (size(dC,1) > 3000)
%         options = optimset('Algorithm','interior-point','Display','none','MaxIter',1e4,'MaxFunEvals',5e6,'GradObj','on','GradConstr','on');
%     else
%         options = optimset('Algorithm','interior-point','Display','none','MaxIter',1e2,'MaxFunEvals',5e6,...
%                   'GradObj','on','GradConstr','on'); %,'HessFcn',@(x,lambda) hessian(x,lambda,sparse(dC),rB));
%     end
% %     Q = fmincon(@(x)residuals(x,sparse(dC),rb),x0,[],[],Aeq,beq,[],[],@(x) constrain_length(x,I,nBonds),options);
%     tic
%     [Q,~,~,~,~,~,H_off] = fmincon(@(x)residuals(x,sparse(dC),rb),x0,[],[],Aeq,beq,[],[],@(x)constrain_length(x,sparse(dC)),options);
%     toc
    options = optimset('Algorithm','interior-point','Display','none','MaxIter',1e3,'MaxFunEvals',5e6,...
                  'GradObj','on','GradConstr','on','Hessian','on','HessFcn',@(x,lambda) hessian(x,lambda,sparse(dC),rb));
 
    Q = fmincon(@(x)residuals(x,sparse(dC),rb),x0,[],[],Aeq,beq,[],[],@(x)constrain_length(x,sparse(dC)),options);
  
%     imagesc(abs(H_off))
%     pause
%     imagesc(abs(H_on))
%     pause
    
    Q = reshape(Q,length(Q)/2,2);

    % Rescale T so that mean is set to one, not the second moment.
    T = sqrt(sum((dC*Q).^2,2));

end

function [ e, de ] = residuals(x, dC, rB)

    Q = reshape(x,length(x)/2,2);
    dQ = dC*Q;
    QBL = sqrt(sum(dQ.^2,2));
    dQ = bsxfun(@rdivide,dQ,QBL);
    
    IP = dot(dQ,rB,2);
    e = .5*mean(IP.^2);
    
    de = bsxfun(@times, IP./QBL, rB - bsxfun(@times,IP,dQ));
%     de = bsxfun(@times,IP,rB);
    de = (dC' * de) / size(dC,1);
    de = de(:);
    
end

function [ H ] = hessian(x,lambda,dC,rB)

    NB = size(dC,1);
    Q = reshape(x,length(x)/2,2);
    dQ = dC*Q;
    QBL = sqrt(sum(dQ.^2,2));
    dQ = bsxfun(@rdivide,dQ,QBL);

    IP = dot(dQ,rB,2);
    II = 1:size(dQ,1);

    C = size(Q,1);
%     Hobj = zeros(length(x),length(x));

    dF_x = sparse(II,II,(rB(:,1) - IP.*dQ(:,1))./QBL) * dC;
    dF_y = sparse(II,II,(rB(:,2) - IP.*dQ(:,2))./QBL) * dC;
    
    % Bilinear term
%     Hobj( 1:C, 1:C ) = dF_x' * dF_x;
%     Hobj( 1:C, (C+1):2*C ) = dF_x' * dF_y;
%     Hobj( (C+1):2*C, 1:C ) = dF_y' * dF_x;
%     Hobj( (C+1):2*C, (C+1):2*C ) = dF_y' *dF_y;

    Hobj = sparse( [dF_x' * dF_x, dF_x' * dF_y; dF_y' * dF_x, dF_y' *dF_y] );
    
    % Curvature term.
    dF_xx = dC' * sparse(II,II,IP.*((3*dQ(:,1).*dQ(:,1).*IP - rB(:,1).*dQ(:,1) - IP - rB(:,1).*dQ(:,1))./QBL.^2)) * dC;
    dF_xy = dC' * sparse(II,II,IP.*((3*dQ(:,1).*dQ(:,2).*IP - rB(:,1).*dQ(:,2) - rB(:,2).*dQ(:,1))./QBL.^2)) * dC;
    dF_yx = dC' * sparse(II,II,IP.*((3*dQ(:,2).*dQ(:,1).*IP - rB(:,2).*dQ(:,1) - rB(:,1).*dQ(:,2))./QBL.^2)) * dC;
    dF_yy = dC' * sparse(II,II,IP.*((3*dQ(:,2).*dQ(:,2).*IP - rB(:,2).*dQ(:,2) - IP - rB(:,2).*dQ(:,2))./QBL.^2)) * dC;

%     Hobj( 1:C, 1:C ) = Hobj( 1:C, 1:C ) + dF_xx;
%     Hobj( 1:C, (C+1):2*C ) = Hobj( 1:C, (C+1):2*C ) + dF_xy;
%     Hobj( (C+1):2*C, 1:C ) = Hobj( (C+1):2*C, 1:C ) + dF_yx;
%     Hobj( (C+1):2*C, (C+1):2*C ) = Hobj( (C+1):2*C, (C+1):2*C ) + dF_yy;
    Hobj = Hobj + sparse( [dF_xx, dF_xy; dF_yx, dF_yy] );
    
    % Nonlinear constraint.
    II = 1:size(dQ,1);
    Hcon = sparse([ [dC' * sparse(II,II,lambda.eqnonlin*(1-dQ(:,1).*dQ(:,1))./QBL) * dC, ...
              dC' * sparse(II,II,-lambda.eqnonlin*dQ(:,1).*dQ(:,2)./QBL) * dC];...
             [dC' * sparse(II,II,-lambda.eqnonlin*dQ(:,2).*dQ(:,1)./QBL) * dC,...
              dC' * sparse(II,II,lambda.eqnonlin*(1-dQ(:,2).*dQ(:,2))./QBL) * dC] ]);
%     Hcon( 1:C, 1:C ) = dC' * bsxfun(@times,(1-dQ(:,1).*dQ(:,1))./QBL,dC);
%     Hcon( 1:C, (C+1):2*C  ) = dC' * bsxfun(@times,(-dQ(:,1).*dQ(:,2))./QBL,dC);
%     Hcon( (C+1):2*C, 1:C ) = dC' * bsxfun(@times,(-dQ(:,2).*dQ(:,1))./QBL,dC);
%     Hcon( (C+1):2*C, (C+1):2*C ) = dC' * bsxfun(@times,(1-dQ(:,2).*dQ(:,2))./QBL,dC);

    H = Hobj/NB + Hcon/NB;
end

function [ c, ceq, dc, dceq ] = constrain_length(x,dC) %I,nbonds) 

    c = []; dc = [];
    Q = reshape(x,length(x)/2,2);
    dQ = dC*Q;
    QBL = sqrt(sum(dQ.^2,2));
    dQ = bsxfun(@rdivide,dQ,QBL);
    
    ceq = mean(QBL) - 1;
    dceq = (dC'*dQ) / size(dC,1);
    dceq = dceq(:);
    
%     ceq = x'*I*x - nbonds;
%     dceq = 2*I*x;
    
end