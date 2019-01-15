function [ Q, T, dC, dV, bulkVerts, bulkCells ] = fitSubTensionGraph( Struct, XRange, YRange, rv )
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

    [ dC, dV, ~, bulkVerts, bulkCells ] = fitDual.subATN.computeSubDiffOperators( Struct, XRange, YRange );  
    
    nCells = size(dC,2);
    nBonds = size(dC,1);
    
    %Log vertex positions.
    if (nargin < 4)
        rv = zeros(length(bulkVerts),2);
        for ii = 1:length(bulkVerts)
            rv(ii,1) = double(Struct.Vdat(bulkVerts(ii)).vertxcoord);
            rv(ii,2) = double(Struct.Vdat(bulkVerts(ii)).vertycoord);
        end
    end
    
%     plot.skel(Struct,'b',0)
%     hold all
%     scatter(rv(:,1),rv(:,2),'r','filled')
%     pause
    
    %Calculate bonds
    rb = dV * rv;
%     rv([17,20],:)
    rb = bsxfun(@rdivide, rb, sqrt(sum(rb.^2,2)));
    
%     % Calculate length constraint matrix.
%     I = [dC'*dC,zeros(nCells,nCells);zeros(nCells,nCells),dC'*dC];
%     I = sparse(I);
    
    % Initialize optimization parameters
    options = optimset('Algorithm','active-set','Display','none','MaxIter',1e4,'MaxFunEvals',5e6,'GradObj','on','GradConstr','on');
    
    %% Generate initial guess for tension net by just taking cell centroid.

    [ x0 ] = fitDual.subATN.seedSubTN( Struct, XRange, YRange );
%     x0 = zeros(nCells,2);
%     for ii = 1:length(bulkCells)
%        x0(ii,:) = Struct.Cdat(bulkCells(ii)).centroid.coord; 
%     end
    x0 = bsxfun(@minus,x0,mean(x0,1));
    x0 = x0 / mean( sqrt( sum( (dC*x0).^2, 2) ) );
    x0 = x0(:);
%     x02 = bsxfun(@minus,x02,mean(x02,1));
%     x02 = x02 / mean( sqrt( sum( (dC*x02).^2, 2) ) );
%     x02 = x02(:);
%     x0 = x02;
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
%     Q = fmincon(@(x)residuals(x,sparse(dC),rb),x0,[],[],Aeq,beq,[],[],@(x) constrain_length(x,I,nBonds),options);
    Q = fmincon(@(x)residuals(x,sparse(dC),rb),x0,[],[],Aeq,beq,[],[],@(x)constrain_length(x,sparse(dC)),options);
    Q = reshape(Q,length(Q)/2,2);
%     [ tri ] = fitDual.returnSubGraph( Struct, XRange, YRange );
%     plot.skel(Struct,'r',0)
%     hold all
%     patch('Faces',tri,'Vertices',reshape(x0,length(x0)/2,2),'FaceColor','none','LineWidth',2,'EdgeColor','b');
%     patch('Faces',tri,'Vertices',Q,'FaceColor','none','LineWidth',2);
%     pause
    % Rescale T so that mean is set to one, not the second moment.
    T = sqrt(sum((dC*Q).^2,2));

end

function [ e, de ] = residuals(x, dC, rB)

    Q = reshape(x,length(x)/2,2);
    dQ = dC*Q;
%     QBL = sqrt(sum(dQ.^2,2));
%     dQ = bsxfun(@rdivide,dQ,QBL);
    
    IP = dot(dQ,rB,2);
    e = .5*mean(IP.^2);
    
%     de = bsxfun(@times, IP./QBL, rB - bsxfun(@times,IP,dQ));
    de = bsxfun(@times,IP,rB);
    de = (dC' * de) / size(dC,1);
    de = de(:);
    
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