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
    
    % Formulate inverse as a minimization problem.
    G = [bsxfun(@times,rb(:,1),dC),bsxfun(@times,rb(:,2),dC)];

    % Fix center of mass of tension graph
    G = [G; [ones(1,nCells),zeros(1,nCells)]; [zeros(1,nCells),ones(1,nCells)]];
    G = sparse(G);
    
    % Calculate length constraint matrix.
    I = [dC'*dC,zeros(nCells,nCells);zeros(nCells,nCells),dC'*dC];
    
    % Initialize optimization parameters
    options = optimset('Algorithm','active-set','Display','iter','MaxIter',1e3,'GradObj','on','GradConstr','on');
    
    %% Generate initial guess for tension net by just taking cell centroid.

    x0 = zeros(2,nCells);
    for ii = 1:nCells
        x0(:,ii) = Struct.Cdat(bulkCells(ii)).centroid.coord';
    end
%     [ x0 ] = fitDual.ATN.seedTN( Struct, extCell );
%     x0 = x0';
    x0 = x0 - mean(x0,2)*ones(1,nCells);
    x0 = reshape(x0,2*nCells,1);

    % Rescale lengths to be consistent with imposed constraint.
    av_len = (x0'*I*x0)/nBonds;
    x0 = x0/sqrt(av_len);

    %% Perform fit
    Q = fmincon(@(x)residuals(x,G),x0,[],[],[],[],[],[],@(x)constrain_length(x,I,nBonds),options);
    Q = reshape(Q,length(Q)/2,2);

    % Rescale T so that mean is set to one, not the second moment.
    T = sqrt(sum((dC*Q).^2,2));

end

function [ e, de ] = residuals(x, G)
    e = x'*(G'*G)*x;
    de = 2*(G'*G)*x;
end

function [c, ceq, dc, dceq] = constrain_length(x,I,nbonds)
    c = []; dc = [];
    ceq = x'*I*x - nbonds;
    dceq = 2*I*x;
end