function [ PN, ERes, r0 ] = returnDual( Struct, mode, extCell )
    % RETURN DUAL. Wrapper function that will be used to fit the dual graph
    % of the observed lattice in Struct using 1 of 3 different fitting
    % algorithms: tension net, small pressure nets, and curvy nets.
    %
    % Parameters
    % ----------
    % Struct. Data structure containing cell array.
    % mode : 
    %   (1) Tension Net (ie straight line net) 
    %   (2) curvy nets to compute forces with pressure differences.
    % extCell : set to 1
    %   Boundary conditions.
    %
    % Returns
    % -------
    % PN : cell of classes
    %   

    if (nargin == 2)
        % By default, the external cell is labeled as cell 1
        extCell = 1;
    end
    
    ERes = zeros(length(Struct),1);
    for t = 1:length(Struct)
        disp(['Fitting dual for time point ', num2str(t), ' / ', num2str(length(Struct))])
        if (mode == 1)
           [ q, theta, tri, i0, i2, ERes(t) ] = fitDual.ATN.returnDual( Struct(t), extCell );
           p = ones(size(theta));
        elseif (mode == 2)
           tic
           [ q, theta, p, tri, i0, i2, ERes(t) ] = fitDual.AFN.returnDual( Struct(t), extCell ); 
           toc
        end
        
        disp('   -> updating pressure net')
        PN{t} = pressure.net( q, theta, p, tri, i0, i2 );
        if (nargout == 3)
            [ ~, ~, ~, ~, ~, r0{t} ] = fitDual.returnGraph(Struct(t), extCell);
        end

    end

end

