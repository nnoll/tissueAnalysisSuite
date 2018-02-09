function [ PN, ERes, r0 ] = returnDual( Struct, mode, extCell )
    % RETURN DUAL. Wrapper function that will be used to fit the dual graph
    % of the observed lattice in Struct using 1 of 3 different fitting
    % algorithms: tension net, small pressure nets, and curvy nets.
    %
    % Inputs.
    % 1. Struct. Data structure containing cell array.
    % 2. mode. (1) TN (2) small pressure net (3) curvy nets.
    % 3. Text. Boundary conditions.

    if (nargin == 2)
        extCell = 0;
    end
    
    ERes = zeros(length(Struct),1);
    for t = 1:length(Struct)
        t
        if (mode == 1)
           [ q, theta, tri, i0, i2, ERes(t) ] = fitDual.ATN.returnDual( Struct(t), extCell );
           p = ones(size(theta));
        elseif (mode == 2)
           [ q, theta ] = fitDual.ATN.returnDual( Struct(t), extCell );
           p = ones(size(theta));
           Q = [q,theta,p];
           [ q, theta, p, tri, i0, ERes(t) ] = fitDual.AFN.returnDual( Struct(t), extCell, Q ); 
        elseif (mode == 3)
           tic
           [ q, theta, p, tri, i0, i2, ERes(t) ] = fitDual.AFN.returnDual( Struct(t), extCell ); 
           toc
        end
        
        PN{t} = pressure.net( q, theta, p, tri, i0, i2 );
        if (nargout == 3)
            [ ~, ~, ~, ~, ~, r0{t} ] = fitDual.returnGraph( Struct(t), extCell );
        end

    end

end

