function [ PN, ERes ] = returnSubDual( Struct, mode, XRange, YRange )
    % RETURN DUAL. Wrapper function that will be used to fit the dual graph
    % of the observed lattice in Struct using 1 of 3 different fitting
    % algorithms: tension net, small pressure nets, and curvy nets.
    %
    % Inputs.
    % 1. Struct. Data structure containing cell array.
    % 2. mode. (1) TN (2) curvy nets.
    
    ERes = zeros(length(Struct),1);
    for t = 1:length(Struct)

        if (mode == 1)
           [ q, theta, tri, i0, ERes(t) ] = fitDual.subATN.returnSubDual( Struct(t), XRange, YRange );
           p = ones(size(theta));
        elseif (mode == 2)
           [ q, theta, p, tri, i0, ERes(t) ] = fitDual.subAFN.returnDual( Struct(t), XRange, YRange ); 
        end
        
        PN{t} = pressure.net( q, theta, p, tri, i0 );

    end

end

