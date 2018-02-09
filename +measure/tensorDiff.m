function [ D ] = tensorDiff( A, B )
    % TENSOR OVLAP
    
    lambda = sqrt(mean( sum(sum((A.*A),1),2) )) / sqrt(mean( sum(sum((B.*B),1),2) ));
    B = lambda*B;
    delta = A-B;

    % Compute normalized distance.
    D = mean( sum(sum((delta.*delta),1),2) );
    D = D / (mean( sqrt(sum(sum(A.*A,1),2)) ) * mean( sqrt(sum(sum(B.*B,1),2)) ));
%     D = sqrt(D);
    
end

