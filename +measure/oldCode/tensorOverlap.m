function [ d, dM, c, lambda ] = tensorOverlap( myosin, stress )
    % TENSOR OVERLAP 
    
    % Compute correlation locally
    cO = myosin(:,:,1).*stress(:,:,1) + 2*myosin(:,:,2).*stress(:,:,2) + myosin(:,:,3).*stress(:,:,3);
    cM = myosin(:,:,1).*myosin(:,:,1) + 2*myosin(:,:,2).*myosin(:,:,2) + myosin(:,:,3).*myosin(:,:,3);
    cS = stress(:,:,1).*stress(:,:,1) + 2*stress(:,:,2).*stress(:,:,2) + stress(:,:,3).*stress(:,:,3);

    c = cO ./ sqrt( cM .* cS );
        
    % Compute global `difference' number.
    lambda = nanmean(cM(:)) ./ nanmean(cO(:));
    d = myosin - lambda*stress;
    dM = d(:,:,1).*d(:,:,1) + 2*d(:,:,2).*d(:,:,2) + d(:,:,3).*d(:,:,3);
    dM = dM ./ ( lambda * sqrt( nanmean(cM(:)) .* nanmean(cS(:)) ) );
    d = nanmean(dM(:));
    
end

