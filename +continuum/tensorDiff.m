function [ D, DTrace ] = tensorDiff( expMyo, sigmaI, mode )
    % TENSOR DIFF
    
    if ( mode == 1 )

        lambdaN = expMyo(:,:,1).*expMyo(:,:,1) + 2*expMyo(:,:,2).*expMyo(:,:,2) + expMyo(:,:,3).*expMyo(:,:,3);
        lambdaD = sigmaI(:,:,1).*sigmaI(:,:,1) + 2*sigmaI(:,:,2).*sigmaI(:,:,2) + sigmaI(:,:,3).*sigmaI(:,:,3);
        
        lambda = sqrt(nanmean(lambdaN(:))) / sqrt(nanmean(lambdaD(:)));
        sigmaI = lambda*sigmaI;
        diff = expMyo - sigmaI;
        
        D = diff(:,:,1).*diff(:,:,1) + 2*diff(:,:,2).*diff(:,:,2) + diff(:,:,3).*diff(:,:,3);
        D = D ./ (expMyo(:,:,1).*expMyo(:,:,1) + 2*expMyo(:,:,2).*expMyo(:,:,2) + expMyo(:,:,3).*expMyo(:,:,3));

    else
        
        traceMyo = expMyo(:,:,1) + expMyo(:,:,3);
        traceSig = sigmaI(:,:,1) + sigmaI(:,:,3);
        
        expMyo(:,:,1) = expMyo(:,:,1) - .5*traceMyo;
        expMyo(:,:,3) = expMyo(:,:,3) - .5*traceMyo;
        sigmaI(:,:,1) = sigmaI(:,:,1) - .5*traceSig;
        sigmaI(:,:,3) = sigmaI(:,:,3) - .5*traceSig;
        
        lambdaN = expMyo(:,:,1).*expMyo(:,:,1) + 2*expMyo(:,:,2).*expMyo(:,:,2) + expMyo(:,:,3).*expMyo(:,:,3);
        lambdaD = sigmaI(:,:,1).*sigmaI(:,:,1) + 2*sigmaI(:,:,2).*sigmaI(:,:,2) + sigmaI(:,:,3).*sigmaI(:,:,3);
        
        lambda = sqrt(nanmean(lambdaN(:))) / sqrt(nanmean(lambdaD(:)));
        sigmaI = lambda*sigmaI;
        diff = expMyo - sigmaI;
        
        D = diff(:,:,1).*diff(:,:,1) + 2*diff(:,:,2).*diff(:,:,2) + diff(:,:,3).*diff(:,:,3);
        D = D ./ (expMyo(:,:,1).*expMyo(:,:,1) + 2*expMyo(:,:,2).*expMyo(:,:,2) + expMyo(:,:,3).*expMyo(:,:,3));
        
        stdM = nanstd(traceMyo(:));
        stdS = nanstd(traceSig(:));
        meanM = nanmean(traceMyo(:));
        meanS = nanmean(traceSig(:));
        DTrace = ((traceMyo-meanM).*(traceSig-meanS)./(stdM.*stdS));
%         scatter(traceMyo(:),traceSig(:))
%         pause
        
    end

end

