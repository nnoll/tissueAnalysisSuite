function [ sigmaI ] = rescaledStress( PN, Struct, xG, yG, ERes, r0, alpha )
%RESCALEDSTRESS Summary of this function goes here
%   Detailed explanation goes here

    Zc = 5;
    for t = 1:length(PN)
        if (~isempty(PN{t}))
            
            for ii = 1:(size(PN{t},1))
                for jj = 1:(size(PN{t},2))
                    Z(ii,jj) = measure.numberOfBulkCells(Struct(t),xG(:,ii),yG(:,jj));
                end
            end

            PN2 = fitDual.rescaleGridDual(PN{t},ERes{t},alpha);
            sigma = zeros(size(PN2,1),size(PN2,2),3);
            for ii = 1:(length(xG))
                for jj = 1:(length(yG))
                    if ( Z(ii,jj) >= Zc && ERes{t}(ii,jj) < 10 )
                        sigma(ii,jj,:) = PN2{ii,jj}.computeBoundaryStress(r0{t}{ii,jj});
                    end
                end
            end

            [ ~, ~, sigmaI{t} ] = measure.interpolateTensor( mean(xG,1), mean(yG,1), sigma, ERes{t} );
        end
    end
end

