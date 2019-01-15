function [ stress, xVec, yVec ] = coarseGrain( scale, PN, xG, yG, ERes )
    % TAKE CONTINUUM LIMIT 
    
    stress = cell(size(PN));
    for t = 1:length(PN)
        stress{t} = zeros(size(PN{t},1),size(PN{t},2),3);

        for ii = 1:size(yG{t},2)
            ii
            for jj = 1:size(xG{t},2)
                if (ERes{t}(ii,jj) < 10) 
                    [ stress{t}(ii,jj,:) ] = PN{t}{ii,jj}.computeBoundaryStress();
                    stress{t}(ii,jj,:) = scale{t}(ii,jj) * stress{t}(ii,jj,:);
                end
            end
        end       
    end
end

