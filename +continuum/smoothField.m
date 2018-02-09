function [ phiS ] = smoothField( phi, x, y, smSc )
%SMOOTHFIELD Summary of this function goes here
%   Detailed explanation goes here

    phiS = cell(length(phi),1);
    [x,y] = meshgrid(x,y);

    D = pdist2([x(:),y(:)],[x(:),y(:)]);
    L = exp(-D.^2/(2*smSc.^2));
    L = bsxfun(@rdivide,L,sum(L,2));

    for t = 1:length(phi)
        phiS{t} = zeros(size(phi{t}));
        dim = size(phi{t});
        nD = length(dim) - 2;
        if (nD == 1)
            for ii = 1:size(phi{t},1)
                smooth = L*squeeze(phi{t}(ii,:))';
                phiS{t}(ii,:,:) = reshape(smooth,[size(phiS{t},2),size(phiS{t},3)]);
            end
        elseif (nD == 2)
            for ii = 1:size(phi{t},1)
                for jj = 1:size(phi{t},2)
                    smooth = L*squeeze(phi{t}(ii,jj,:));
                    phiS{t}(ii,jj,:,:) = smooth;
                end
            end
        end
    end


end

