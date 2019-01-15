function [ ovLap ] = tensorOvlap( expMyo, sigmaI, mode )
    % TENSOR OVLAP 
    
    if ( nargin < 3)
        mode = 1;
    end
    
    if ( mode == 1 )
%         ovLap = expMyo(:,:,1).*sigmaI(:,:,1) + 2*expMyo(:,:,2).*sigmaI(:,:,2) + expMyo(:,:,3).*sigmaI(:,:,3);
%         ovLap = ovLap ./ sqrt(expMyo(:,:,1).*expMyo(:,:,1) + 2*expMyo(:,:,2).*expMyo(:,:,2) + expMyo(:,:,3).*expMyo(:,:,3));
%         ovLap = ovLap ./ sqrt(sigmaI(:,:,1).*sigmaI(:,:,1) + 2*sigmaI(:,:,2).*sigmaI(:,:,2) + sigmaI(:,:,3).*sigmaI(:,:,3));
        ovLap = zeros(size(expMyo,1),size(expMyo,2));
          for ii = 1:size(expMyo,1)
              for jj = 1:size(expMyo,2)
                  A = [expMyo(ii,jj,1),expMyo(ii,jj,2);expMyo(ii,jj,2),expMyo(ii,jj,3)];
                  B = [sigmaI(ii,jj,1),sigmaI(ii,jj,2);sigmaI(ii,jj,2),sigmaI(ii,jj,3)];
                  if (sum(isnan(A(:))) == 0 && sum(isnan(B(:))) == 0)
                      [Va,Da] = eig(A);
                      [Vb,Db] = eig(B);
                      ovLap(ii,jj) = Va(:,2)' * Vb(:,2);
                  else
                      ovLap(ii,jj) = nan;
                  end
              end
          end
%         ovLap = ovLap ./ (2*(expMyo(:,:,1).*expMyo(:,:,3) - expMyo(:,:,2).*expMyo(:,:,2));
    else
        ovLap = expMyo(:,:,1).*sigmaI(:,:,3) - 2*expMyo(:,:,2).*sigmaI(:,:,2) + expMyo(:,:,3).*sigmaI(:,:,1);
        ovLap = ovLap ./ sqrt(expMyo(:,:,1).*expMyo(:,:,1) + 2*expMyo(:,:,2).*expMyo(:,:,2) + expMyo(:,:,3).*expMyo(:,:,3));
        ovLap = ovLap ./ sqrt(sigmaI(:,:,1).*sigmaI(:,:,1) + 2*sigmaI(:,:,2).*sigmaI(:,:,2) + sigmaI(:,:,3).*sigmaI(:,:,3));
    end

end

