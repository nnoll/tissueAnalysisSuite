% alphaVec = linspace(0,.1,10);
% alphaVec = [0,.01,.1];
alphaVec = 10;
T = 29;

clear c lambda sigma
for kk = 1:length(alphaVec)
%     ETmp = ERes{T}(3:21,3:21);
%     r0Tmp = r0{T}(3:21,3:21);
    [scale, goodTiles] = fitDual.rescaleGridDual(PN{T},ERes{T},alphaVec(kk));

%     sigma = zeros(size(PN2,1),size(PN2,2),3);
%     for ii = 1:size(PN2,1)
%         for jj = 1:size(PN2,2)
%             if ( ismember(ii + size(PN2,1)*(jj-1),goodTiles) )
%                 sigma(ii,jj,:) = PN2{ii,jj}.computeBoundaryStress(r0{T}{ii,jj});
%             end
%         end
%     end
% 
%     [ ~, ~, Stmp ] = measure.interpolateTensor( mean(xG,1), mean(yG,1), sigma, goodTiles );
%     [c(kk,:),lambda(kk,:)] = maxDelta( time_series, Stmp );
    
end