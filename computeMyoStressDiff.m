clear diff
% alphaVec = linspace(0,10,30);
alphaVec = 1.7;
deltaVec = 0:20;
for ii = 1:length(alphaVec)
    [ expMyo ] = measure.embryoMyo( time_series, alphaVec(ii) );
    for jj = 1:length(deltaVec)
        for t = timePts
            [ D ] = continuum.tensorDiff( expMyo{t+deltaVec(jj)}, sigma{t} );
            diff(t,jj,ii) = nanmean(D(:));
        end
    end
end
