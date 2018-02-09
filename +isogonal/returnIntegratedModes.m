function [ thetaInt ] = returnIntegratedModes( theta, bulkCells, cPair, tStart )
%RETURNINTEGRATEDMODES 

    thetaInt = theta{tStart};
    lastInd = 1:length(bulkCells{tStart});
    for t = (tStart+1):length(theta)
        trackedCells = cPair{t}(bulkCells{t-1}(lastInd));
        lastInd = find(ismember(trackedCells,bulkCells{t}));
        for ii = 1:length(lastInd)
            thetaInt(lastInd(ii)) = theta{t}(bulkCells{t}==trackedCells(lastInd(ii)));
        end
    end

end

