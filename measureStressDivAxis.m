delta = 5;
clear pStressAngle stressAngle longAngle oStressAngle tStressAngle
ii = 1;
ind = [];
for c = 1:length(mitCells)
    if ( length(mitCells(c).stressAxis) > delta-1 );
        tStressAngle(ii) = abs(cos(mitCellsT(c).stressAxis(end-delta+1) - mitCellsT(c).divAxis));
%         tStressAngle(ii) = mitCells(c).stressAxis(end-delta+1);
        longAngle(ii) = abs(cos(mitCells(c).longAxis(end-delta+1) - mitCells(c).divAxis));
        pStressAngle(ii) = abs(cos(mitCells(c).stressAxis(end-delta+1) - mitCells(c).divAxis));
        oStressAngle(ii) = abs(cos(mitCellsO(c).stressAxis(end-delta+1) - mitCellsO(c).divAxis));

%         tStressAngle(ii) = abs(cos(mitCellsT(c).stressAxis(end-delta+1) - mitCellsT(c).divAxis));
%         longAngle(ii) = abs(cos(mitCells(c).longAxis(end-delta+1) - mitCells(c).divAxis));
%         pStressAngle(ii) = abs(cos(mitCells(c).stressAxis(end-delta+1) - mitCells(c).divAxis));
%         oStressAngle(ii) = abs(cos(mitCellsO(c).stressAxis(end-delta+1) - mitCellsO(c).divAxis));
        ind = [ind,c];
        ii = ii + 1;
    end
end