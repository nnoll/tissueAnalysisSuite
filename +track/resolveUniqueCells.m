function [ cells ] = resolveUniqueCells( Struct, cPair )
    %RESOLVE UNIQUE CELLS 

    cells(1) = struct('label',[],'t0',[]);

    for c = 1:length(Struct(1).Cdat)
        cells(c).label = c;
        cells(c).t0 = 1;
    end
    
    lastID0 = 1:length(Struct(1).Cdat);
    for t = 2:length(cPair)
        trackedLabels = find(cPair{t} > 0);
        label0 = cPair{t}(trackedLabels);
        
        newLabels = zeros(length(cells),1);
        newLabels(lastID0(label0)) = trackedLabels;
        
        for c = find(newLabels)'
            cells(c).label = [cells(c).label,newLabels(c)];
        end
        
        lastID = zeros(length(Struct(t).Cdat),1);
        lastID(trackedLabels) = lastID0(label0);
        
        % Add new ids
        C = length(cells);
        n = C + 1;
        for utL = find(cPair{t} == 0)'
            cells(n).label = utL;
            cells(n).t0 = t;
            lastID(utL) = n;
            n = n + 1;
        end
        lastID0 = lastID;
    end

end

