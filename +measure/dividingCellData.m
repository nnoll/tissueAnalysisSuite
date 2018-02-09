function [ mitCells ] = dividingCellData( mitCells, Struct, L )
    % DIVIDING CELL DATA 
    
    for c = 1:length(mitCells)
        
        % Tracked times
        delta = length(mitCells(c).label);
        tS = mitCells(c).t0 - delta + 1;
        t0 = mitCells(c).t0;
        
        % Division Axis
        PixelIdxList = L(:,:,t0) == mitCells(c).label(end);
        iCells = L(:,:,t0+1);
        iCells = iCells(PixelIdxList);
        iCells = iCells(iCells>1);
        daughter1 = mode(iCells);
        daughter2 = mode(iCells(iCells~=daughter1));
        
        axis =  (L(:,:,t0+1) == daughter1) + (L(:,:,t0+1) == daughter2);
        rgb(:,:,1) = axis;
        rgb(:,:,2) = L(:,:,t0) == mitCells(c).label(end);
        axis = bwmorph(axis,'bridge') - axis;
        rgb(:,:,3) = 0; %axis;

        [y,x] = find(axis);
        alpha = mean( (y-mean(y)) .* (x-mean(x)) ) / mean( (x-mean(x)) .* (x-mean(x)) );
        mitCells(c).divAxis = atan(alpha) + pi/2;        
%         beta = mean(y) - alpha*mean(x);

        % Long Axis 
        mitCells(c).longAxis = zeros(delta,1);
        
        % Area
        mitCells(c).area = zeros(delta,1);

        % Eccentricity
        mitCells(c).ecc = zeros(delta,1);
        
        % Pressure
        mitCells(c).pressure = zeros(delta,1);

        % Stress
        mitCells(c).stressAxis = zeros(delta,1);
        mitCells(c).stressAnis = zeros(delta,1);

        ii = 1;
        for t = tS:t0
            S = regionprops(L(:,:,t),'Area','Orientation','MajorAxisLength','MinorAxisLength');
            
            % Long Axis
            mitCells(c).longAxis(ii) = -S(mitCells(c).label(ii)).Orientation;
            if (mitCells(c).longAxis(ii) < 0)
                mitCells(c).longAxis(ii) = mitCells(c).longAxis(ii) + 180; 
            end
            mitCells(c).longAxis(ii) = pi*mitCells(c).longAxis(ii)/180;
            
            % Eccentricity
            mitCells(c).ecc(ii) = (S(mitCells(c).label(ii)).MajorAxisLength-S(mitCells(c).label(ii)).MinorAxisLength)/S(mitCells(c).label(ii)).MajorAxisLength;

            % Area
            mitCells(c).area(ii) = S(mitCells(c).label(ii)).Area;

            % Pressure
%             if (~isempty(Struct(t).Cdat(mitCells(c).label(ii)).pressure))
%                 P = mean([Struct(t).Cdat.pressure]);
%                 mitCells(c).pressure(ii) = Struct(t).Cdat(mitCells(c).label(ii)).pressure/P;
%             else
%                 mitCells(c).pressure(ii) = 1;
%             end
            
            % Stress
            if (~isempty(Struct(t).Cdat(mitCells(c).label(ii)).smooth_stress))
                P = mean([Struct(t).Cdat.pressure]);
                stress = Struct(t).Cdat(mitCells(c).label(ii)).smooth_stress;
                [V,D] = eig(stress);
                mitCells(c).stressAxis(ii) = atan(V(2,1)/V(1,1));
                mitCells(c).stressAnis(ii) = abs(D(2,2) - D(1,1))/P;
                mitCells(c).pressure(ii) = (D(2,2) + D(1,1))/P;
            end
            ii = ii + 1;
        end
        
    end


end

