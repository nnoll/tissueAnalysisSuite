function [ ] = bondlengths( Struct, varargin )
% BONDLENGTHS  Plot the bonds colored by the length of each bond. This is
% unfinished. Please finish writing it!
% 
% Parameters
% ----------
% Struct : struct
%   The segmentation and mechanics for this experiment
% 
% mode : (0 or 1) unused currently

    if ~isempty(varargin)
        % a variable input argument has been passed! Expect the bondlengths
        % which have already been computed to be the variable.
        bondlengths = varargin ;
    end
    
    % Collate bonds into linesegs
    xyxy = zeros(length(bdat), 4) ;
    vectors = zeros(length(bdat), 2) ;
    bondxy = zeros(length(bdat), 2) ;
    for bondi = 1:length(bdat)
        ij = bdat(bondi).verts ;
        x1 = vdat(ij(1)).vertxcoord ;
        y1 = vdat(ij(1)).vertycoord ;
        x2 = vdat(ij(2)).vertxcoord ;
        y2 = vdat(ij(2)).vertycoord ;
        xyxy(bondi, 1:2) = [x1, y1];
        xyxy(bondi, 3:4) = [x2, y2];
        vectors(bondi, :) = xyxy(bondi, 3:4) - xyxy(bondi, 1:2) ;
        bondxy(bondi, :) = [0.5 * (x1 + x2), 0.5 * (y1 + y2)] ;
    end
    
    break
    
    T = zeros(length(Struct.Bdat),1);
    dV = zeros(length(Struct.Bdat),2);

    badBonds = [];
    if (mode == 0)
        for b = 1:length(Struct.Bdat)
            if (~isempty(Struct.Bdat(b).tension))
                T(b) = Struct.Bdat(b).tension;
                dV(b,:) = Struct.Bdat(b).verts;
            else
                badBonds = [badBonds,b];
            end
        end
    else
        for b = 1:length(Struct.Bdat)
            if (~isempty(Struct.Bdat(b).actual_tension))
                T(b) = Struct.Bdat(b).actual_tension;
                dV(b,:) = Struct.Bdat(b).verts;
            else
                badBonds = [badBonds,b];
            end
        end
    end
    
    T(badBonds) = [];
    dV(badBonds,:) = [];
    
    % Get limits of tensions to plot
    Tmax = prctile(T,98);
    Tmin = prctile(T,2);

    T = (T-Tmin)/(Tmax-Tmin);
    T(T>1) = 1;
    T(T<0) = 0;
    
    % Convert tension to color.
    cmap = hot(256);
    x = linspace(0,1,256);
    Tcolor(:,1) = interp1(x,cmap(:,1),T);
    Tcolor(:,2) = interp1(x,cmap(:,2),T);
    Tcolor(:,3) = interp1(x,cmap(:,3),T);
    
    % Plot the tensions
    hold on
    Vlist = Struct.Vdat;
    for b = 1:size(dV,1)
        v = dV(b,1);
        nv = dV(b,2);
        plot([Vlist(v).vertxcoord,Vlist(nv).vertxcoord],[Vlist(v).vertycoord,Vlist(nv).vertycoord],'Color',Tcolor(b,:),'LineWidth',2);
    end
    
end