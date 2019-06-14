function [ ] = tension( Struct, mode )
% TENSION  Plot the bonds colored by the tension in each bond
% 
% Parameters
% ----------
% Struct : struct
%   The segmentation and mechanics for this experiment
% 
% mode : (0 or 1) If zero, plots Struct.Bdat.tension, but if nonzero plots
% Struct.Bdat.actual_tension

    if (nargin == 1)
        mode = 0;
    end
    
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
