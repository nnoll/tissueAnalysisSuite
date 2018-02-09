function [ ] = curvedTension( Struct, mode )
%INFERCURVE Summary of this function goes here
%   Detailed explanation goes here

    T = zeros(length(Struct.Bdat),1);
    dV = zeros(length(Struct.Bdat),2);

    badBonds = [];
    goodBonds = 1:length(Struct.Bdat);
    for b = 1:length(Struct.Bdat)
        if (~isempty(Struct.Bdat(b).tension))
            if (nargin == 1 || mode == 0 )
                T(b) = Struct.Bdat(b).tension;
            else
                T(b) = Struct.Bdat(b).chem(1);
            end
            dV(b,:) = Struct.Bdat(b).verts;
        else
            badBonds = [badBonds,b];
        end
    end
        
    T(badBonds) = [];
    goodBonds(badBonds) = [];
    dV(badBonds,:) = [];
    
    Tmax = prctile(T,90);
    Tmin = prctile(T,1);

    T = (T-Tmin)/(Tmax-Tmin);
    T(T>1) = 1;
    T(T<0) = 0;
    
    % Convert tension to color.
    cmap = hot(256);
    x = linspace(0,1,256);
    
    Tcolor(:,1) = interp1(x,cmap(:,1),T);
    Tcolor(:,2) = interp1(x,cmap(:,2),T);
    Tcolor(:,3) = interp1(x,cmap(:,3),T);
    
    hold on
    for b = 1:size(dV,1)
        bond = goodBonds(b);
        r = double([[Struct.Vdat(Struct.Bdat(bond).verts).vertxcoord]; ...
                    [Struct.Vdat(Struct.Bdat(bond).verts).vertycoord]]);
        r = r';
        rho = Struct.Bdat(bond).rBar;
        S = sign( (r(1,1) - rho(1))*(r(2,2) - rho(2)) - (r(1,2) - rho(2))*(r(2,1) - rho(1)) );
        if ( sign(S) < 0 )
            r = r([2,1],:);
        end
        R = Struct.Bdat(bond).radius;

        theta = atan2(r(:,2)-rho(2),r(:,1)-rho(1));
        theta(theta<0) = theta(theta<0) + 2*pi;
        theta = sort(theta);
        if (theta(2) - theta(1) > pi)
            theta(2) = theta(2)-2*pi;
            theta = theta([2,1]);
        end
        theta = linspace(theta(1),theta(2),20);
        x = rho(1) + R*cos(theta);
        y = rho(2) + R*sin(theta);
        plot(x,y,'LineWidth',2,'Color',Tcolor(b,:))
    end
    
end

