function [ d0, bCells, rBX, rBY, mask, iCells ] = getCurvedGraph( Struct, extCell )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    [ d0, ~, ~, ~, iCells ] = fitDual.ATN.computeDiffOperators( Struct, extCell );  
    [ involvedBonds ] = generate.bondMap(Struct);
    involvedBonds = involvedBonds{1};

    bCells = zeros(length(involvedBonds),2);
    badBonds = [];
    nPix = zeros(length(involvedBonds),1);
    for b = 1:length(involvedBonds)
        if (involvedBonds(b) > 0)
            c1 = find(d0(b,:)==1);
            c2 = find(d0(b,:)==-1);
            nPix(b) = length(Struct.Bdat(involvedBonds(b)).pix);
%             if (length(Struct.Bdat(involvedBonds(b)).pix) > z)
%             z = length(Struct.Bdat(involvedBonds(b)).pix);
%             end
            bCells(b,:) = [c1,c2];
        else
            badBonds = [badBonds,b];
        end
    end
    z
    involvedBonds(badBonds) = [];
    bCells(badBonds,:) = [];
    d0(badBonds,:) = [];
    nPix(badBonds) = [];
    
    z = round(median(nPix));
    rBX = zeros(length(involvedBonds),z+2);
    rBY = zeros(length(involvedBonds),z+2);
    for ii = 1:length(involvedBonds)
%         ind = Struct.Bdat(involvedBonds(ii)).pix;
%         [y,x] = ind2sub(dim,ind);
%         x = double(x);
%         y = double(y);
%         v1 = Struct.Bdat(involvedBonds(ii)).verts(1);
%         v2 = Struct.Bdat(involvedBonds(ii)).verts(2);
%         rBX(ii,1:(length(x)+2)) = [double(x);Struct.Vdat(v1).vertxcoord;Struct.Vdat(v2).vertxcoord];
%         rBY(ii,1:(length(x)+2)) = [double(y);Struct.Vdat(v1).vertycoord;Struct.Vdat(v2).vertycoord];
        if (length(Struct.Bdat(b).verts) == 2 && sum(Struct.Bdat(b).verts > 0) == 2)
           if (Struct.Bdat(b).radius < inf)
               r = double([[Struct.Vdat(Struct.Bdat(b).verts).vertxcoord]; ...
                    [Struct.Vdat(Struct.Bdat(b).verts).vertycoord]]);
               r = r';
               rho = Struct.Bdat(b).rBar;
               S = sign( (r(1,1) - rho(1))*(r(2,2) - rho(2)) - (r(1,2) - rho(2))*(r(2,1) - rho(1)) );
               if ( sign(S) < 0 )
                   r = r([2,1],:);
               end
               R = Struct.Bdat(b).radius;

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
               plot(x,y,'LineWidth',2,'Color',color)
 
           else
               r = [[Struct.Vdat(Struct.Bdat(b).verts).vertxcoord]; ...
                    [Struct.Vdat(Struct.Bdat(b).verts).vertycoord]];
               r = r';
               plot(r(:,1),r(:,2),'LineWidth',2,'Color',color)
           end
    end

    mask = (rBX ~= 0);
end

