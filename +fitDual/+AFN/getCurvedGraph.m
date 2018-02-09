function [ d0, bCells, rBX, rBY, mask, iCells ] = getCurvedGraph( Struct, extCell )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    dim = size(Struct.labelMat);
    [ d0, ~, ~, ~, iCells ] = fitDual.ATN.computeDiffOperators( Struct, extCell );  
    [ involvedBonds ] = generate.bondMap(Struct);
    involvedBonds = involvedBonds{1};

    z = 0;
    bCells = zeros(length(involvedBonds),2);
    badBonds = [];
    nPix = zeros(length(involvedBonds),1);
    for b = 1:length(involvedBonds)
        if (involvedBonds(b) > 0)
            c1 = find(d0(b,:)==1);
            c2 = find(d0(b,:)==-1);
            nPix(b) = length(Struct.Bdat(involvedBonds(b)).pix);
            bCells(b,:) = [c1,c2];
        else
            badBonds = [badBonds,b];
        end
    end
    
    involvedBonds(badBonds) = [];
    bCells(badBonds,:) = [];
    d0(badBonds,:) = [];
    nPix(badBonds) = [];
    
    z = round(median(nPix));
    rBX = zeros(length(involvedBonds),z);
    rBY = zeros(length(involvedBonds),z);
    for ii = 1:length(involvedBonds)
       b = involvedBonds(ii);
       
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
           theta = linspace(theta(1),theta(2),z);
           rBX(ii,:) = rho(1) + R*cos(theta);
           rBY(ii,:) = rho(2) + R*sin(theta);
       else
           r1 = double([Struct.Vdat(Struct.Bdat(b).verts(1)).vertxcoord; ...
                 Struct.Vdat(Struct.Bdat(b).verts(1)).vertycoord]);
           r2 = double([Struct.Vdat(Struct.Bdat(b).verts(2)).vertxcoord; ...
                 Struct.Vdat(Struct.Bdat(b).verts(2)).vertycoord]);    
           delta = r2 - r1;
           s = linspace(0,1,z);
           rBX(ii,:) = r1(1) + s*delta(1);
           rBY(ii,:) = r1(2) + s*delta(2);
       end
       
    end

    mask = (rBX ~= 0);
    
end

