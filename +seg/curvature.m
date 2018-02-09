function [ Struct ] = curvature( Struct, dim )
%EXTRACT_CURVATURE. This algorithm takes in a data structure WITH bond
%information, and computes the curvature of each bond in the image. Data is
%stored and then output in an updated data structure.

%%Inputs
%1. Struct - data structure

%%Outputs
%1. bStruct - data structure with curvature

for t = 1:length(Struct)

    for b = 1:length(Struct(t).Bdat)
        v = Struct(t).Bdat(b).verts(1); 
        nv = Struct(t).Bdat(b).verts(2);
        
        % Get vertex positions and bond vector.
        rv = double([Struct(t).Vdat(v).vertxcoord; Struct(t).Vdat(v).vertycoord]);
        rnv = double([Struct(t).Vdat(nv).vertxcoord; Struct(t).Vdat(nv).vertycoord]);

        % Dilate pixels associated with watershedded edge.
        Pix = Struct(t).Bdat(b).pix;
%         Pix = horzcat(Pix', [rv(2) + dim(1)*(rv(1)-1),rnv(2) + dim(1)*(rnv(1)-1)]);
        [ dilPix ] = dilatePix(Pix, 0, dim);
        
        %Put bond pixs into a vector.
        [by, bx] = ind2sub(dim,dilPix);                
        B = vertcat(bx,by);
        
        %Fit circle to observed data points.
        nB = [0,1;-1,0] * (rv - rnv);
        D = sqrt(sum(nB.^2));
        nB = nB / D;
        x0 = .5*(rv + rnv);
% 
%         Par = seg.CircleFitByPratt(B');
%         if (isnan(Par(1)) || isinf(Par(1)))
%             Struct(t).Bdat(b).radius = inf;
%             Struct(t).Bdat(b).rBar = inf*ones(1,2);
%         else
%             Struct(t).Bdat(b).radius = Par(3);
%             Struct(t).Bdat(b).rBar = Par(1:2)';
%         end

        [ y, Ecircle ] = seg.fitToCircle( B', nB, x0, D );
        
        delta = bsxfun(@minus,B,x0);
        lineDistance = mean( (delta(1,:) .* nB(1) + delta(2,:) .* nB(2)).^2 );
        if (Ecircle < lineDistance && length(bx) > 3)
            Struct(t).Bdat(b).radius = sqrt( y^2 + (D/2)^2 );
            Struct(t).Bdat(b).rBar = x0 + y * nB;
            Struct(t).Bdat(b).fitEnergy = Ecircle;
        else
            Struct(t).Bdat(b).radius = inf;
            Struct(t).Bdat(b).rBar = inf*ones(size(x0));
            Struct(t).Bdat(b).fitEnergy = lineDistance;
        end
%         plot.skel(Struct,'k',0)
%         scatter(rv(1),rv(2))
%         hold all
%         scatter(rnv(1),rnv(2))
%         scatter(Struct(t).Bdat(b).rBar(1),Struct(t).Bdat(b).rBar(2))
%         pause
    end

end

end

function [ F ] = fit_circle( yc, xc, nB, xdata, ydata )
    Rc = xc + yc * nB;
    F = mean( abs(xdata.^2 - 2*Rc(1)*xdata + ydata.^2 - 2*Rc(2)*ydata) );
end

function [ dilPix ] = dilatePix(pix, d, dim)
    nPix = length(pix);
    dilPix = zeros(1,(2*d+1)^2*nPix);
    k = 1;
    for x = -d:d
        for y = -(d-abs(x)):(d-abs(x))
            dilPix(nPix*(k-1)+1:k*nPix) = pix + y + dim(1)*x;
            k = k + 1;
        end
    end
    dilPix = unique(dilPix(find(dilPix)));
end

