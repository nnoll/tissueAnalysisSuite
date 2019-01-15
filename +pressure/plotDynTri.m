function [ mov ] = plotDynTri( PN )
%PLOTDYNTRI Summary of this function goes here
%   Detailed explanation goes h

    for t = 1:length(PN)
        clf
        PN{t}.plotTri();
        axis tight
        XLim = get(gca,'XLim');
        YLim = get(gca,'YLim');

        if (t == 1)
            Xr = XLim;
            Yr = YLim;
        else
            if (XLim(1) < Xr(1))
                Xr(1) = XLim(1);
            end
            
            if (XLim(2) > Xr(2))
                Xr(2) = XLim(2);
            end
            
            if (YLim(1) < Yr(1))
                Yr(1) = YLim(1);
            end
            
            if (YLim(2) > Yr(2))
                Yr(2) = YLim(2);
            end
        end
            
    end

    for t = 1:length(PN)
        clf
        PN{t}.plotTri();
        set(gca,'XLim',Xr,'YLim',[Yr(1),500])
        mov(t) = getframe();
    end

end

