function [ mov ] = plotVector( v, x, y )
%PLOTVECTOR 

    [x,y] = meshgrid(x,y);
    for t = 1:length(v)
        quiver(x(:),y(:),v{t}(1,:)',v{t}(2,:)',0,'b','LineWidth',2);
        set(gca,'XLim',[1,1738],'YLim',[1,2050])
        mov(t) = getframe();
    end
end

