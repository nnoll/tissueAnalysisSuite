function [ v, xM, yM, mov ] = PIV( dat, nX, nY )

    % Define a grid over the image coordinates.
    [ xG, yG ] = measure.cellGrid( dat(:,:,1), nX, nY );
    xG = round(xG);
    yG = round(yG);
    
    xG(xG > size(dat,2)) = size(dat,2);
    yG(yG > size(dat,1)) = size(dat,1);
    
    v = cell(size(dat,3)-1,1);
    for t = 2:size(dat,3)
        v{t-1} = zeros(2,size(yG,2),size(xG,2));
        
        % Do correlation analysis for each grid element.
        for ii = 1:size(yG,2)
            for jj = 1:size(xG,2)
                v{t-1}(:,ii,jj) = continuum.computeShift(dat(yG(1,ii):yG(2,ii),xG(1,jj):xG(2,jj),t-1),dat(yG(1,ii):yG(2,ii),xG(1,jj):xG(2,jj),t));
            end
        end
        
    end

    xM = mean(xG,1);
    yM = mean(yG,1);
    [x,y] = meshgrid(xM,yM);
    
    if (nargout > 3)
        for t = 1:length(v)
            clf
            imshow(dat(:,:,t))
            set(gcf,'units','normalized','outerposition',[0 0 1 1])
            hold on
            quiver(x(:),y(:),v{t}(1,:)',v{t}(2,:)','LineWidth',3,'Color','r')
            mov(t) = getframe(gcf);
        end
    end

end

