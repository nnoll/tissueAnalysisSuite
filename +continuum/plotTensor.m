function [  ] = plotTensor( Phi, xG, yG )
    % PLOT STRESS 
   
    [XG,YG] = meshgrid(xG,yG);
    radius = zeros(size(Phi,1),size(Phi,2));
    XP = zeros(size(Phi,1),size(Phi,2));
    YP = zeros(size(Phi,1),size(Phi,2));
    aniso = zeros(size(Phi,1),size(Phi,2));
    
    for ii = 1:size(Phi,3)
        for jj = 1:size(Phi,4)
            Sigma = .5*(Phi(:,:,ii,jj) + Phi(:,:,ii,jj)');
            radius(ii,jj) = Phi(1,1,ii,jj) + Phi(2,2,ii,jj);

            [V,D] = eig(Sigma);
            aniso(ii,jj) = abs(D(2,2) - D(1,1));
            XP(ii,jj) = V(1,2);
            YP(ii,jj) = V(2,2);
        end
    end
    
    radius = radius - min(radius(:));
    radius = radius / mean(radius(:));
    dX = xG(2) - xG(1);
    dY = yG(2) - yG(1);

    R0 = .2*sqrt(dX^2 + dY^2);
    radius = R0 * radius;

    RP = sqrt( XP.^2 + YP.^2 );
    XP = XP ./ RP;
    YP = YP ./ RP;
    
    aniso = aniso / mean(aniso(:));
    XP = R0 * aniso .* XP;
    YP = R0 * aniso .* YP;
    
    viscircles([XG(:),YG(:)],radius(:),'LineWidth',2,'DrawBackgroundCircle',0,'EdgeColor',[0,0,1]);
    hold on
    plot([XG(:)+XP(:),XG(:)-XP(:)]',[YG(:)+YP(:),YG(:)-YP(:)]','LineWidth',3,'Color',[1,0,0])
    set(gca,'XLim',[min(xG),max(xG)],'YLim',[min(yG),max(yG)])

end

