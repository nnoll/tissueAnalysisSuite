function [  ] = plotStress( A, dX, dY )
    % PLOT STRESS 

    AS = A(1:dY:size(A,1),1:dX:size(A,2),:);
    yG = 1:dY:size(A,1);
    xG = 1:dX:size(A,2);
    
    if (nargin > 4)
        yG = 1.0*(yG * Ly) / size(A,1);
        xG = 1.0*(xG * Lx) / size(A,2);
    end
    
    [XG,YG] = meshgrid(xG,yG);
    lambdaA1 = zeros(size(AS,1),size(AS,2));
    lambdaA2 = zeros(size(AS,1),size(AS,2));
    
    A1 = zeros(size(AS,1),size(AS,2),2);
    
    badPts = zeros(size(AS,1),size(AS,2));
    for ii = 1:size(AS,1)
        for jj = 1:size(AS,2)
            if (~isnan(AS(ii,jj,1)) && ~isnan(AS(ii,jj,2)) && ~isnan(AS(ii,jj,3)) )
            
                sigmaA = [AS(ii,jj,1),AS(ii,jj,2);AS(ii,jj,2),AS(ii,jj,3)];
                [VA,DA] = eig(sigmaA);

                lambdaA1(ii,jj) = DA(1,1);
                lambdaA2(ii,jj) = DA(2,2);

                A1(ii,jj,:) = VA(:,2);
            else
                badPts(ii,jj) = 1;
            end
        end
    end
    
    if (nargin > 4)
        R0 = .04*(dX^2 + dY^2)*(Lx*Ly/(size(A,1)*size(A,2)));
    else
        R0 = .04*(dX^2 + dY^2);
    end
    
    area_A = abs(lambdaA1.*lambdaA2);
    scaleA = R0 / mean(area_A(badPts==0));
    lambdaA1 = sqrt(scaleA)*lambdaA1;
    lambdaA2 = sqrt(scaleA)*lambdaA2;
    lambdaA1(badPts==1) = 0;
    lambdaA2(badPts==1) = 0;
    lambdaA1 = lambdaA1(:);
    lambdaA2 = lambdaA2(:);

    aM1(:,1) = XG(:) + bsxfun(@times,lambdaA2,reshape(A1(:,:,1),size(lambdaA1)));
    aM1(:,2) = YG(:) + bsxfun(@times,lambdaA2,reshape(A1(:,:,2),size(lambdaA2)));
    
    aM2(:,1) = XG(:) - bsxfun(@times,lambdaA2,reshape(A1(:,:,1),size(lambdaA1)));
    aM2(:,2) = YG(:) - bsxfun(@times,lambdaA2,reshape(A1(:,:,2),size(lambdaA1)));
    
    t = linspace(0,2*pi,50);
    
    XA = lambdaA2 * cos(t);
    YA = lambdaA1 * sin(t);

    wA = atan2( aM2(:,2) - aM1(:,2), aM2(:,1) - aM1(:,1) );

    xA = bsxfun(@plus,XG(:),bsxfun(@times,XA,cos(wA)) - bsxfun(@times,YA,sin(wA)));
    yA = bsxfun(@plus,YG(:),bsxfun(@times,XA,sin(wA)) + bsxfun(@times,YA,cos(wA)));
    
    hold on
    plot(xA',yA','Color',[1,0,0],'LineWidth',2)

    if (nargin < 4)
        set(gca,'XLim',[1,size(A,2)],'YLim',[1,size(A,1)])
    else
        set(gca,'XLim',[1,Lx],'YLim',[1,Ly])
    end
    
    axis equal
    
end

