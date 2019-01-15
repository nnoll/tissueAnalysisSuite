function [] = compareTensor( A, B, dX, dY, Lx, Ly )
%COMPARETENSOR 

    AS = A(1:dY:size(A,1),1:dX:size(A,2),:);
    BS = B(1:dY:size(B,1),1:dX:size(B,2),:);

    yG = 1:dY:size(A,1);
    xG = 1:dX:size(A,2);
    
    if (nargin > 4)
        yG = 1.0*(yG * Ly) / size(A,1);
        xG = 1.0*(xG * Lx) / size(A,2);
    end
    
    [XG,YG] = meshgrid(xG,yG);
    lambdaA1 = zeros(size(AS,1),size(AS,2));
    lambdaA2 = zeros(size(AS,1),size(AS,2));
    
    lambdaB1 = zeros(size(BS,1),size(BS,2));
    lambdaB2 = zeros(size(BS,1),size(BS,2));

    A1 = zeros(size(AS,1),size(AS,2),2);
    B1 = zeros(size(BS,1),size(BS,2),2);
    
    badPts = zeros(size(AS,1),size(AS,2));
    for ii = 1:size(AS,1)
        for jj = 1:size(AS,2)
            if (~isnan(AS(ii,jj,1)) && ~isnan(AS(ii,jj,2)) && ~isnan(AS(ii,jj,3)) && ...
                ~isnan(BS(ii,jj,1)) && ~isnan(BS(ii,jj,2)) && ~isnan(BS(ii,jj,3)) )
            
                sigmaA = [AS(ii,jj,1),AS(ii,jj,2);AS(ii,jj,2),AS(ii,jj,3)];
                [VA,DA] = eig(sigmaA);
                
                sigmaB = [BS(ii,jj,1),BS(ii,jj,2);BS(ii,jj,2),BS(ii,jj,3)];
                [VB,DB] = eig(sigmaB);
                
                lambdaA1(ii,jj) = DA(1,1);
                lambdaA2(ii,jj) = DA(2,2);
                lambdaB1(ii,jj) = DB(1,1);
                lambdaB2(ii,jj) = DB(2,2);

                A1(ii,jj,:) = VA(:,2);
                B1(ii,jj,:) = VB(:,2);
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
    
    area_B = abs(lambdaB1.*lambdaB2);
    scaleB = R0 / mean(area_B(badPts==0));
    lambdaB1 = sqrt(scaleB)*lambdaB1;
    lambdaB2 = sqrt(scaleB)*lambdaB2;
    lambdaB1(badPts==1) = 0;
    lambdaB2(badPts==1) = 0;
    lambdaB1 = lambdaB1(:);
    lambdaB2 = lambdaB2(:);

    aM1(:,1) = XG(:) + bsxfun(@times,lambdaA2,reshape(A1(:,:,1),size(lambdaA1)));
    aM1(:,2) = YG(:) + bsxfun(@times,lambdaA2,reshape(A1(:,:,2),size(lambdaA2)));
    
    aM2(:,1) = XG(:) - bsxfun(@times,lambdaA2,reshape(A1(:,:,1),size(lambdaA1)));
    aM2(:,2) = YG(:) - bsxfun(@times,lambdaA2,reshape(A1(:,:,2),size(lambdaA1)));
    
    bM1(:,1) = XG(:) + bsxfun(@times,lambdaB2,reshape(B1(:,:,1),size(lambdaA1)));
    bM1(:,2) = YG(:) + bsxfun(@times,lambdaB2,reshape(B1(:,:,2),size(lambdaA1)));
    
    bM2(:,1) = XG(:) - bsxfun(@times,lambdaB2,reshape(B1(:,:,1),size(lambdaA1)));
    bM2(:,2) = YG(:) - bsxfun(@times,lambdaB2,reshape(B1(:,:,2),size(lambdaA1)));
    
    t = linspace(0,2*pi,50);
    
    XA = .75*lambdaA2 * cos(t);
    YA = .75*lambdaA1 * sin(t);
    
    XB = lambdaB2 * cos(t);
    YB = lambdaB1 * sin(t);
    
    wA = atan2( aM2(:,2) - aM1(:,2), aM2(:,1) - aM1(:,1) );
    wB = atan2( bM2(:,2) - bM1(:,2), bM2(:,1) - bM1(:,1) );

    xA = bsxfun(@plus,XG(:),bsxfun(@times,XA,cos(wA)) - bsxfun(@times,YA,sin(wA)));
    yA = bsxfun(@plus,YG(:),bsxfun(@times,XA,sin(wA)) + bsxfun(@times,YA,cos(wA)));
    
    xB = bsxfun(@plus,XG(:),bsxfun(@times,XB,cos(wB)) - bsxfun(@times,YB,sin(wB)));
    yB = bsxfun(@plus,YG(:),bsxfun(@times,XB,sin(wB)) + bsxfun(@times,YB,cos(wB)));
    
    hold on
    plot(xA',yA','Color',[0,0,0],'LineWidth',2)
    plot(xB',yB','Color',[0,.5,1],'LineWidth',2)

    if (nargin < 4)
        set(gca,'XLim',[1,size(A,2)],'YLim',[1,size(A,1)])
    else
        set(gca,'XLim',[1,Lx],'YLim',[1,Ly])
    end
    
    axis equal
    
end

