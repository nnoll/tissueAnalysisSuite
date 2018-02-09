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
    A2 = zeros(size(AS,1),size(AS,2),2);
    
    B1 = zeros(size(BS,1),size(BS,2),2);
    B2 = zeros(size(BS,1),size(BS,2),2);
    
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

                A1(ii,jj,:) = VA(:,1);
                A2(ii,jj,:) = VA(:,2);
                B1(ii,jj,:) = VB(:,1);
                B2(ii,jj,:) = VB(:,2);
            else
                badPts(ii,jj) = 1;
            end
        end
    end
    
    if (nargin > 4)
        R0 = .05*(dX^2 + dY^2)*(Lx*Ly/(size(A,1)*size(A,2)));
    else
        R0 = .05*(dX^2 + dY^2);
    end
    
    area_A = abs(lambdaA1.*lambdaA2);
    scaleA = R0 / mean(area_A(badPts==0));
    lambdaA1 = sqrt(scaleA)*lambdaA1;
    lambdaA2 = sqrt(scaleA)*lambdaA2;
    lambdaA1(badPts==1) = 0;
    lambdaA2(badPts==1) = 0;

    area_B = abs(lambdaB1.*lambdaB2);
    scaleB = R0 / mean(area_B(badPts==0));
    lambdaB1 = sqrt(scaleB)*lambdaB1;
    lambdaB2 = sqrt(scaleB)*lambdaB2;
    lambdaB1(badPts==1) = 0;
    lambdaB2(badPts==1) = 0;
    
    Ax1 = XG + lambdaA1 .* A1(:,:,1) + lambdaA2 .* A2(:,:,1);
    Ax2 = XG - lambdaA1 .* A1(:,:,1) + lambdaA2 .* A2(:,:,1);
    Ax3 = XG - lambdaA1 .* A1(:,:,1) - lambdaA2 .* A2(:,:,1);
    Ax4 = XG + lambdaA1 .* A1(:,:,1) - lambdaA2 .* A2(:,:,1);
    
    Ay1 = YG + lambdaA1 .* A1(:,:,2) + lambdaA2 .* A2(:,:,2);
    Ay2 = YG - lambdaA1 .* A1(:,:,2) + lambdaA2 .* A2(:,:,2);
    Ay3 = YG - lambdaA1 .* A1(:,:,2) - lambdaA2 .* A2(:,:,2);
    Ay4 = YG + lambdaA1 .* A1(:,:,2) - lambdaA2 .* A2(:,:,2);

    Bx1 = XG + lambdaB1 .* B1(:,:,1) + lambdaB2 .* B2(:,:,1);
    Bx2 = XG - lambdaB1 .* B1(:,:,1) + lambdaB2 .* B2(:,:,1);
    Bx3 = XG - lambdaB1 .* B1(:,:,1) - lambdaB2 .* B2(:,:,1);
    Bx4 = XG + lambdaB1 .* B1(:,:,1) - lambdaB2 .* B2(:,:,1);
    
    By1 = YG + lambdaB1 .* B1(:,:,2) + lambdaB2 .* B2(:,:,2);
    By2 = YG - lambdaB1 .* B1(:,:,2) + lambdaB2 .* B2(:,:,2);
    By3 = YG - lambdaB1 .* B1(:,:,2) - lambdaB2 .* B2(:,:,2);
    By4 = YG + lambdaB1 .* B1(:,:,2) - lambdaB2 .* B2(:,:,2);
    
    hold on
    plot([Ax1(:),Ax2(:),Ax3(:),Ax4(:),Ax1(:)]',[Ay1(:),Ay2(:),Ay3(:),Ay4(:),Ay1(:)]','LineWidth',2,'Color',[255,0,0]/255)
    plot([Bx1(:),Bx2(:),Bx3(:),Bx4(:),Bx1(:)]',[By1(:),By2(:),By3(:),By4(:),By1(:)]','LineWidth',2,'Color',[0,0,255]/255)

    if (nargin < 4)
        set(gca,'XLim',[1,size(A,2)],'YLim',[1,size(A,1)])
    else
        set(gca,'XLim',[1,Lx],'YLim',[1,Ly])
    end
    
    axis equal
    
end

