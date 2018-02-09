function [ L, Struct, match, Dis ] = generateWeightedVoronoi( r, p, theta, dim )
 % GENERATE WEIGHTED VORONOI 

%     r = round(r);
%     
    [X,Y] = meshgrid(1:dim(2),1:dim(1));
    
    X = X(:);
    Y = Y(:);

    helpx = [ones(size(r,1),1), -2*r(:,1), r(:,1).^2 ];
    helpX = [X.^2 , X, ones(length(X),1)];
    
    helpy = [ones(size(r,1),1), -2*r(:,2), r(:,2).^2 ];
    helpY = [Y.^2 , Y, ones(length(Y),1)];

    D = helpx * helpX' + helpy * helpY';
    D = bsxfun(@plus,D,theta);
    D = bsxfun(@times,p,D);
    D = sqrt(abs(min(D,[],1)));
    D = reshape(D,[dim(1),dim(2)]);
 
%     D = inf*ones(dim);
%     for n = 1:size(r,1)
%         bw = zeros(dim);
%         bw(r(n,2),r(n,1)) = 1;
%         Dtmp = sqrt(p(n))*sqrt(bwdist(bw).^2 + theta(n));
%         D = min(D,Dtmp);
%     end
    
    L = watershed(D);
%     imshow(imdilate(L==0,strel('disk',2)),[])
%     pause
    [ L, Struct ] = seg.generate_structs(L,0,0,1);

%     rgb(:,:,1) = (imdilate(L==0,strel('disk',2)));
    [ L, del ] = seg.removeBadCells(Struct,L,1);
%     rgb(:,:,2) = (imdilate(L==0,strel('disk',2)));
%     rgb(:,:,3) = 0;
%     
%     imshow(uint8(255*rgb));
%     pause
    Dis = D;
%     figure(1)
%     imshow(D,[]);
%     hold all
%     scatter(r(:,1),r(:,2),'b','filled')
    
%     figure(2)
%     imshow(plot.imoverlay(mat2gray(D),imdilate(L==0,strel('disk',3)),[0,0,1]));
    
    % Relabel L
    if (del == 1)
        [ L, Struct ] = seg.generate_structs(L,0,0,1);
    end

    S = regionprops(L,'Centroid');
    Rc = vertcat(S(:).Centroid);
    D = pdist2(Rc,r);
    [match] = track.munkres(D);

end

