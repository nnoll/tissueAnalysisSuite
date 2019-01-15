function [ tri, L, Struct, q, p, theta, contin ] = calculateTri( q, p, theta, b0, alpha )
    % CALCULATE TRI 
    
    if (nargin == 4)
        alpha = 200;
    end
    
    qT = alpha*q;
    thetaT = alpha^2*theta;
    
    XMin = min(qT(:,1)); XMax = max(qT(:,1));
    YMin = min(qT(:,2)); YMax = max(qT(:,2));

    dimX = round(XMax - XMin) + 50;
    dimY = round(YMax - YMin) + 50;

    qT = bsxfun(@plus,qT,[XMax+25,YMax+25]);

    [ L, Struct, indx, D ] = pressure.generateWeightedVoronoi( qT, p, thetaT, [dimY,dimX] );
    
    q = q(indx,:);
    p = p(indx);
    theta = theta(indx);

    if (all(ismember(find(b0),find(indx))))
        contin = 0;
    else
        contin = 1;
    end
    
    tri = zeros(length(Struct.Vdat),3);
    if (contin == 0)
        for t = 1:size(tri,1)
            nCells = Struct.Vdat(t).ncells; 
            if (length(nCells) ~= 3)
                contin = 1;
                t
                break
            end
            Rq = q(nCells,:);
            delta = bsxfun(@minus,Rq,mean(Rq,1));
            angle = mod(atan2(delta(:,2),delta(:,1)),2*pi);
            [~,ind] = sort(angle);
            tri(t,:) = nCells(ind);
        end

        % Find cells not indexed in triangulation and remove.
        iCells = unique(tri(:));
        delCells = 1:size(q,1);
        delCells = delCells(~ismember(delCells,iCells));

        changeIndx = zeros(size(tri));
        for ii = 1:length(delCells)
           changeIndx = changeIndx - (tri > delCells(ii));
        end
        tri = tri + changeIndx;

        q(delCells,:) = [];
        p(delCells) = [];
        theta(delCells) = [];
    end
    
%     figure(3)
%     PN = pressure.net(bsxfun(@plus,alpha*q,[XMax+25,YMax+25]),alpha^2*theta,p,tri);
%     imshow(D,[])
%     hold all
%     PN.plotCurvedPrimal()
%     pause
    
end

