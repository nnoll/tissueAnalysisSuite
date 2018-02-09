function [ p ] = seedPressure( q, bCells, t1, t2, r1, r2 )
%SEEDPRESSURE 

    L1 = zeros(size(bCells,1),size(q,1));
    L2 = zeros(size(bCells,1),size(q,1));

    nB = size(bCells,1);
    L1( (1:nB)' + nB*(bCells(:,1)-1) ) = dot(q(bCells(:,1),:) - r1,t1,2);
    L1( (1:nB)' + nB*(bCells(:,2)-1) ) = -dot(q(bCells(:,2),:) - r1,t1,2);

    L2( (1:nB)' + nB*(bCells(:,1)-1) ) = dot(q(bCells(:,1),:) - r2,t2,2);
    L2( (1:nB)' + nB*(bCells(:,2)-1) ) = -dot(q(bCells(:,2),:) - r2,t2,2);
    
    scale = mean(sqrt(dot(q,q,2)));
    b = [zeros(2*size(L1,1),1);scale];
    L = sparse([L1;L2;scale*ones(1,size(q,1))/size(q,1)]);

    p = L \ b;
    p = p / mean(p);

end

