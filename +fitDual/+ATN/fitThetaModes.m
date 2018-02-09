function [ Q, theta, tri ] = fitThetaModes( Q, Struct, extCell )
    % FIT THETA MODES 

    Rot = [0,-1;1,0];
    [ tri, ~, ~, ~, ~, rV ] = fitDual.returnGraph( Struct, extCell );
    
    [ r0 ] = fitDual.ATN.returnVertexPositions( Q, zeros(size(Q,1),1), tri );
    delta = rV - r0;
    
    [ tri ] = fitDual.orderTri( Q, tri );
    
    q1 = Q(tri(:,1),:);
    q2 = Q(tri(:,2),:);
    q3 = Q(tri(:,3),:);
    
    t1 = q3 - q2;
    t2 = q1 - q3;
    t3 = q2 - q1;
    
    S0 = .5*abs(t1(:,1).*t2(:,2) - t1(:,2).*t2(:,1));
    L = zeros(2*size(tri,1),size(Q,1));
    
    t1 = t1*Rot';
    t2 = t2*Rot';
    t3 = t3*Rot';
    
    nT = size(tri,1);
    
    for t = 1:size(tri,1)
        
        % X Component
        L(t,tri(t,1)) = t1(t,1)/(4*S0(t));
        L(t,tri(t,2)) = t2(t,1)/(4*S0(t));
        L(t,tri(t,3)) = t3(t,1)/(4*S0(t));
        
        % Y Component
        L(t+nT,tri(t,1)) = t1(t,2)/(4*S0(t));
        L(t+nT,tri(t,2)) = t2(t,2)/(4*S0(t));
        L(t+nT,tri(t,3)) = t3(t,2)/(4*S0(t));
        
    end

    L = sparse([L;ones(1,size(L,2))]);
    theta = L \ [delta(:);0];
    deltaR = L*theta;
    deltaR = deltaR(1:end-1);
    deltaR = reshape(deltaR,length(deltaR)/2,2);
    dMag = sqrt(dot(deltaR,deltaR,2))/mean(sqrt(sum(delta.^2,2)));
    badVerts = find(dMag > 15);
    if (~isempty(badVerts))
        L([badVerts,badVerts+nT],:) = [];
        b = delta; b(badVerts,:) = [];
        b = [b(:);0];
        theta = pinv(full(L)) * b;
        deltaR = L*theta;
        deltaR = deltaR(1:end-1);
        deltaR = reshape(deltaR,length(deltaR)/2,2);
        dMag = sqrt(dot(deltaR,deltaR,2))/mean(sqrt(sum(delta.^2,2)));
    end
    
end

