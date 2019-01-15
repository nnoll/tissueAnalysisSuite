function [ Q, rV, theta, tri ] = fitThetaModes( Q, p, Struct )
    % FIT THETA MODES 

    [ tri, ~, ~, ~, ~, rV ] = fitATN.returnGraph( Struct );
    [ tri ] = fitAFN.returnPTri(Q,p,tri);
    
    p1 = p(tri(:,1));
    p2 = p(tri(:,2));
    p3 = p(tri(:,3));
    
    q1 = Q(tri(:,1),:);
    q2 = Q(tri(:,2),:);
    q3 = Q(tri(:,3),:);
    
    pq1 = bsxfun(@times,p1,q1);
    pq2 = bsxfun(@times,p2,q2);
    pq3 = bsxfun(@times,p3,q3);
    
    t1 = pq3 - pq2;
    t2 = pq1 - pq3;
    t3 = pq2 - pq1;
    
    eta = bsxfun(@times,p1,t1) + bsxfun(@times,p2,t2) + bsxfun(@times,p3,t3);
    delta = dot(eta,rV,2) + .5*(p3-p2).*p1.*sum(q1.^2,2) + .5.*(p1-p3).*p2.*sum(q2.^2,2) + .5.*(p2-p1).*p3.*sum(q3.^2,2);
    
    L = zeros(size(tri,1),size(Q,1));
    L((1:size(tri,1)) + size(tri,1)*(tri(:,1)'-1)) = .5*(p2-p3).*p1;
    L((1:size(tri,1)) + size(tri,1)*(tri(:,2)'-1)) = .5*(p3-p1).*p2;
    L((1:size(tri,1)) + size(tri,1)*(tri(:,3)'-1)) = .5*(p1-p2).*p3;

    theta = pinv(L) * delta;

end

