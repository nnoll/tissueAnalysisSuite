function [ lambda, delta ] = scaleRelation( q1, q2 )
    % SCALE RELATION 

    A = [q2(:),[ones(size(q1,1),1);zeros(size(q1,1),1)],[zeros(size(q1,1),1);ones(size(q1,1),1)]];
    b = A \ q1(:);
    lambda = b(1);
    delta = b(2:3);
    
end

