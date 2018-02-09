function [ branchPoints ] = findBranchPoints( skel )
    
    branchPoints = zeros(size(skel));
    rows = 2:(size(skel,1)-1);
    cols = 2:(size(skel,2)-1);
    branchPoints(rows,cols) = skel(rows+1,cols) + skel(rows-1,cols) + ...
                              skel(rows,cols+1) + skel(rows,cols-1);
    branchPoints = branchPoints .* skel;
    branchPoints = branchPoints >= 3;

end

