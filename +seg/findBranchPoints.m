function [ branchPoints, bxy ] = findBranchPoints( skel )

% Original version by Nick Noll
% branchPoints = zeros(size(skel));
% rows = 2:(size(skel,1)-1);
% cols = 2:(size(skel,2)-1);
% branchPoints(rows,cols) = skel(rows+1,cols) + skel(rows-1,cols) + ...
%                           skel(rows,cols+1) + skel(rows,cols-1);
% branchPoints = branchPoints .* skel;
% branchPoints = branchPoints >= 3;

% New version (NPM / Matlab builtin)
branchPoints = bwmorph(skel, 'branchpoints') ;

if nargout > 1
    [y,x] = find(branchPoints) ;
    bxy = [x(:), y(:)] ;
end

end

