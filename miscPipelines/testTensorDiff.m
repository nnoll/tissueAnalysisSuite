A = expMyo{50};
% theta = linspace(0,pi,100);
theta = pi/2;

for ii = 1:length(theta)
    c = cos(theta(ii));
    s = sin(theta(ii));

    B = zeros(size(A));

    B(:,:,1) = c^2*A(:,:,1) - 2*c*s*A(:,:,2) + s^2*A(:,:,3);
    B(:,:,2) = s*c*(A(:,:,1)-A(:,:,3)) + A(:,:,2)*(c^2-s^2);
    B(:,:,3) = c^2*A(:,:,3) - 2*c*s*A(:,:,2) + s^2*A(:,:,1);

    [ D ] = continuum.tensorDiff( A, B );
    ovlap(ii) = sqrt(nanmean(D(:).^2));
end