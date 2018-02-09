clear all
N = 11;
alpha = .03;
scale = 100;

% betaVec = linspace(0,1,10);
betaVec = 1;
niter = 1;
for n = 1:niter
    for ii = 1:length(betaVec)
        [ chi(ii,n), chiC(ii,n), chiCV(ii,n), kappa, L, Struct ] = returnCompatibility( N, alpha, betaVec(ii), scale );
    end
end

