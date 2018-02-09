N = 11;
alpha = .03;
scale = 200;
mode = 1;

betaVec = linspace(0,1,8);
% betaVec = linspace(.1,1,8);
% betaVec = .1;
% gammaVec = linspace(0,.05,8);
gammaVec = linspace(0,.1,8);
niter = 50;

cOld = zeros(length(betaVec),length(gammaVec),niter);
cNew = zeros(length(betaVec),length(gammaVec),niter);

for n = 1:niter
    n
    for ii = 1:length(betaVec)
        for jj = 1:length(gammaVec)
            if (betaVec(ii) == 0)
                mode = 1;
            else
                mode = 3;
            end
            [ cOld(ii,jj,n), cNew(ii,jj,n) ] = testInverseMethods( N, alpha, betaVec(ii), gammaVec(jj), scale, mode );
        end
    end
end
