function [ stress, xVec, yVec ] = takeContinuumLimit( Struct, smScale, scale, PN, xG, yG, ERes )
    % TAKE CONTINUUM LIMIT 
    
%     xVec = linspace(1,1738,218);
    xVec = linspace(740,1000,218);
%     yVec = linspace(1,2050,257);
    yVec = linspace(1020,1180,257);

    [X,Y] = meshgrid(xVec,yVec);
    X = X(:);   
    Y = Y(:);

    stress = cell(size(PN));
    for t = 1:length(PN)
        stress{t} = zeros(257,218,3);

        for ii = 14 %1:size(yG{t},2)
            for jj = [12,13,14] %1:size(xG{t},2)
                if (ERes{t}(ii,jj) < 10) 
                    
                    rEst = PN{t}{ii,jj}.computePrimalVerts();
                    [ ~, ~, ~, bulkVerts ] = fitDual.returnSubGraph( Struct(t), xG{t}(:,jj), yG{t}(:,ii) );
                    [ extVerts ] = fitDual.returnExtVerts( Struct(t), bulkVerts );
                    [ r0 ] = generate.embedSubGraph( Struct(t), bulkVerts, extVerts, Struct(t).embed  );
                    [ ~, dV ] = fitDual.subATN.computeSubDiffOperators( Struct(t), xG{t}(:,jj), yG{t}(:,ii) );
                    
%                     dV = PN{t}{ii,jj}.d1; 
                    goodBonds = ((sum(abs(PN{t}{ii,jj}.d1),1)) == 2);
%                     dV = sparse(dV);
                    
                    rV = zeros(size(r0));
                    rV(1:size(rEst,1),:) = rEst;
                    rV((size(rEst,1)+1):size(rV,1),:) = r0((size(rEst,1)+1):size(rV,1),:);
                    D = sqrt(sum((dV*rV).^2,2));

                    rB = PN{t}{ii,jj}.d0*PN{t}{ii,jj}.q; 
                    rB = rB*[0,-1;1,0];
                    rB = bsxfun(@rdivide,rB,sqrt(sum(rB.^2,2)));
                    nB = rB * [0,-1;1,0];

                    R = .5*abs(PN{t}{ii,jj}.d1')*rEst;
                    badBonds = find(sum(abs(PN{t}{ii,jj}.d1),1) == 1);
                    extVerts = zeros(size(badBonds));
                    for nn = 1:length(badBonds)
                        extVerts(nn) = find(PN{t}{ii,jj}.d1(:,badBonds(nn)));
                    end
                    S = -PN{t}{ii,jj}.d1(extVerts + size(PN{t}{ii,jj}.d1,1)*(badBonds-1));
                    R(badBonds,:) = rEst(extVerts,:) + bsxfun(@times,S'.*D(badBonds)/2,rB(badBonds,:));
%                     scatter(R(:,1),R(:,2),'r','filled')
%                     hold on
%                     scatter(R(badBonds,1),R(badBonds,2),'go')
%                     PN{t}{ii,jj}.plotPrimal();
%                     pause
%                     
                    T = PN{t}{ii,jj}.returnTension();
                    T = scale{t}(ii,jj)*T;
%                     T(~goodBonds) = .5*T(~goodBonds);
%                     T = T(goodBonds);
%                     R = R(goodBonds,:);
%                     D = D(goodBonds);
%                     rB = rB(goodBonds,:);
%                     nB = nB(goodBonds,:);
%                     
                    y = bsxfun(@minus,nB(:,1)*X',nB(:,1).*R(:,1)) + bsxfun(@minus,nB(:,2)*Y',nB(:,2).*R(:,2));
                    x = bsxfun(@minus,rB(:,1)*X',rB(:,1).*R(:,1)) + bsxfun(@minus,rB(:,2)*Y',rB(:,2).*R(:,2));

                    I = (exp(-y.^2/(2*smScale^2))/(2*sqrt(2*pi)*smScale)) .* ...
                        (erf( bsxfun(@minus,D/2,x) / sqrt(2*smScale^2) ) - erf( bsxfun(@minus,-D/2,x) / sqrt(2*smScale^2) ) );

%                     I( :, (Y<=yG{t}(1,ii) | Y>=yG{t}(2,ii)) | (X<=xG{t}(1,jj) | X>=xG{t}(2,jj)) ) = 0;
                    sigma1 = sum(bsxfun(@times,T.*rB(:,1).*rB(:,1),I),1)';
                    sigma2 = sum(bsxfun(@times,T.*rB(:,1).*rB(:,2),I),1)';
                    sigma3 = sum(bsxfun(@times,T.*rB(:,2).*rB(:,2),I),1)';

                    sigma1 = reshape(sigma1,257,218);
                    sigma2 = reshape(sigma2,257,218);
                    sigma3 = reshape(sigma3,257,218);

                    stress{t}(:,:,1) = stress{t}(:,:,1) + sigma1;
                    stress{t}(:,:,2) = stress{t}(:,:,2) + sigma2;
                    stress{t}(:,:,3) = stress{t}(:,:,3) + sigma3;

                end
            end
        end
        
%         stress{t} = generate.symmetrizedStress(stress{t});
        
    end
end

