function [ stress, xVec, yVec ] = takeContinuumLimit( Struct, smScale, scale, PN, xG, yG, ERes )
    % TAKE CONTINUUM LIMIT 
    
    xVec = linspace(1,1738,218);
%     xVec = linspace(100,1600,218);
    yVec = linspace(1,2050,257);
%     yVec = linspace(100,1950,257);

    [X,Y] = meshgrid(xVec,yVec);

    stress = cell(size(PN));
    
    for t = 1:length(PN)
        stress{t} = zeros(257,218,3);

        TB = cell(length(Struct(t).Bdat),1);
        rBB = cell(length(Struct(t).Bdat),1);
        RB = cell(length(Struct(t).Bdat),1);
        DB = cell(length(Struct(t).Bdat),1);
        
        bCells = zeros(length(Struct(t).Bdat),2);
        for b = 1:length(Struct(t).Bdat)
            if (length(Struct(t).Bdat(b).cells)==2)
                bCells(b,:) = Struct(t).Bdat(b).cells;
            end
        end
        bCells = sort(bCells,2);
    
        for ii = 1:size(yG{t},2)
            for jj = 1:size(xG{t},2)
                if (ERes{t}(ii,jj) < 10) 
                    
                    [ bInd ] = pressure.matchBonds(PN{t}{ii,jj},bCells);
                    rEst = PN{t}{ii,jj}.computePrimalVerts();
                    [ ~, ~, ~, bulkVerts ] = fitDual.returnSubGraph( Struct(t), xG{t}(:,jj), yG{t}(:,ii) );
                    [ extVerts ] = fitDual.returnExtVerts( Struct(t), bulkVerts );
                    [ r0 ] = generate.embedSubGraph( Struct(t), bulkVerts, extVerts, Struct(t).embed  );
                    [ ~, dV ] = fitDual.subATN.computeSubDiffOperators( Struct(t), xG{t}(:,jj), yG{t}(:,ii) );

                    rV = zeros(size(r0));
                    rV(1:size(rEst,1),:) = rEst;
                    rV((size(rEst,1)+1):size(rV,1),:) = r0((size(rEst,1)+1):size(rV,1),:);
                    D = sqrt(sum((dV*rV).^2,2));

                    rB = PN{t}{ii,jj}.d0*PN{t}{ii,jj}.q; 
                    rB = rB*[0,-1;1,0];
                    rB = bsxfun(@rdivide,rB,sqrt(sum(rB.^2,2)));

                    R = .5*abs(PN{t}{ii,jj}.d1')*rEst;
                    badBonds = find(sum(abs(PN{t}{ii,jj}.d1),1) == 1);
                    extVerts = zeros(size(badBonds));
                    for nn = 1:length(badBonds)
                        extVerts(nn) = find(PN{t}{ii,jj}.d1(:,badBonds(nn)));
                    end
                    S = -PN{t}{ii,jj}.d1(extVerts + size(PN{t}{ii,jj}.d1,1)*(badBonds-1));
                    R(badBonds,:) = rEst(extVerts,:) + bsxfun(@times,S'.*D(badBonds)/2,rB(badBonds,:));

                    T = PN{t}{ii,jj}.returnTension();
                    T = scale{t}(ii,jj)*T;
                    
                    for b = 1:length(bInd)
                        if (bInd(b) > 0)
                           TB{bInd(b)} = [TB{bInd(b)},T(b)];
                           DB{bInd(b)} = [DB{bInd(b)},D(b)];
                           rBB{bInd(b)} = [rBB{bInd(b)};rB(b,:)];
                           RB{bInd(b)} = [RB{bInd(b)};R(b,:)];
                        end
                    end
                end
            end
        end
       
        for b = 1:length(TB)
            nInfers = length(TB{b});
            if (nInfers > 0)
                for n = 1:nInfers
                    T = TB{b}(n)/nInfers;
                    D = DB{b}(n);
                    rB = rBB{b}(n,:);
                    R = RB{b}(n,:);
                    nB = rB*[0,-1;1,0];
                    
                    y = bsxfun(@minus,nB(:,1)*X,nB(:,1).*R(:,1)) + bsxfun(@minus,nB(:,2)*Y,nB(:,2).*R(:,2));
                    x = bsxfun(@minus,rB(:,1)*X,rB(:,1).*R(:,1)) + bsxfun(@minus,rB(:,2)*Y,rB(:,2).*R(:,2));
                    
                    I = (exp(-y.^2/(2*smScale^2))/(2*sqrt(2*pi)*smScale)) .* ...
                        (erf( bsxfun(@minus,D/2,x) / sqrt(2*smScale^2) ) - erf( bsxfun(@minus,-D/2,x) / sqrt(2*smScale^2) ) );
                   
                    sigma1 = T.*rB(:,1).*rB(:,1)*I;
                    sigma2 = T.*rB(:,1).*rB(:,2)*I;
                    sigma3 = T.*rB(:,2).*rB(:,2)*I;

                    stress{t}(:,:,1) = stress{t}(:,:,1) + sigma1;
                    stress{t}(:,:,2) = stress{t}(:,:,2) + sigma2;
                    stress{t}(:,:,3) = stress{t}(:,:,3) + sigma3;
                end
            end
        end
%         stress{t} = generate.symmetrizedStress(stress{t});
        
    end
end

