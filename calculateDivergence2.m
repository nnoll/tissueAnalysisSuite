Tdiv = cell(length(PN),1);
Mdiv = cell(length(PN),1);
Sdiv = cell(length(PN),1);

for t = timePts
    t
    bCells = vertcat(tStruct(t).Bdat.cells);
    bCells = sort(bCells,2);
    
    Tdiv{t} = zeros(size(PN{t}));
    Mdiv{t} = zeros(size(PN{t}));
    Sdiv{t} = zeros(size(PN{t}));
    
    for ii = 1:size(yG{t},2)
        for jj = 1:size(xG{t},2)
            if (ERes{t}(ii,jj) < 10)
                
                [ ~, ~, ~, bulkVerts ] = fitDual.returnSubGraph( tStruct(t), xG{t}(:,jj), yG{t}(:,ii) );
                [ extVerts ] = fitDual.returnExtVerts( tStruct(t), bulkVerts );

                [ r0 ] = generate.embedSubGraph( tStruct(t), bulkVerts, extVerts, tStruct(t).embed  );

                rEst = PN{t}{ii,jj}.computePrimalVerts();
                dV = PN{t}{ii,jj}.d1;              
                dV = sparse(dV);
                
                rB = dV'*rEst;
                rB = bsxfun(@rdivide,rB,sqrt(sum(rB.^2,2)));
                
                Tvec = PN{t}{ii,jj}.d0*PN{t}{ii,jj}.q;
                Tvec = Tvec * [0,-1;1,0];
                Tvec = scale{t}(ii,jj)*Tvec;
                
                bInd = pressure.matchBonds(PN{t}{ii,jj},bCells);
                M = zeros(size(bInd));
                M(bInd>0) = [tStruct(t).Bdat(bInd(bInd>0)).myo];
 
                S = zeros(size(bInd));
                S(bInd>0) = [tStruct(t).Bdat(bInd(bInd>0)).stress];
                
                Mvec = bsxfun(@times,M',rB);
                Svec = bsxfun(@times,S',rB);

                Tnet = dV*Tvec;
                Mnet = dV*Mvec;
                Snet = dV*Svec;

%                 clf
%                 quiver(rEst(:,1),rEst(:,2),Mnet(:,1),Mnet(:,2))
%                 hold on
%                 quiver(rEst(:,1),rEst(:,2),Snet(:,1),Snet(:,2))
%                 pause
%                 
                TMag = sqrt(sum(Tnet.^2,2));
                MMag = sqrt(sum(Mnet.^2,2));
                SMag = sqrt(sum(Snet.^2,2));
                
                Tdiv{t}(ii,jj) = nanmedian(TMag);
                Mdiv{t}(ii,jj) = nanmean(MMag);
                Sdiv{t}(ii,jj) = nanmedian(SMag);
            end
        end
    end
end