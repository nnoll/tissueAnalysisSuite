function [ chi, chiEst ] = compareComptVelocity( Struct, vPair, cPair, extCell )
%COMPARECOMPTPRESSURE Summary of this function goes here
%   Detailed explanation goes here

    for t = 2:length(Struct)
        
        [ ~, bulkCells, ~, intVerts ] = fitDual.returnGraph( Struct(t-1), extCell ); 

        chi{t-1} = zeros(length(bulkCells),1);
        chiEst{t-1} = zeros(length(bulkCells),1);
        n = 1;
        
        clear R0 R
        R0(:,1) = double([Struct(t-1).Vdat.vertxcoord]); 
        R0(:,2) = double([Struct(t-1).Vdat.vertycoord]);
        R(:,1) = double([Struct(t).Vdat.vertxcoord]);
        R(:,2) = double([Struct(t).Vdat.vertycoord]);

        deltaR = zeros(length(intVerts),2);
        [vTrack,cTrack] = track.propagateVertex(vPair,cPair,Struct,[t-1,t]);

        for v = 1:length(vTrack)
            if (vTrack(v) > 0 && ismember(vTrack(v),intVerts))
                ind = (intVerts == vTrack(v)); 
                deltaR(ind,:) = R(v,:) - R0(vTrack(v),:); 
            end
        end
        
        emptyData = find(sum(deltaR.^2,2)==0);
        if (~isempty(emptyData))
            deltaR(emptyData,:) = isogonal.fitMissingDisplacements(Struct,cTrack,[t-1,t],R0(intVerts(emptyData),:));
        end
    
        bVerts = [Struct(t-1).Bdat.verts];
        kappa = zeros(length(Struct(t-1).Bdat),1);
        for b = 1:length(Struct(t-1).Bdat)

            v1 = Struct(t-1).Bdat(b).verts(1);
            v2 = Struct(t-1).Bdat(b).verts(2);

            ind1 = find(intVerts == v1);
            ind2 = find(intVerts == v2);
            if (~isempty(ind1) && ~isempty(ind2))
                deltaB = deltaR(ind1,:) - deltaR(ind2,:);
                
                rb = double([Struct(t-1).Vdat(v1).vertxcoord;Struct(t-1).Vdat(v1).vertycoord]) - ...
                     double([Struct(t-1).Vdat(v2).vertxcoord;Struct(t-1).Vdat(v2).vertycoord]);
                D = sqrt(sum(rb.^2));
                rb = rb / D;
                
                rbRot = [0,-1;1,0] * rb;
                theta = dot(deltaB,rbRot)/D;
                epsilon = dot(deltaB,rb)/D;
                kappa(b) = theta / (1+epsilon);
            end
        end

        for c = bulkCells

            cverts = Struct(t-1).Cdat(c).nverts;
            % Order verts
            rC = double([Struct(t-1).Vdat(cverts).vertxcoord;Struct(t-1).Vdat(cverts).vertycoord])';
            deltaRC = bsxfun(@minus,rC,mean(rC,1));

            theta = mod(atan2(deltaRC(:,2),deltaRC(:,1)),2*pi);
            [~,ind] = sort(theta);
            cverts = cverts(ind);
            rC = rC(ind,:);

            for ii = 1:length(cverts)

                p =  mod(ii,length(cverts)) + 1;
                m = mod(ii-2,length(cverts)) + 1;

                bondP = find( ((bVerts(1,:) == cverts(ii)) & (bVerts(2,:) == cverts(p))) | ...
                              ((bVerts(1,:) == cverts(p)) & (bVerts(2,:) == cverts(ii))) );
                
                bondM = find( ((bVerts(1,:) == cverts(ii)) & (bVerts(2,:) == cverts(m))) | ...
                              ((bVerts(1,:) == cverts(m)) & (bVerts(2,:) == cverts(ii))) );
                          
                delta1 = rC(p,:) - rC(ii,:); 
                delta1 = delta1/sqrt(sum(delta1.^2));
                delta2 = rC(m,:) - rC(ii,:); 
                delta2 = delta2/sqrt(sum(delta2.^2));

                vE = Struct(t-1).Vdat(cverts(ii)).nverts(~ismember(Struct(t-1).Vdat(cverts(ii)).nverts,cverts));

                bondE = find( ((bVerts(1,:) == cverts(ii)) & (bVerts(2,:) == vE)) | ...
                              ((bVerts(1,:) == vE) & (bVerts(2,:) == cverts(ii))) );

                deltaE = double([Struct(t-1).Vdat(vE).vertxcoord,Struct(t-1).Vdat(vE).vertycoord]);
                deltaE = deltaE - rC(ii,:); 
                deltaE = deltaE/sqrt(sum(deltaE.^2));

                angle1 = acos(dot(deltaE,delta1));
                angle2 = acos(dot(deltaE,delta2));

                chi{t-1}(n) = chi{t-1}(n) + log(sin(angle1)) - log(sin(angle2));
                
                if (~isempty(bondE))
                    chiEst{t-1}(n) = chiEst{t-1}(n) + cot(angle1)*(kappa(bondP)-kappa(bondE)) - cot(angle2)*(kappa(bondE)-kappa(bondM));
                end
                
            end
            n = n + 1;
        end
    end

end

