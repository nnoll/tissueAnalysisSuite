function [ chi, chiEst ] = compareComptPressure( Struct, extCell )
%COMPARECOMPTPRESSURE Summary of this function goes here
%   Detailed explanation goes here

    for t = 1:length(Struct)
        
        [ ~, bulkCells, ~, bulkVerts ] = fitDual.returnGraph( Struct(t), extCell ); 

        chi{t} = zeros(length(bulkCells),1);
        chiEst{t} = zeros(length(bulkCells),1);
        n = 1;

        bVerts = [Struct(t).Bdat.verts];
        kappa = zeros(length(Struct(t).Bdat),1);
        for b = 1:length(Struct(t).Bdat)
            R = Struct(t).Bdat(b).radius;
            v1 = Struct(t).Bdat(b).verts(1);
            v2 = Struct(t).Bdat(b).verts(2);
            rb = double([Struct(t).Vdat(v1).vertxcoord;Struct(t).Vdat(v1).vertycoord]) - ...
                 double([Struct(t).Vdat(v2).vertxcoord;Struct(t).Vdat(v2).vertycoord]);
            r = sqrt(sum(rb.^2));
            rBar = Struct(t).Bdat(b).rBar;
            if (size(rBar,1) == 1)
                rBar = rBar';
            end
            S = sign( dot(double([Struct(t).Vdat(v1).vertxcoord;Struct(t).Vdat(v1).vertycoord])-rBar,[0,-1;1,0]*rb) );
            kappa(b) = -S*r./R;
        end

        for c = bulkCells
            
            cverts = Struct(t).Cdat(c).nverts;

            % Order verts
            rC = double([Struct(t).Vdat(cverts).vertxcoord;Struct(t).Vdat(cverts).vertycoord])';
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
                          
                if (bVerts(1,bondP) == cverts(ii))
                    Sp = -1;
                else
                    Sp = 1;
                end
                
                bondM = find( ((bVerts(1,:) == cverts(ii)) & (bVerts(2,:) == cverts(m))) | ...
                              ((bVerts(1,:) == cverts(m)) & (bVerts(2,:) == cverts(ii))) );
                          
                if (bVerts(1,bondM) == cverts(ii))
                    Sm = -1;
                else
                    Sm = 1;
                end
                
                delta1 = rC(p,:) - rC(ii,:); 
                delta1 = delta1/sqrt(sum(delta1.^2));
                delta2 = rC(m,:) - rC(ii,:); 
                delta2 = delta2/sqrt(sum(delta2.^2));

                vE = Struct(t).Vdat(cverts(ii)).nverts(~ismember(Struct(t).Vdat(cverts(ii)).nverts,cverts));

                bondE = find( ((bVerts(1,:) == cverts(ii)) & (bVerts(2,:) == vE)) | ...
                              ((bVerts(1,:) == vE) & (bVerts(2,:) == cverts(ii))) );
                          
                if (bVerts(1,bondE) == cverts(ii))
                    Se = -1;
                else
                    Se = 1;
                end      
                
                deltaE = double([Struct(t).Vdat(vE).vertxcoord,Struct(t).Vdat(vE).vertycoord]);
                deltaE = deltaE - rC(ii,:); 
                deltaE = deltaE/sqrt(sum(deltaE.^2));

                angle1 = acos(dot(deltaE,delta1));
                angle2 = acos(dot(deltaE,delta2));

                chi{t}(n) = chi{t}(n) + log(sin(angle1)) - log(sin(angle2));
                
                if (~isempty(bondE))
                    chiEst{t}(n) = chiEst{t}(n) + cot(angle1)*(Sp*kappa(bondP)-Se*kappa(bondE)) - cot(angle2)*(Se*kappa(bondE)-Sm*kappa(bondM));
                end
                
            end
            n = n + 1;
        end
    end

end

