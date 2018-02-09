function [ chi, chiC, chiCV ] = buildControlDistribution( Struct, extCell )
    % BUILD CONTROL DISTRIBUTION 
    
    for t = 1:length(Struct)
        
        [ ~, bulkCells, ~, bulkVerts ] = fitDual.returnGraph( Struct(t), extCell ); 

        % Build angle distribution
        latticeAngles = zeros(length(bulkVerts),3);
        n = 1;
        for v = bulkVerts
            rV = double([Struct(t).Vdat(v).vertxcoord;Struct(t).Vdat(v).vertycoord]);
            
            nv = Struct(t).Vdat(v).nverts;
            rNV1 = double([Struct(t).Vdat(nv(1)).vertxcoord;Struct(t).Vdat(nv(1)).vertycoord]);
            rNV2 = double([Struct(t).Vdat(nv(2)).vertxcoord;Struct(t).Vdat(nv(2)).vertycoord]);
            rNV3 = double([Struct(t).Vdat(nv(3)).vertxcoord;Struct(t).Vdat(nv(3)).vertycoord]);
            
            delta1 = rNV1 - rV; delta1 = delta1/sqrt(sum(delta1.^2));
            delta2 = rNV2 - rV; delta2 = delta2/sqrt(sum(delta2.^2));
            delta3 = rNV3 - rV; delta3 = delta3/sqrt(sum(delta3.^2));
            
            cos12 = dot(delta1,delta2);
            cos23 = dot(delta2,delta3);
            cos31 = dot(delta3,delta1);
            
            latticeAngles(n,:) = acos([cos12,cos23,cos31]);
            n = n + 1;
        end
        
        chi{t} = zeros(length(bulkCells),1);
        chiC{t} = zeros(length(bulkCells),1);
        chiCV{t} = zeros(length(bulkCells),1);
        
        controlAngles = [];
        n = 1;
        
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
                
                delta1 = rC(p,:) - rC(ii,:); delta1 = delta1/sqrt(sum(delta1.^2));
                delta2 = rC(m,:) - rC(ii,:); delta2 = delta2/sqrt(sum(delta2.^2));
                
                vE = Struct(t).Vdat(cverts(ii)).nverts(~ismember(Struct(t).Vdat(cverts(ii)).nverts,[cverts]));
                
                if (length(vE) == 2)
                    vE = vE(1);
                end
                deltaE = double([Struct.Vdat(vE).vertxcoord,Struct.Vdat(vE).vertycoord]);
                deltaE = deltaE - rC(ii,:); deltaE = deltaE/sqrt(sum(deltaE.^2));
                
                angle1 = acos(dot(deltaE,delta1));
                angle2 = acos(dot(deltaE,delta2));
                chi{t}(n) = chi{t}(n) + log(sin(angle1)) - log(sin(angle2));
                
                cellAngle = acos((dot(delta1,delta2)));
                contin = 1;
                while (contin == 1)
                    randAngle1 = randsample(latticeAngles(:),1);
                    randAngle2 = 2*pi - randAngle1 - cellAngle;
                    contin = (randAngle1 > pi) || (randAngle2 > pi);
                end
                
                controlAngles = [controlAngles,[randAngle1]];
                
                randNumber = rand(1);
                
                if (randNumber < .5)
                    chiC{t}(n) = chiC{t}(n) + log(sin(randAngle1)) - log(sin(randAngle2));
                else
                    chiC{t}(n) = chiC{t}(n) + log(sin(randAngle2)) - log(sin(randAngle1));
                end
                
                randV = randsample(size(latticeAngles,1),1);
                randInd = randsample(1:3,2);
                chiCV{t}(n) = chiCV{t}(n) + log(sin(latticeAngles(randV,randInd(1)))) - log(sin(latticeAngles(randV,randInd(2))));
                
            end
            
            n = n + 1;
            
        end
        
        cdfplot(latticeAngles(:))
        hold all
        cdfplot(controlAngles(:))
        [~,p] = kstest2(latticeAngles(:),controlAngles(:))
%         angle = {latticeAngles(:),controlAngles};
% % %         plot.nhist(angle,'binFactor',3,'LineWidth',2);
% %         histfit(latticeAngles(:),20)
% %         hold all
% %         histfit(controlAngles,20)
%          plot.niceHist( angle, 42 )
%         
    end

end

