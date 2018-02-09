function [ nStruct ] = removeFourFold( Struct )
%REMOVEFOURFOLD 
    
    nStruct = Struct;
    for t = 1:length(Struct)
        
        NV = length(Struct(t).Vdat);
        NB = length(Struct(t).Bdat);
        
        for v = 1:length(Struct(t).Vdat)
           if (length(nStruct(t).Vdat(v).nverts) == 4 && ~ismember(1,nStruct(t).Vdat(v).ncells)) % Four-fold vertex, replace!
               
               nVerts = nStruct(t).Vdat(v).nverts;
               nBonds = nStruct(t).Vdat(v).bond;
               
               r1 = double([nStruct(t).Vdat(nVerts(1)).vertxcoord;nStruct(t).Vdat(nVerts(1)).vertycoord]);
               r2 = double([nStruct(t).Vdat(nVerts(2)).vertxcoord;nStruct(t).Vdat(nVerts(2)).vertycoord]);
               r3 = double([nStruct(t).Vdat(nVerts(3)).vertxcoord;nStruct(t).Vdat(nVerts(3)).vertycoord]);
               r4 = double([nStruct(t).Vdat(nVerts(4)).vertxcoord;nStruct(t).Vdat(nVerts(4)).vertycoord]);
               
               rV = double([nStruct(t).Vdat(v).vertxcoord;nStruct(t).Vdat(v).vertycoord]);
               
               R = [r1,r2,r3,r4];
               Rcom = mean(R,2);
               R = bsxfun(@minus,R,Rcom);
               
               I = R*R';
               [V,~] = eig(I);

               axis = V(:,2); 
               axis = axis / sqrt(sum(axis.^2));
               
%                clf
%                plot.skel(nStruct(t),'k',0)
%                hold all
%                scatter(rV(1),rV(2),'b','filled')
%                scatter(R(1,:)+Rcom(1),R(2,:)+Rcom(2),'g','filled')
%                quiver(rV(1),rV(2),axis(1),axis(2),50,'Color','c')
%                quiver(rV(1),rV(2),-axis(1),-axis(2),50,'Color','c')
%                pause
               
                rV1 = rV + axis/2;
                rV2 = rV - axis/2;
                
                [~,ind] = sort( (R(1,:)*axis(1) + R(2,:)*axis(2)) );
                posVerts = zeros(size(ind));
                posVerts(ind(3:4)) = 1;
                posVerts = posVerts == 1;
                negVerts = ~posVerts;
                
                nStruct(t).Vdat(v).vertxcoord = rV1(1);
                nStruct(t).Vdat(v).vertycoord = rV1(2);
                nStruct(t).Vdat(v).nverts = [nVerts(posVerts),NV+1];
                nStruct(t).Vdat(v).fourfold = 0;
                nStruct(t).Vdat(v).old_fourfold = 1;
                
                posCell = nStruct(t).Vdat(v).ncells;
                for nv = nVerts(posVerts)
                   posCell = posCell(ismember(posCell,nStruct(t).Vdat(nv).ncells)); 
                end
                
                nStruct(t).Vdat(NV+1).vertxcoord = rV2(1);
                nStruct(t).Vdat(NV+1).vertycoord = rV2(2);
                nStruct(t).Vdat(NV+1).nverts = [nVerts(negVerts),v];
                nStruct(t).Vdat(NV+1).old_fourfold = 1;

                negCell = nStruct(t).Vdat(v).ncells;
                for nv = nVerts(negVerts)
                   nStruct(t).Vdat(nv).nverts(nStruct(t).Vdat(nv).nverts == v) = NV+1;
                   negCell = negCell(ismember(negCell,nStruct(t).Vdat(nv).ncells)); 
                end
                
                jointCells = nStruct(t).Vdat(v).ncells;
                jointCells = jointCells(~ismember(jointCells,[posCell,negCell]));
                
                nStruct(t).Vdat(NV+1).ncells = [negCell,jointCells];
                nStruct(t).Vdat(v).ncells = [posCell,jointCells];
                
                % Update bond data structure
                negBonds = nBonds(negVerts);
                posBonds = nBonds(posVerts);
                
                if (negBonds(1) > 1)
                    nStruct(t).Bdat(negBonds(1)).verts(nStruct(t).Bdat(negBonds(1)).verts==v) = NV + 1;
                end
                
                if (negBonds(2) > 1)
                    nStruct(t).Bdat(negBonds(2)).verts(nStruct(t).Bdat(negBonds(2)).verts==v) = NV + 1;
                end
                
                nStruct(t).Bdat(NB+1).verts = [v,NV+1];
                nStruct(t).Bdat(NB+1).cells = jointCells;
                nStruct(t).Bdat(NB+1).pix = []; %round(rV(2)) + (round(rV(1))-1)*dim(1);
                nStruct(t).Bdat(NB+1).radius = inf;
                nStruct(t).Bdat(NB+1).rBar = [inf;inf];
                
                nStruct(t).Vdat(v).bond = [posBonds,NB+1];
                nStruct(t).Vdat(NV+1).bond = [negBonds,NB+1];
                
                NV = NV + 1;
                NB = NB + 1;
           end
        end
        
        %% Update cell data structure
        for c = 1:length(nStruct(t).Cdat)
            nStruct(t).Cdat(c).numv = 0;
            nStruct(t).Cdat(c).nverts = [];
            nStruct(t).Cdat(c).bonds = [];
            nStruct(t).Cdat(c).all_threefold = [];
        end
        
        for v = 1:length(nStruct(t).Vdat)
           for nc = nStruct(t).Vdat(v).ncells 
                nStruct(t).Cdat(nc).nverts = [nStruct(t).Cdat(nc).nverts,v];
           end
        end
        
        for c = 1:length(nStruct(t).Cdat)
            nStruct(t).Cdat(c).numv = length(nStruct(t).Cdat(c).nverts);
        end
        
        for v = 1:length(nStruct(t).Vdat)
            nStruct(t).Vdat(v).vertxcoord = double(nStruct(t).Vdat(v).vertxcoord);
            nStruct(t).Vdat(v).vertycoord = double(nStruct(t).Vdat(v).vertycoord);
        end
    end
    
    nStruct = seg.threefold_cell(nStruct);

end

