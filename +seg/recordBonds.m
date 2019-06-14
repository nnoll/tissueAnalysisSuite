function [ Struct ] = recordBonds( Struct, L )
    % Populate Struct.Bdat

    for t = 1:length(Struct)

        LDat = squeeze(L(:,:,t)); 
        bDat = (LDat == 0);

        rv = [Struct(t).Vdat.vertxcoord; Struct(t).Vdat.vertycoord];
        rv = double(round(rv));
        verts = zeros(size(bDat));
        verts(rv(2,:) + size(LDat,1)*(rv(1,:)-1)) = 1;

        % Inspect for fourfold
        if (any(strcmp('fourfold',fieldnames(Struct(t).Vdat))))
            nFourFold = sum([Struct(t).Vdat.fourfold]);
        else
            nFourFold = 0;
        end

        % Prepare bDat
        bDat(verts == 1) = 0;
        bEnd = bDat .* imdilate(verts, strel('disk',1));

        if (nFourFold > 0)
            for x = -1:1
                for y = -1:1
                    verts(rv(2,end-nFourFold+1:end)+y+size(LDat,1)*(rv(1,end-nFourFold+1:end)-1+x)) = 1;
                end
            end
        end

        bDat(verts == 1) = 0;
        bEnd = (bEnd .* bDat) + (bwmorph(bDat,'endpoints') .* bDat);

        bIndx = find(bEnd);
        clear rE
        [rE(2,:),rE(1,:)] = ind2sub(size(bDat),bIndx);
        D = pdist2(rE',rv');
        bCC = bwconncomp(bDat,4);
        bL = labelmatrix(bCC);
        endLabels = bL(bIndx);

        % This timepoints' bond data is given a struct
        Struct(t).Bdat(bCC.NumObjects) = struct('pix',[],'verts',[],'cells',[]);

        nB = 1;
        for b = 1:bCC.NumObjects
            endPoints = find(endLabels == b);

            v1 = 0; v2 = 0;
            if (length(endPoints) == 2)
                [~,v1] = min(D(endPoints(1),:));
                [~,v2] = min(D(endPoints(2),:));
                if (v1 == v2)
                    [sort1,ind1] = sort(D(endPoints(1),:));
                    [sort2,ind2] = sort(D(endPoints(2),:));
                    if (abs(sort1(1) - sort1(2)) <= sqrt(3))
                        v1 = ind1(2);
                    elseif (abs(sort2(1) - sort2(2)) <= sqrt(3))
                        v2 = ind2(2);
                    end
                end
            % elseif (length(endPoints) == 1)
            %     [~,tmp] = sort(D(endPoints,:),'ascend');
            %     v1 = tmp(1);
            %     v2 = tmp(2);
            end

            if (v1 > 0 && v2 > 0 && ismember(v2,Struct(t).Vdat(v1).nverts) && ~all(ismember([v1,v2],Struct(t).Cdat(1).nverts)))
                Struct(t).Bdat(nB).pix = bCC.PixelIdxList{b};
                Struct(t).Bdat(nB).verts = [v1;v2];
                Struct(t).Bdat(nB).cells = Struct(t).Vdat(v1).ncells(ismember(Struct(t).Vdat(v1).ncells,Struct(t).Vdat(v2).ncells));
                nB = nB + 1;
            end
        end

        if (nB <= bCC.NumObjects)
            Struct(t).Bdat(nB:bCC.NumObjects) = [];
        end

        bverts = [Struct(t).Bdat.verts];
        nB = length(Struct(t).Bdat);

        for v = 1:length(Struct(t).Vdat)

            vbonds1 = (bverts(1,:)==v);
            vbonds2 = (bverts(2,:)==v);

            kk = 1;
            for nv = Struct(t).Vdat(v).nverts

                nvbonds1 = (bverts(1,:)==nv);
                nvbonds2 = (bverts(2,:)==nv);

                bond1 = find(vbonds1.*nvbonds2);
                bond2 = find(vbonds2.*nvbonds1);

                if (~isempty(bond1))
                    Struct(t).Vdat(v).bond(kk:kk+length(bond1)-1) = bond1;
                    kk = kk + length(bond1);
                elseif (~isempty(bond2))
                    Struct(t).Vdat(v).bond(kk:kk+length(bond2)-1) = bond2;
                    kk = kk + length(bond2);
                elseif (~ismember(v,Struct(t).Cdat(1).nverts) && ~ismember(nv,Struct(t).Cdat(1).nverts))
                    nB = nB + 1;
                    [x,y] = plot.bresenham(double(Struct(t).Vdat(v).vertxcoord), ...
                            double(Struct(t).Vdat(v).vertycoord),double(Struct(t).Vdat(nv).vertxcoord),...
                            double(Struct(t).Vdat(nv).vertycoord));
                    Struct(t).Bdat(nB).pix = y + size(L,1)*(x-1);
                    Struct(t).Bdat(nB).verts = [v;nv];
                    Struct(t).Bdat(nB).cells = Struct(t).Vdat(v).ncells(ismember(Struct(t).Vdat(v).ncells,Struct(t).Vdat(nv).ncells));
                    bverts = horzcat(bverts,[v;nv]);

                    Struct(t).Vdat(v).bond(kk) = size(bverts,2);
                    vbonds1 = horzcat(vbonds1,1);
                    vbonds2 = horzcat(vbonds2,0);
                    kk = kk + 1;
                else
                    Struct(t).Vdat(v).bond(kk) = 0;
                    kk = kk + 1;
                end
            end

        end
        
        % For each cell, record the bonds
        for c = 2:length(Struct(t).Cdat)
            cverts = Struct(t).Cdat(c).nverts;
            if (length(cverts) > 1)
                rc = double([Struct(t).Vdat(cverts).vertxcoord;Struct(t).Vdat(cverts).vertycoord]);
                rc = rc';
                rc = bsxfun(@minus,rc,mean(rc,1));
                [theta] = mod(atan2(rc(:,2),rc(:,1)),2*pi);
                [~,ind] = sort(theta);
                cverts = cverts(ind);
                Struct(t).Cdat(c).orderednverts = cverts;
                if ( ~ismember(c,Struct(t).Cdat(1).ncells) && Struct(t).Cdat(c).all_threefold == 1 )
                    Struct(t).Cdat(c).bonds = zeros(1,length(cverts)-1);
                    for ii = 1:length(cverts)-1
                        if (~isempty(find(Struct(t).Vdat(cverts(ii)).nverts == cverts(ii+1))))
                            Struct(t).Cdat(c).bonds(ii) = Struct(t).Vdat(cverts(ii)).bond(find(Struct(t).Vdat(cverts(ii)).nverts == cverts(ii+1)));
                        else
                            Struct(t).Cdat(c).bonds(ii) = 0;
                        end
                    end
                end
            end
        end
    end

end
