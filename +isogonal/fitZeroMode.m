function [ theta, bulkCells, R, rgb, F ] = fitZeroMode( Struct, vPair, cPair, L, t, displayFlag )
    %CONFORMAL_FLUCT This function compares a time series of vertex positions
    %to a reference lattice (previous time point) and fits the vertex movement 
    %to its conformal mode via least squares. Will return the time trajectory 
    %of such modes for each cell as well as the fraction of explained variance 
    %by this mode.

    %INPUTS:
    %1. Struct - Recorded observables for vertices, cells, and bonds.
    %2. v_LUT - Tracked vertex identities in time.
    %3. [t(1),t(2)] - initial and final time-point
    %4. L - watershedded image

    %OUTPUTS:
    %1. theta - Value of fitted mode at all times
    %2. R - Fraction of explained variance at all times

    %% Obtain `best-fit' cellular angles that obey compatability.
    if (~displayFlag)
        rgb = 0;
    end

    [ Q, S, dC, dV, bulkVerts, bulkCells ] = isogonal.getTriGeo( Struct(t(1)) );

    intIndx = find(sum(abs(dV),1) > 1);
    intVerts = bulkVerts(intIndx);
    
    %% Get displacement vectors of bulk vertices
    R0(:,1) = [Struct(t(1)).Vdat.vertxcoord]; R0(:,2) = [Struct(t(1)).Vdat.vertycoord];
    R(:,1) = [Struct(t(2)).Vdat.vertxcoord]; R(:,2) = [Struct(t(2)).Vdat.vertycoord];

    deltaR = zeros(length(intVerts),2);
    [vTrack,cTrack] = track.propagateVertex(vPair,cPair,Struct,t);

    for v = 1:length(vTrack)
        if (vTrack(v) > 0 && ismember(vTrack(v),intVerts))
            ind = (intVerts==vTrack(v)); 
            deltaR(ind,:) = R(v,:) - R0(vTrack(v),:); 
        end
    end
    
    emptyData = find(sum(deltaR.^2,2)==0);
    if (~isempty(emptyData))
        deltaR(emptyData,:) = isogonal.fitMissingDisplacements(Struct,cTrack,t,R0(intVerts(emptyData),:));
    end
    
    % Subtract off overall translation.
    Rtrans = mean(deltaR,1);
    deltaR = bsxfun(@minus,deltaR,Rtrans);
    
    %% Construct the structure factor matrix to invert.
    F = zeros(2*length(intVerts),length(bulkCells));
    
    % Build tension vectors.
%     RR0 = R0;
    R0 = R0(bulkVerts,:);
    TB = (dC*Q)*[0,-1;1,0];
    rB = dV*R0;
    rotS = sign(dot(rB,TB,2));
    TB = bsxfun(@times,rotS,TB);
    
    for v = 1:length(intVerts)
        conBonds = find(dV(:,intIndx(v)));
        conCells = (sum(abs(dC(conBonds,:)),1)>0);
        for b = conBonds'
            bondCells = abs(dC(b,:));
            isoCells = conCells;
            isoCells(bondCells==1) = 0;
            F(v,isoCells) = (dV(b,intIndx(v)) * TB(b,1))/S(intIndx(v));
            F(v+length(intVerts),isoCells) = (dV(b,intIndx(v)) * TB(b,2))/S(intIndx(v));
        end
    end
    
%     delVerts = [emptyData;emptyData+length(intVerts)];
%     Fsmall = F;
%     Fsmall(delVerts,:) = [];
    Delta = deltaR(:);
%     Delta(delVerts) = [];
    theta = pinv(F)*Delta;
    rEst = F*theta;
    
%     deltaR(emptyData,:) = [];
    Rtrans = ones(size(deltaR,1),1)*Rtrans;
    rEst = rEst + Rtrans(:);
    rEst = reshape(rEst,length(rEst)/2,2);
    deltaR = deltaR + Rtrans;
    
%     scatter(RR0(:,1),RR0(:,2))
%     hold all
%     quiver(RR0(intVerts,1),RR0(intVerts,2),deltaR(:,1),deltaR(:,2));
%     quiver(RR0(intVerts,1),RR0(intVerts,2),rEst(:,1),rEst(:,2));
% 
%     scatter(RR0(intVerts(emptyData),1),RR0(intVerts(emptyData),2),'r.')
%     pause

    nStruct = Struct(t(1));
    nnStruct = Struct(t(1));
    for ii = 1:length(intVerts)
        nStruct.Vdat(intVerts(ii)).vertxcoord = nStruct.Vdat(intVerts(ii)).vertxcoord + deltaR(ii,1);
        nStruct.Vdat(intVerts(ii)).vertycoord = nStruct.Vdat(intVerts(ii)).vertycoord + deltaR(ii,2);
        nnStruct.Vdat(intVerts(ii)).vertxcoord = nnStruct.Vdat(intVerts(ii)).vertxcoord + rEst(ii,1);
        nnStruct.Vdat(intVerts(ii)).vertycoord = nnStruct.Vdat(intVerts(ii)).vertycoord + rEst(ii,2);
    end
    
    if (displayFlag)
        rgb(:,:,1) = plot.skel(Struct(t(1)),size(L,2),size(L,1),1);
        rgb(:,:,2) = plot.skel(nnStruct,size(L,2),size(L,1),1);
        rgb(:,:,3) = plot.skel(nStruct,size(L,2),size(L,1),1);
    end

    resid = sqrt(sum((rEst-deltaR).^2,2));

    %% Calculate bond length
    R = 1 - median(resid)/mean(sqrt(sum(deltaR.^2,2)));
    
    end


