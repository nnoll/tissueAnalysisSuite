function [ N, b0 ] = computeNormal( Struct )
    % COMPUTE NORMAL 

    rC = [Struct.Cdat.centroid];
    rC = vertcat(rC.coord);
    
    embed = Struct.embed;
    rE = zeros(size(rC,1),3);

    [Xg,Yg] = meshgrid(1:size(embed,2),1:size(embed,1));
    rE(:,1) = interp2(Xg,Yg,embed(:,:,1),rC(:,1),rC(:,2));
    rE(:,2) = interp2(Xg,Yg,embed(:,:,2),rC(:,1),rC(:,2));
    rE(:,3) = interp2(Xg,Yg,embed(:,:,3),rC(:,1),rC(:,2));

    [ Tri, bC, eC, b0 ] = fitDual.returnGraph(Struct,1);

    iC = [bC,eC];
    N = zeros(length(b0),3);
    for v = 1:length(b0)
        N(v,:) = cross( rE(iC(Tri(v,2)),:) - rE(iC(Tri(v,1)),:), rE(iC(Tri(v,3)),:) - rE(iC(Tri(v,1)),:) );
        N(v,:) = N(v,:) / sqrt(sum(N(v,:).^2));
    end
    
end

