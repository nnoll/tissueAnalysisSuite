function [ rE, r0 ] = embedGraph( Struct, embed  )
    % EMBED SUB GRAPH 

    iVerts = 1:length(Struct.Vdat);
    
    r0 = zeros(length(iVerts),2);
    for ii = 1:length(iVerts)
        r0(ii,1) = double(Struct.Vdat(iVerts(ii)).vertxcoord);
        r0(ii,2) = double(Struct.Vdat(iVerts(ii)).vertycoord);
    end
    
    % Compute positions in R3
    rE = zeros(length(iVerts),3);

    if (isa(embed,'function_handle'))
        rE = embed(r0(:,1),r0(:,2));
    else
        [Xg,Yg] = meshgrid(1:size(embed,2),1:size(embed,1));
        rE(:,1) = interp2(Xg,Yg,embed(:,:,1),r0(:,1),r0(:,2));
        rE(:,2) = interp2(Xg,Yg,embed(:,:,2),r0(:,1),r0(:,2));
        rE(:,3) = interp2(Xg,Yg,embed(:,:,3),r0(:,1),r0(:,2));
    end
  
end

