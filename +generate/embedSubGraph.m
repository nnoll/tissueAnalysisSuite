function [ rv, E, rE ] = embedSubGraph( Struct, bulkVerts, extVerts, embed  )
    % EMBED SUB GRAPH 

    iVerts = [bulkVerts,extVerts];
    
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
%         rE(:,1) = embed(r0G(:,2)+size(embed,1)*(r0G(:,1)-1));
%         rE(:,2) = embed(r0G(:,2)+size(embed,1)*(r0G(:,1)-1)+size(embed,1)*size(embed,2));
%         rE(:,3) = embed(r0G(:,2)+size(embed,1)*(r0G(:,1)-1)+2*size(embed,1)*size(embed,2));
        [Xg,Yg] = meshgrid(1:size(embed,2),1:size(embed,1));
        rE(:,1) = interp2(Xg,Yg,embed(:,:,1),r0(:,1),r0(:,2));
        rE(:,2) = interp2(Xg,Yg,embed(:,:,2),r0(:,1),r0(:,2));
        rE(:,3) = interp2(Xg,Yg,embed(:,:,3),r0(:,1),r0(:,2));
    end
    
    rE = bsxfun(@minus,rE,mean(rE,1));
    [~,S,V] = svd(rE);
%     Vg = 100*V;
%     scatter3(rE(:,1),rE(:,2),rE(:,3))
%     hold all
%     quiver3(0,0,0,Vg(1,1),Vg(2,1),Vg(3,1),0)
%     quiver3(0,0,0,Vg(1,2),Vg(2,2),Vg(3,2),0)
%     quiver3(0,0,0,Vg(1,3),Vg(2,3),Vg(3,3),0)
%     axis equal
%     pause
    
    rv(:,1) = rE*V(:,1)/sqrt(sum(V(:,1).^2));
    rv(:,2) = rE*V(:,2)/sqrt(sum(V(:,2).^2));
    z = rE*V(:,3)/sqrt(sum(V(:,3).^2));
    E = 2*S(3,3)/(S(1,1)+S(2,2));
    Rm = mean(r0,1);
    r0 = bsxfun(@minus,r0,Rm);
    
    % Fit rotation between points
    C = r0'*rv;
    [Ur,~,Vr] = svd(C);
    R = Ur*Vr';
%     if (abs(det(Ur*Vr') -1 ) < 1e-3)
%         R = Ur*eye(2)*Vr';
%     else
%         R = Ur*diag([1,-1])*Vr'; 
%     end
    rv = rv*R';
    
%     r0 = bsxfun(@plus,r0,Rm);
%     hold on
%     for v = bulkVerts
%        nVerts = Struct.Vdat(v).nverts;
%        nVerts = nVerts(ismember(nVerts,iVerts));
%        for nv = nVerts
%            ind = [find(iVerts==v),find(iVerts==nv)];
%            plot3(rv(ind,1),rv(ind,2),z(ind)+.1,'LineWidth',2,'Color','b')
%        end
%     end
%     scatter3(rv(:,1),rv(:,2),z+.1,'b','filled')
% 
%     for v = bulkVerts
%        nVerts = Struct.Vdat(v).nverts;
%        nVerts = nVerts(ismember(nVerts,iVerts));
%        for nv = nVerts
%            ind = [find(iVerts==v),find(iVerts==nv)];
%            plot(rv(ind,1),rv(ind,2),'LineWidth',2,'Color','k')
%        end
%     end
%     scatter(rv(:,1),rv(:,2),'k','filled')

    rv = bsxfun(@plus,rv,Rm);
   
%     rv
%     plot.skel(Struct,'k',0);
%     scatter(r0(:,1),r0(:,2),'b','filled')
%     hold on
%     scatter(rv(:,1),rv(:,2),'c','filled')
%     rE([17,20],:)
%     r0([17,20],:)
%     rv([17,20],:)
%     quiver(r0(:,1),r0(:,2),rv(:,1)-r0(:,1),rv(:,2)-r0(:,2),0)
%     pause
end

