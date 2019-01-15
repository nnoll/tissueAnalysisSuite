function [ dT ] = calculateDivTensor( Mesh, T )

    % Flatten tensor T
    TF = zeros(size(T,1),9);
    TF(:,1) = T(:,1,1);
    TF(:,2) = T(:,1,2);
    TF(:,3) = T(:,1,3);
    TF(:,4) = T(:,2,1);
    TF(:,5) = T(:,2,2);
    TF(:,6) = T(:,2,3);
    TF(:,7) = T(:,3,1);
    TF(:,8) = T(:,3,2);
    TF(:,9) = T(:,3,3);
    
    X = Mesh.circumcenters;
    tri = (1:size(Mesh.Triangulation,1))';
    Rbary = Mesh.cartToBary(tri, X);
    Rref = Mesh.baryToRef(Rbary);
    interpData = {uint32(tri),Rref};
    
    INTERP_2 = interpTGrad(Mesh.X,uint32(Mesh.Triangulation),[],[],interpData,uint32(Mesh.Triangulation),TF);

    % Unpack gradient.
    tGrad = zeros(3,3,size(X,1),3);
    for ii = 1:3
        tGrad(1,1,:,ii) = INTERP_2.DATA{1,ii};
        tGrad(1,2,:,ii) = INTERP_2.DATA{2,ii};
        tGrad(1,3,:,ii) = INTERP_2.DATA{3,ii};
        tGrad(2,1,:,ii) = INTERP_2.DATA{4,ii};
        tGrad(2,2,:,ii) = INTERP_2.DATA{5,ii};
        tGrad(2,3,:,ii) = INTERP_2.DATA{6,ii};
        tGrad(3,1,:,ii) = INTERP_2.DATA{7,ii};
        tGrad(3,2,:,ii) = INTERP_2.DATA{8,ii};
        tGrad(3,3,:,ii) = INTERP_2.DATA{9,ii};
    end
    
    dT = squeeze(tGrad(1,:,:,1)) + squeeze(tGrad(2,:,:,2)) + squeeze(tGrad(3,:,:,3));

end

