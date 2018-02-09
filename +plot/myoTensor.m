function [ ] = myoTensor( Struct, L, smoothSize )

    [ ~, ~, ~, ~, iCells ] = fitDual.ATN.computeDiffOperators( Struct, 1 );  
    Rc = zeros(length(iCells),2);
    radius = zeros(length(iCells),1);
    pAxis = zeros(length(iCells),2);
    sAniso = zeros(length(iCells),1);
    stress = zeros(2,2,length(iCells));
    
    for c = 1:length(iCells)
        Rc(c,:) = Struct.Cdat(iCells(c)).centroid.coord;
        stress(:,:,c) = Struct.Cdat(iCells(c)).myo;
    end
    
    S = regionprops(L,'Area');
    A = [S(2:end).Area];
    R0 = .5*mean(sqrt(A/pi));
    if (nargin > 2 && smoothSize > 0)
        Dc = pdist2(Rc,Rc);
        SmK = exp(-Dc.^2/(2*(smoothSize*R0)^2));
        SmK = bsxfun(@rdivide,SmK,sum(SmK,2));
        stress(1,1,:) = SmK*squeeze(stress(1,1,:));
        stress(1,2,:) = SmK*squeeze(stress(1,2,:));
        stress(2,1,:) = SmK*squeeze(stress(2,1,:));
        stress(2,2,:) = SmK*squeeze(stress(2,2,:));
    end
    
    for c = 1:length(iCells)
        radius(c) = trace(stress(:,:,c));
        [V,D] = eig(stress(:,:,c));
        pAxis(c,:) = V(:,2);
        pAxis(c,:) = pAxis(c,:)/sqrt(sum(pAxis(c,:).^2));
        sAniso(c) = (D(2,2) - D(1,1))/(D(1,1)+D(2,2));
    end
    
    % Scale radius to order of cell size for visibility.
    pAxis = bsxfun(@times,pAxis,R0*sAniso/mean(sAniso));
    radius = R0*radius/mean(radius);
    
    viscircles(Rc,radius,'LineWidth',2,'DrawBackgroundCircle',0,'EdgeColor',[0,0,1]);
    hold on
    plot([Rc(:,1)+pAxis(:,1),Rc(:,1)-pAxis(:,1)]',[Rc(:,2)+pAxis(:,2),Rc(:,2)-pAxis(:,2)]','LineWidth',3,'Color',[1,0,0])
    
end

