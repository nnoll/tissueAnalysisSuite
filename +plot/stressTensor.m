function [ ] = stressTensor( Struct, L, smoothSize )

    [ ~, ~, ~, ~, iCells ] = fitDual.ATN.computeDiffOperators( Struct, 1 );  

    Rc = zeros(length(iCells),2);
    pAxis = zeros(length(iCells),2);
    sMajor = zeros(length(iCells),1);
    sMinor = zeros(length(iCells),1);
    stress = zeros(2,2,length(iCells));
    
    S = regionprops(L,'Area');
    A = [S.Area];
    R0 = .3*nanmean(sqrt(A(2:end)/pi));

    for c = 1:length(iCells)
        Rc(c,:) = Struct.Cdat(iCells(c)).centroid.coord;
        if ( ~isempty(Struct.Cdat(iCells(c)).stress))
            stress(:,:,c) = Struct.Cdat(iCells(c)).stress .* A(iCells(c));
        end
    end
    
    stress(isnan(stress)) = 0;
    stress(isinf(stress)) = 0;
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
        [V,D] = eig(stress(:,:,c));
        pAxis(c,:) = V(:,2);
        pAxis(c,:) = pAxis(c,:)/sqrt(sum(pAxis(c,:).^2));
        sMajor(c) = D(2,2);
        sMinor(c) = D(1,1);
    end
    
    t = linspace(0,2*pi,100);
    
    P = sqrt(median(sMajor.*sMinor));
    sMajor = (sMajor / P);
    sMajor(sMajor>3) = 3;
    sMinor = (sMinor / P);
    sMinor(sMinor<0) = 0;
    sMinor(sMinor>3) = 3;
    
    sMinor = R0*sMinor;
    sMajor = R0*sMajor;
    
    xM1 = Rc + bsxfun(@times,sMajor,pAxis);
    xM2 = Rc - bsxfun(@times,sMajor,pAxis);
    
    X = sMajor * cos(t);
    Y = sMinor * sin(t);
    
    w = atan2( xM2(:,2) - xM1(:,2), xM2(:,1) - xM1(:,1) );
    
    x = bsxfun(@plus,Rc(:,1),bsxfun(@times,X,cos(w)) - bsxfun(@times,Y,sin(w)))';
    y = bsxfun(@plus,Rc(:,2),bsxfun(@times,X,sin(w)) + bsxfun(@times,Y,cos(w)))';
    
    plot(x,y,'r','LineWidth',2)

end

