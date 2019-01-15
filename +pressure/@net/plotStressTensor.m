function [ ] = plotStressTensor( this, k )
    % PLOT STRESS TENSOR 

    sigma = this.computeStressTensor();
    
    if (nargin > 1 || ((nargin == 2) && k > 0))
        [L,bndryCells] = this.returnSmoothingOperator();
        bulkCells = 1:size(sigma,1);
        bulkCells(bndryCells) = [];
        L = expm(L)^k;
        sigma(bulkCells,:) = L * sigma(bulkCells,:);
    end
    
    r = this.computePrimalVerts();
    [ faces ] = this.computeFaces( r );
    
    Rc = zeros(size(this.d0,2),2);
    pAxis = zeros(size(this.d0,2),2);
    sMajor = zeros(size(this.d0,2),1);
    sMinor = zeros(size(this.d0,2),1);
    stress = zeros(2,2,size(this.d0,2));

    R0 = .2*mean(sqrt(sum((this.d0*this.q).^2,2)));
    for c = 1:size(this.d0,2)
        nverts = faces(c,~isnan(faces(c,:)));
        Rc(c,:) = mean(r(nverts,:),1);
        stress(:,:,c) = [sigma(c,1),sigma(c,2);sigma(c,2),sigma(c,3)];
    end
    
   for c = 1:size(this.d0,2)
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

