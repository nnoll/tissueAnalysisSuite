function [ sigma, Z, Phi ] = computeStress( Struct )
    % TAKE CONTINUUM LIMIT 

    [ rE, r0 ] = generate.embedGraph( Struct, Struct.embed  );
    [ ~, dV, ~, i0 ] = fitDual.ATN.computeDiffOperators( Struct, 1 );
    rE = rE(i0,:); r0 = r0(i0,:);

    i1 = generate.bondMap(Struct);
    dV = sparse(dV);
    
    T = zeros(length(i1{1}));
    for b = 1:length(i1{1})
       if (i1{1}(b)>0)
           T(b) = Struct.Bdat(i1{1}(b)).tension;
       else
           T(b) = 1;
       end
    end
    
    T = T / mean(T);
    rB = dV*rE;
    rB = bsxfun(@rdivide,rB,sqrt(sum(rB.^2,2)));
    
    rBProj = dV*r0;
    rBProj = bsxfun(@rdivide,rBProj,sqrt(sum(rBProj.^2,2)));
    
    R = .5*abs(dV)*rE;
    D = pdist2(R,R);
    L = exp( -D.^2 / (2*50^2));
    L = bsxfun(@rdivide,L,sum(L,2));
    
    stress = zeros(size(rB,1),2,2);
    for ii = 1:2
        for jj = 1:2
            stress(:,ii,jj) = T .* rBProj(:,ii) .* rBProj(:,jj);
        end
    end
    Rp = .5*abs(dV)*r0; % Position within tangent plane
    
    % Coarse-grain in embedding space.
    for ii = 1:2
        for jj = 1:2
            stress(:,ii,jj) = L*stress(:,ii,jj);
        end
    end
    
    
    % Interpolate stress into tangent plane.
    [Z,Phi] = meshgrid(linspace(1,1738,218),linspace(1,2050,257));
    sigma = zeros(size(Z,1),size(Z,2),3);
    sigma(:,:,1) = griddata(Rp(:,1),Rp(:,2),stress(:,1,1),Z,Phi,'natural');
    sigma(:,:,2) = griddata(Rp(:,1),Rp(:,2),stress(:,1,2),Z,Phi,'natural');
    sigma(:,:,3) = griddata(Rp(:,1),Rp(:,2),stress(:,2,2),Z,Phi,'natural');

end

