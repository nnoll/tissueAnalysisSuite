function [ Fnet, rE ] = computeForceBalance( Struct )
    
    % COMPUTE FORCE BALANCE 
    [ rE ] = generate.embedGraph( Struct, Struct.embed  );
    [ ~, dV, ~, i0 ] = fitDual.ATN.computeDiffOperators( Struct, 1 );
    rE = rE(i0,:);
    [ N, b0 ] = embedSpace.computeNormal( Struct );
    
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
    
    Tvec = bsxfun(@times,T,rB);
    
    Fnet = -dV'*Tvec;
    Fnet(1:length(b0),:) = Fnet(1:length(b0),:) - bsxfun(@times,dot(Fnet(1:length(b0),:),N,2),N);
end

