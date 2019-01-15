function [ I ] = maxABIC(x, A, b, B, g, m0, n, m )
    % MAX ABIC 
    
    S = sparse([ [A,b]; [sqrt(x)*B,sqrt(x)*g] ]);
    [~,R,E] = qr(S,'matrix');
    Rdiag = R(1:(m+1),(size(R,2)-m):size(R,2));
    Rdiag = E*Rdiag*E';
    Rdiag = full(diag(Rdiag));
    E = Rdiag(end)^2;
    
    lambda = Rdiag(1:end-1);
    lambda = abs(lambda(lambda~=0));
    detFull = sum(log(lambda));    

    I = -.5*((size(B,1)-m0)*log(x) - (n-m0+1)*(log(E)) - detFull); 
    
end

