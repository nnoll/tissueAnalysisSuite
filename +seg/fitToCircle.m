function [ y, E ] = fitToCircle( rData, nB, x0, D )
%FITTOCIRCLE 

    delta = bsxfun(@minus,rData,x0');
    IP = delta(:,1)*nB(1) + delta(:,2)*nB(2);
    L0 = D/2;
    
    A = 2*sum(IP.^2);
    B = sum(  (sum(delta.^2,2) - (D/2)^2 ) .* IP );
    
    y0 = B/A;
    optimset = optimoptions('fminunc','Display','none','Algorithm','quasi-newton');

%     f = @(x) sum( ((sum(bsxfun(@minus,delta,x*nB').^2,2) - (x^2+L0^2)).^2) ) / (2*sum((sum(bsxfun(@minus,delta,x*nB').^2,2))) - (x^2+L0^2));
    f = @(x) mean( ( sqrt(sum(bsxfun(@minus,delta,x*nB').^2,2)) - sqrt((x^2+L0^2)) ).^2 );
    
    if (~isnan(y0))
        [y,E] = fminunc(f,y0,optimset);
    else
        [y,E]= fminunc(f,0,optimset);
    end
   
end

