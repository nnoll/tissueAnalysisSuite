function [ E, dE ] = sphericalEnergy( X, Lbar )
%SPHERICALENERGY Summary of this function goes here
%   Detailed explanation goes here

    az = X(:,1);
    el = X(:,2);
    rv = zeros(size(X,1),3);
    rv(:,1) = cos(az).*sin(el);
    rv(:,2) = sin(az).*sin(el);
    rv(:,3) = cos(el);
    Tri = convhulln(rv);
    [ d0 ] = generate.primalBonds( triangulation(Tri,rv) );
    
    L0 = .5*Lbar*ones(size(d0,1),1);
    rB = d0*rv; 
    rB = bsxfun(@rdivide,rB,sqrt(sum(rB.^2,2)));
%     elB = .5*abs(d0)*el;
%     azB = .5*abs(d0)*az;
%     delta = cos(elB).*cos(azB).*rB(:,1) + cos(elB).*sin(azB).*rB(:,2) - sin(elB).*rB(:,3);
%     delta = delta.^2;
    delta = rB(:,2).*rB(:,1);
    
    L0 = L0 + .1*delta;
%     clf
%     for b = 1:size(d0,1)
%        v = find(d0(b,:));
%        plot3(rv(v,1),rv(v,2),rv(v,3),'Color',[delta(b)/max(delta),0,0])
%        hold on
%     end
%     pause
    L = sqrt(sum((d0*rv).^2,2));
    E = .5*sum((L-L0).^2);
    
    if (nargout > 1)
       rB = d0*rv;
       D = L; D(D==0) = 1;
       rB = bsxfun(@rdivide,rB,D);
       T = bsxfun(@times,L-L0,rB);
       Fnet = d0'*T;
       dE(:,1) = cos(az) .* Fnet(:,2) - sin(az) .* Fnet(:,1);
       dE(:,2) = cos(el) .* cos(az) .* Fnet(:,1) + cos(el) .* sin(az) .* Fnet(:,2) - sin(el) .*Fnet(:,3);
    end
    
end

