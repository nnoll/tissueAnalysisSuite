function [ theta, q, p ] = fitThetaModes( q, p, Struct, extCell )
    % FIT THETA MODES 

    [ tri, bulk0, ext0, ~, ~ ] = fitDual.returnGraph( Struct, extCell );
    i0 = [ bulk0, ext0 ];
    
    [ d0, ~, ~, ~, ~, involvedBonds ] = fitDual.AFN.returnBonds( Struct, i0 );
%     pTri = fitDual.AFN.returnPTri(q,p,tri);

    R = zeros(size(d0,1),1);
    Qsq = sum((d0*q).^2,2);
    Rflat = zeros(size(d0,1),1);
    
    x = bsxfun(@rdivide,d0*bsxfun(@times,p,q),d0*p);
    for b = 1:size(d0,1)
        v1 = Struct.Bdat(involvedBonds(b)).verts(1);
        v2 = Struct.Bdat(involvedBonds(b)).verts(2);
        r1 = double([Struct.Vdat(v1).vertxcoord,Struct.Vdat(v1).vertycoord]);
        r2 = double([Struct.Vdat(v2).vertxcoord,Struct.Vdat(v2).vertycoord]);
        R(b) = mean( [sqrt(sum((r1-x(b,:)).^2)) , sqrt(sum((r2-x(b,:)).^2))] ).^2;
        Rflat(b) = p(d0(b,:)==1)*p(d0(b,:)==-1)*Qsq(b);
    end
    
    dP = d0*p;
    R = R.*dP.^2;

    A = bsxfun(@times,dP,d0);
%     L = [A;[ones(1,size(d0,2)),0]];
%     b = [R;0];
    b = Rflat - R;
    theta = [A;ones(1,size(A,2))] \ [b;0];
%     x0 = sparse(L) \ b;
%     optimset = optimoptions('lsqlin','Display','final','MaxIter',500,'Algorithm','active-set'); 
%     c = zeros(size(Rflat));
%     theta = theta \ b;
%     [theta,~,~,eF] = lsqlin(full(A),b,full(A),Rflat,[],[],[],[],[],optimset);

%     L = -sparse(bsxfun(@times,dP,d0));
%     b = R - Rflat;
%     
%     A = bsxfun(@times,dP,d0);
%     c = Rflat;
%     energyFunc = @(x) radiusTheta(x,b,L);
%     nonlinFunc = @(x) fitDual.AFN.thetaCon( x, q, p, pTri );
%     
%     Aeq = ones(1,size(L,2))/size(L,2);
%     beq = 0;
%     optimset = optimoptions('fmincon','Display','iter','Algorithm','active-set','TolFun',1e-5, ... 
%                             'MaxFunEvals',5e6,'MaxIter',5e2,'GradObj','on','GradConstr','on');
%                         
%     [theta,ERes] = fmincon(energyFunc,L\b,A,c,Aeq,beq,[],[],[],optimset);

%     lambda = theta(end);
%     theta = theta(1:(end-1));
%     q = bsxfun(@times,sqrt(lambda),q);
    
end

% 
% function [E,dE] = radiusTheta(x,b,L)
%     
%     E = sum( .5*(L*x - b).^2 );
%     dE = L'*(L*x-b);
%     
% end

