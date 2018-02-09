function [ q, p ] = estimateSeed( Struct, extCell )
    % ESTIMATE SEED 
    
    [ ~, bulk0, ext0, ~, rC ] = fitDual.returnGraph( Struct, extCell );
    i0 = [bulk0,ext0];
    nBulk = length(bulk0);
    [ d0, rho ] = fitDual.AFN.returnBondCentroid( Struct, i0 );
    
    L = [d0,zeros(size(d0)),-bsxfun(@times,rho(:,1),d0); ...
         zeros(size(d0)),d0,-bsxfun(@times,rho(:,2),d0)];
    b = zeros(size(L,1),1);
    
    bConX = zeros(length(ext0),size(L,2));
    bConX((1:length(ext0))+length(ext0)*((nBulk+1):length(i0))) = 1; 
    bConY = zeros(length(ext0),size(L,2));
    bConY((1:length(ext0))+length(ext0)*((length(i0)+nBulk+1):2*length(i0))) = 1; 
    L = [L;bConX(1,:);bConY(1,:)];
    conX = rC((nBulk+1):length(i0),1) - mean(rC((nBulk+1):length(i0),1));
    conY = rC((nBulk+1):length(i0),2) - mean(rC((nBulk+1):length(i0),2));
    b = [b;conX(1);conY(1)];
    
    
    % Outside q's at cell centroids. 
    nB = size(L,1);
    nC = size(d0,2);
    L = [L; [ones(1,size(d0,2))/nB,zeros(1,size(d0,2)),zeros(1,size(d0,2))]; ...
            [zeros(1,size(d0,2)),ones(1,size(d0,2))/nB,zeros(1,size(d0,2))]; ...
            [zeros(1,size(d0,2)),zeros(1,size(d0,2)),ones(1,size(d0,2))/nC]];
    b = [b;0;0;1];
    
    x = sparse(L) \ b;

    x = reshape(x,length(x)/3,3);
    q = x(:,1:2);
    p = x(:,3);
    scale = mean(p);
    p = p/scale;
    q = bsxfun(@rdivide,q,p);
    
%     rhoEst = bsxfun(@rdivide,d0*q,d0*p);
%     scatter(rho(:,1),rho(:,2),'b')
%     hold on
%     scatter(rhoEst(:,1),rhoEst(:,2),'r')
%     quiver(rho(:,1),rho(:,2),rhoEst(:,1)-rho(:,1),rhoEst(:,2)-rho(:,2),0)
%     pause

    
end

