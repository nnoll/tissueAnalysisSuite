function [ E, dE ] = comptEnergy( x, rv0, rvExt, kappa, compt, B1, B2 )
    % COMPT ENERGY 

    Rot = [0,1;-1,0];
    
    r = reshape(x,length(x)/2,2);
    V = length(x)/2;
    
    r0 = [r;rvExt];
    
    rB1 = full(B1*r0);
    rB2 = full(B2*r0);
        
    LB1 = sqrt(sum(rB1.^2,2));
    LB2 = sqrt(sum(rB2.^2,2));
    
    rB1 = bsxfun(@rdivide,rB1,LB1);
    rB2 = bsxfun(@rdivide,rB2,LB2);
    
    ssinTheta = (rB1(:,1).*rB2(:,2)-rB1(:,2).*rB2(:,1));
    S = sign(ssinTheta);
    sinTheta = S.*ssinTheta;
    
    silInd = (sinTheta < 1e-10);
    sinTheta(silInd) = 1e-10;

    chi = compt*log( sinTheta );

    E = kappa*mean(sum((r-rv0).^2,2)) + sum( (chi - mean(chi)).^2 ) / length(chi);
     
    if (nargout == 2)
        
        % Compute gradient of log(sinTheta)
        dSinTheta1 = bsxfun(@rdivide, (rB2*Rot' - bsxfun(@times,ssinTheta,rB1)),ssinTheta.*LB1);
        dSinTheta2 = bsxfun(@rdivide, (-rB1*Rot' - bsxfun(@times,ssinTheta,rB2)),ssinTheta.*LB2);
        dSinTheta1(silInd,:) = 0;
        dSinTheta2(silInd,:) = 0;
        
        dSinTheta_r = zeros(size(B1,1),size(B1,2),2);
        dSinTheta_r(:,:,1) = bsxfun(@times,B1,dSinTheta1(:,1)) + bsxfun(@times,B2,dSinTheta2(:,1));
        dSinTheta_r(:,:,2) = bsxfun(@times,B1,dSinTheta1(:,2)) + bsxfun(@times,B2,dSinTheta2(:,2));
        
        dE = squeeze(sum(dSinTheta_r,1));
        dCompt_r(:,:,1) = compt*dSinTheta_r(:,:,1);
        dCompt_r(:,:,2) = compt*dSinTheta_r(:,:,2);
           
        dE(:,1) = 2/size(compt,1) * sum( bsxfun(@times,chi-mean(chi),dCompt_r(:,:,1) ), 1);
        dE(:,2) = 2/size(compt,1) * sum( bsxfun(@times,chi-mean(chi),dCompt_r(:,:,2) ), 1);

        % Add in effects of springs to initial positions.
        dE = dE(1:V,:);
        dE = dE + (2*kappa*(r-rv0)/V);   
        dE = dE(:);
        
    end
    
end

