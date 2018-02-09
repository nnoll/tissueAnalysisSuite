function [ Xg, Yg, sigmaG ] = interpolateCellHM( sigma, i0, L )
    % INTERPOLATE CELL HM 

    S = regionprops(L,'Centroid');
    Rc = vertcat(S.Centroid);
    Rc = Rc(i0,:);
    
    for ii = 1:3
        F = scatteredInterpolant(Rc(:,1),Rc(:,2),sigma(:,ii),'natural');
        [Yg,Xg] = meshgrid(linspace(min(Rc(:,2)),max(Rc(:,2)),257),linspace(min(Rc(:,1)),max(Rc(:,1)),218));
        Zg = F(Xg,Yg)';
        if (ii == 1 || ii == 3)
            Zg(1:128,:) = .5*(Zg(1:128,:) + Zg(257:-1:130,:));     
            Zg(257:-1:130,:) = Zg(1:128,:);
        else
            Zg(1:128,:) = .5*(Zg(1:128,:) - Zg(257:-1:130,:)); 
            Zg(257:-1:130,:) = -Zg(1:128,:);
        end
        sigmaG(:,:,ii) = medfilt2(Zg,[50,50],'symmetric');
    end
    
    
end

