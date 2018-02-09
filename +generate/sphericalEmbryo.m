function [ emb ] = sphericalEmbryo( N )
%SPHERICALEMBRYO Summary of this function goes here
%   Detailed explanation goes here

%     [v,Tri] = read_ply('./sphere.ply');
    v = generate.RandSampleSphere(N);
    Tri = convhulln(v);
    [ d0 ] = generate.primalBonds( triangulation(Tri,v) );

    az = atan2(v(:,2),v(:,1));
    el = acos(v(:,3));
    x0 = [az,el];
    Lbar = mean( sqrt( sum( (d0*v).^2 ,2)));
    
    energyFun = @(x) generate.sphericalEnergy(x,Lbar);
    optimset = optimoptions('fmincon','Display','iter','Algorithm','active-set','GradObj','on','MaxIter',1000);
    
    [x,E] = fmincon(energyFun,x0,[],[],[],[],[-pi*ones(length(az),1);zeros(length(az),1)],...
                   [pi*ones(length(az),1);pi*ones(length(az),1)],[],optimset);
    rv = zeros(size(x0,1),3);
    rv(:,1) = cos(x(:,1)).*sin(x(:,2));
    rv(:,2) = sin(x(:,1)).*sin(x(:,2));
    rv(:,3) = cos(x(:,2));
    Tri = convhulln(rv);
    
    emb = generate.embryo(rv,Tri);
end

