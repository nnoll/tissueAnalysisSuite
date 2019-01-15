function [ PN, qF, Struct ] = straightVM( Struct, dim )

    % CURVED VERTEX MODEL 
    for t = 1:length(Struct)
       
        [ d0, dV, ~, iVerts, iCells ] = fitDual.ATN.computeDiffOperators( Struct(t), 1 );  
        [ involvedBonds ] = generate.bondMap(Struct(t));
        involvedBonds = involvedBonds{1};
        
        z = 0;
        badBonds = [];
        for b = 1:length(involvedBonds)
            if (involvedBonds(b) > 0)
                c1 = find(iCells == Struct(t).Bdat(involvedBonds(b)).cells(1));
                c2 = find(iCells == Struct(t).Bdat(involvedBonds(b)).cells(2));
                if (length(Struct(t).Bdat(involvedBonds(b)).pix) > z)
                    z = length(Struct(t).Bdat(involvedBonds(b)).pix);
                end
            else
                badBonds = [badBonds,b];
            end
        end
        involvedBonds(badBonds) = [];
        d0(badBonds,:) = [];
                
        rB = zeros(length(involvedBonds),2,z+2);
        for ii = 1:length(involvedBonds)
            ind = Struct(t).Bdat(involvedBonds(ii)).pix;
            [y,x] = ind2sub(dim,ind);
            x = double(x);
            y = double(y);
            v1 = Struct(t).Bdat(involvedBonds(ii)).verts(1);
            v2 = Struct(t).Bdat(involvedBonds(ii)).verts(2);
            rB(ii,1,1:(length(x)+2)) = [double(x);Struct(t).Vdat(v1).vertxcoord;Struct(t).Vdat(v2).vertxcoord];
            rB(ii,2,1:(length(x)+2)) = [double(y);Struct(t).Vdat(v1).vertycoord;Struct(t).Vdat(v2).vertycoord];
        end
        
        mask = permute(rB(:,1,:) ~= 0,[1,3,2]);
        [ q, theta, tri, ~, ERes ] = fitDual.ATN.returnDual( Struct(t), 1 );
%         [ q, theta, tri ] = fitDual.ATN.obtainStartPoint( Struct(t), 1 );
        if (ERes > 100)
            theta = zeros(size(theta));
        end
        qvec = d0*q;
        QL = sqrt(sum(qvec.^2,2));
        qvec = bsxfun(@rdivide,qvec,QL);
        rho = .5*abs(d0)*q + .5*bsxfun(@times,(d0*theta)./QL,qvec);
        for b = 1:size(d0,1)
            clf
            plot.skel(Struct(t),'k',0)
            hold on
            x = squeeze(rB(b,1,:));
            y = squeeze(rB(b,2,:));
            scatter(x(:),y(:),'g','filled')
            scatter(rho(b,1),rho(b,2),'b','filled')
            PN = pressure.net(q,theta,ones(size(theta)),tri);
            PN.plotPrimal('r')
            pause(.5)
        end
        q0 = [q,theta];

%         Aeq = [[zeros(1,2*size(q,1)),ones(1,size(q,1))/size(q,1)]; ...
%                [ones(1,size(q,1))/size(q,1),zeros(1,2*size(q,1))] ; ...
%                [zeros(1,size(q,1)),ones(1,size(q,1))/size(q,1),zeros(1,size(q,1))]];
%         beq = [mean(q0(:,3));mean(q0(:,1));mean(q0(:,2))];
        Aeq = [zeros(1,2*size(q,1)),ones(1,size(q,1))/size(q,1)]; 
        beq = mean(q0(:,3));
        optimset = optimoptions('fmincon','Display','iter','MaxFunEvals',1e6,'MaxIter',5e3,'TolFun',1e-6,...
                   'Algorithm','sqp','GradObj','on','GradConstr','on','DerivativeCheck','off');
        qF = fmincon(@(q) generate.segEnergy(q, d0, rB, mask),q0,[],[],Aeq,beq,[],[],[],optimset);
        PN = pressure.net(qF(:,1:2),qF(:,3),ones(size(qF(:,3))),tri,iCells);

    end

end

function [ c, ceq, dc, dceq ] = ensurePositive(Q,bCells,d0) 

    q = Q(:,1:2);
    theta = Q(:,3);
    p = Q(:,4);
    
    ceq = [];
    c =  ((d0*p) .* (d0*theta)) - p(bCells(:,1)).*p(bCells(:,2)).*sum((d0*q).^2,2);
    
    if (nargout > 2)
        dP = d0*p;
        dT = d0*theta;
        dQ = d0*q;
        QL = sum(dQ.^2,2);
        
        % X gradient
        gX = -2*(p(bCells(:,1)).*p(bCells(:,2)).*dQ(:,1));
        gX = bsxfun(@times,gX,d0);
        
        % Y gradient
        gY = -2*(p(bCells(:,1)).*p(bCells(:,2)).*dQ(:,2));
        gY = bsxfun(@times,gY,d0);

        % T gradient
        gTh = bsxfun(@times,dP,d0);
        
        % P gradient
        gP = bsxfun(@times,dT,d0);
        gP = gP - bsxfun(@rdivide,bsxfun(@times,(QL.*p(bCells(:,1)).*p(bCells(:,2))),abs(d0))',p)';
        
        dc = [gX,gY,gTh,gP]';
        dceq = [];
    end
end

