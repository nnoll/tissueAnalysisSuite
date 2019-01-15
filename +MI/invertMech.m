function [ T, P, dC, dV, iVerts, iCells, S ] = invertMech( Struct, extCell, mode )
    % INVERT MECH 

    Rot = [0,-1;1,0];
    [ dC, dV, ~, iVerts, iCells ] = fitDual.ATN.computeDiffOperators( Struct, extCell );  

    %Log vertex positions.
    rv = zeros(length(iVerts),2);
    for ii = 1:length(iVerts)
        rv(ii,1) = double(Struct.Vdat(iVerts(ii)).vertxcoord);
        rv(ii,2) = double(Struct.Vdat(iVerts(ii)).vertycoord);
    end
    
    rc = zeros(length(iCells),2);
    for ii = 1:length(iCells)
       rc(ii,1) = double(Struct.Cdat(iCells(ii)).centroid.coord(1)); 
       rc(ii,2) = double(Struct.Cdat(iCells(ii)).centroid.coord(2));
    end
    
    %Obtain bond vectors.
    rb = dV * rv;
    D = sqrt(sum(rb.^2,2));
    nb = rb*Rot;
    rb = bsxfun(@rdivide, rb, D);
    
    adj = dV'*dV;
    Z = diag(adj);
    bulkVerts = Z > 1;
    
    % Build MI matrix
    mTension = [bsxfun(@times,dV(:,bulkVerts),rb(:,1))';bsxfun(@times,dV(:,bulkVerts),rb(:,2))'];
    
    if (nargin <= 2 || mode == 1)
        
        mTension = sparse([mTension;ones(1,size(dV,1))/size(dV,1)]);
        T = mTension \ [zeros(size(mTension,1)-1,1);1];
        T = abs(T);
        P = ones(length(iCells),1);
    elseif (mode == 2)
        
        rC =  zeros(length(iCells),2);
        for ii = 1:length(iCells)
            rC(ii,:) = Struct.Cdat(iCells(ii)).centroid.coord;
        end
        
        nC = dC * rC;
        S = sign(dot(nC,nb,2));

        dC = bsxfun(@times,S,dC);
        mPressure = [bsxfun(@times,abs(dV(:,bulkVerts)),nb(:,1))'*dC;bsxfun(@times,abs(dV(:,bulkVerts)),nb(:,2))'*dC];
        m = [mTension,mPressure];
               
        % Impose curvature condition
        curv = dC;
        badBonds = [];
        bcells = sort(vertcat(Struct.Bdat.cells),2);
        I = eye(size(dC,1));
        
        for b = 1:size(dC,1)
            cells = [find(dC(b,:)==1),find(dC(b,:)==-1)];
            bond = find( (bcells(:,1) == iCells(cells(1))) .* (bcells(:,2) == iCells(cells(2))) );
            if (~isempty(bond))
                R = Struct.Bdat(bond).radius;
                
                if (R < inf)
                    
                    I(b,:) = I(b,:)/R;
                    rBar = Struct.Bdat(bond).rBar;
                    
                    if (size(rBar,2) == 1)
                        rBar = rBar';
                    end
                    
                    d1 = sum( (rBar - Struct.Cdat(iCells(cells(1))).centroid.coord).^2);
                    d2 = sum( (rBar - Struct.Cdat(iCells(cells(2))).centroid.coord).^2);

                    if (d1 < d2)
                        curv(b,:) = -curv(b,:);
                    else
                        curv(b,:) = curv(b,:);
                    end
                    
                else
                    I(b,:) = 0;
                end
            else
                badBonds = [badBonds,b];
            end
            
        end
        
        curv(badBonds,:) = [];
        I(badBonds,:) = [];

        mCurv = [I,-curv];
        m = [m;mCurv];
        
        bndryCond = zeros(2,size(m,2));
        bndryCond(1,(1:size(mTension,2))) = 1;
        bndryCond(2,(size(mTension,2)+1):size(bndryCond,2)) = 1;
        
        m = ([m;bndryCond]);
        b = [zeros(size(m,1)-2,1);size(mTension,2);size(dC,2)];
        
        Mech = m \ b; %pinv(m) * b;
        
        T = Mech(1:size(mTension,2));
        P = Mech((size(mTension,2)+1):end);
        
    else
        [ dV, dC, t1, t2, YL, flatEdges ] = MI.returnBonds( Struct );

        bulkVerts = sum(abs(dV),1) > 2;
        dVpos = dV(:,bulkVerts) == 1;
        dVneg = dV(:,bulkVerts) == -1;
        
        Fbal = [bsxfun(@times,dVpos,t1(:,1))';bsxfun(@times,dVpos,t1(:,2))'] + [bsxfun(@times,dVneg,t2(:,1))';bsxfun(@times,dVneg,t2(:,2))'];
        Fbal = [Fbal;ones(1,size(Fbal,2))/size(Fbal,2)];
        
%         [ bInd ] = generate.bondMap( Struct );
%         bInd = bInd{1};
%     
%         T = zeros(size(bInd));
%         for b = 1:length(T)
%            T(b) = Struct.Bdat(bInd(b)).actual_tension; 
%         end
%         
%         p = zeros(length(iCells),1);
%         for c = 1:length(iCells)
%            p(c) = Struct.Cdat(iCells(c)).actual_pressure; 
%         end
        
%         Fnet = Fbal*T;
%         Fnet = Fnet(1:end-1);
%   
%         Fnet = reshape(Fnet,length(Fnet)/2,2);
%         plot.curveSkel(Struct,'k');
%         hold on
%         scatter(rv(bulkVerts,1),rv(bulkVerts,2),'r','filled')
%         quiver(rv(bulkVerts,1),rv(bulkVerts,2),Fnet(:,1),Fnet(:,2))
%         pause
        
        b = zeros(size(Fbal,1),1);
        b(end) = 1;
        
        T = sparse(Fbal) \ b;
        Tyl = T; Tyl(flatEdges==1) = 0;
        
        S = sign( dot( (bsxfun(@times,T,t1) + bsxfun(@times,T,t2))*[0,1;-1,0], dV*rv, 2) );
        S = S .* sign(dot(dC*rc,dV*rv*[0,-1;1,0],2));
        
        YL = bsxfun(@times,S,YL);
        
%         Test = -YL*p;
        YL = [-YL;ones(1,size(YL,2))/size(YL,2)];
%         
%         rB = .5 * abs(dV) * rv;
%         plot.curveSkel(Struct,'k')
%         hold on
%         scatter(rB(:,1),rB(:,2),[],T./Test,'filled');
%         rc = [Struct.Cdat.centroid];
%         rc = vertcat(rc.coord);
%         scatter(rc(:,1),rc(:,2),'g','filled')
%         pause
%         
%         scatter(T,Test)
%         pause
        
        Tyl = [Tyl;1];
        P = sparse(YL) \ Tyl;

%         P = ones(length(iCells),1);
%         cdfplot(P)
%         pause
    end
    
end

