function [ T, P ] = invertSubMech( Struct, mode, XRange, YRange )
    % INVERT MECH 

    Rot = [0,-1;1,0];
    [ dC, dV, ~, iVerts, iCells ] = fitDual.subATN.computeSubDiffOperators( Struct, XRange, YRange );  

    %Log vertex positions.
    rv = zeros(length(iVerts),2);
    for ii = 1:length(iVerts)
        rv(ii,1) = double(Struct.Vdat(iVerts(ii)).vertxcoord);
        rv(ii,2) = double(Struct.Vdat(iVerts(ii)).vertycoord);
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
    
    if (nargin <= 3 || mode == 1)
        
        mTension = sparse([mTension;ones(1,size(dV,1))]);
        T = mTension \ [zeros(size(mTension,1)-1,1);size(dV,1)];
        T = abs(T);
        P = ones(length(iCells),1);
        
    else
        
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

        mCurv = [I,curv];
        m = [m;mCurv];
        
        bndryCond = zeros(2,size(m,2));
        bndryCond(1,(1:size(mTension,2))) = 1;
        bndryCond(2,(size(mTension,2)+1):size(bndryCond,2)) = 1;
        
        m = sparse([m;bndryCond]);
        b = [zeros(size(m,1)-2,1);size(mTension,2);size(dC,2)];
        
        Mech = m \ b;

        T = Mech(1:size(mTension,2));
        P = Mech((size(mTension,2)+1):end);
        
    end
    
end

