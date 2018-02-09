function [ Struct ] = stressTensor( Struct, L )

    for t = 1:length(Struct)
        [ dC, dV, ~, iVerts, iCells ] = fitDual.ATN.computeDiffOperators( Struct(t), 1 );  
        P = zeros(length(iCells),1);
        T = zeros(size(dC,1),1);

        for c = 1:length(iCells)
            P(c) = Struct(t).Cdat(iCells(c)).pressure; 
        end

        bVerts = zeros(length(Struct(t).Bdat),2);
        for b = 1:length(Struct(t).Bdat)
           bVerts(b,1:length(Struct(t).Bdat(b).verts)) = Struct(t).Bdat(b).verts;
        end

        for b = 1:length(T)
            verts = iVerts(dV(b,:)~=0);
            
            if (length(verts) == 2)
                ind = find(ismember(bVerts,verts,'rows'));
            else
                ind = [];
            end
            
            if (~isempty(ind) && ~isempty(Struct(t).Bdat(ind).tension))
                T(b) = Struct(t).Bdat(ind).tension;
            else
                T(b) = 1;
            end
            
        end

        rv = zeros(length(iVerts),2);
        for ii = 1:length(iVerts)
            rv(ii,1) = double(Struct(t).Vdat(iVerts(ii)).vertxcoord);
            rv(ii,2) = double(Struct(t).Vdat(iVerts(ii)).vertycoord);
        end

        rb = dV * rv;
        D = sqrt(sum(rb.^2,2));
        rb = bsxfun(@rdivide, rb, D);
        Rot = [0,-1;1,0];
        nb = rb*Rot';
        dP = dC * P;

        sigmaB = zeros(size(rb,1),3);
        sigmaB(:,1) = rb(:,1) .* T .* rb(:,1);
        sigmaB(:,2) = rb(:,1) .* T .* rb(:,2);
        sigmaB(:,3) = rb(:,2) .* T .* rb(:,2);

        sigmaP = zeros(size(rb,1),3);
        sigmaP(:,1) = nb(:,1) .* dP.*D .* nb(:,1);
        sigmaP(:,2) = nb(:,1) .* dP.*D .* nb(:,2);
        sigmaP(:,3) = nb(:,2) .* dP.*D .* nb(:,2);

        sigma = sparse(abs(dC)') * sigmaB + .5*dC'*sigmaP;
        S = regionprops(L(:,:,t),'Area');
        A = [S(iCells).Area];
        sigma = bsxfun(@rdivide,sigma,A');

        for c = 1:length(iCells)
            Struct(t).Cdat(iCells(c)).stress = [sigma(c,1),sigma(c,2);sigma(c,2),sigma(c,3)]; 
        end
    end
    
end

