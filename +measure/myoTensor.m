function [ Struct ] = myoTensor( Struct, L )
%MYOTENSOR Summary of this function goes here
%   Detailed explanation goes here

    for t = 1:length(Struct)
        [ dC, dV, ~, iVerts, iCells ] = fitDual.ATN.computeDiffOperators( Struct(t), 1 );  
        M = zeros(size(dC,1),1);

        bVerts = zeros(length(Struct(t).Bdat),2);
        for b = 1:length(Struct(t).Bdat)
           bVerts(b,1:length(Struct(t).Bdat(b).verts)) = Struct(t).Bdat(b).verts;
        end

        for b = 1:length(M)
            verts = iVerts(dV(b,:)~=0);
            
            if (length(verts) == 2)
                ind = find(ismember(bVerts,verts,'rows'));
            else
                ind = [];
            end
            
            if (~isempty(ind) && ~isempty(Struct(t).Bdat(ind).chem))
                M(b) = Struct(t).Bdat(ind).chem;
            else
                M(b) = 1;
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

        sigmaB = zeros(size(rb,1),3);
        sigmaB(:,1) = rb(:,1) .* M .* rb(:,1);
        sigmaB(:,2) = rb(:,1) .* M .* rb(:,2);
        sigmaB(:,3) = rb(:,2) .* M .* rb(:,2);

        sigma = sparse(abs(dC)') * sigmaB;
        S = regionprops(L(:,:,t),'Area');
        A = [S(iCells).Area];
        sigma = bsxfun(@rdivide,sigma,A');

        for c = 1:length(iCells)
            Struct(t).Cdat(iCells(c)).myo = [sigma(c,1),sigma(c,2);sigma(c,2),sigma(c,3)]; 
        end
    end
    
end

