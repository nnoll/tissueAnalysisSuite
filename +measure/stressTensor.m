function [ Struct ] = stressTensor( Struct, L, mode )
%STRESSTENSOR - One line description of what the function or script performs (H1 line)
%Optional file header info (to give more details about the function than in the H1 line)
%Optional file header info (to give more details about the function than in the H1 line)
%Optional file header info (to give more details about the function than in the H1 line)
%
% Syntax:  [Struct] = function_name(Struct, L, mode)
%
% Inputs:
%    input1 - Description
%    input2 - Description
%    input3 - Description
%
% Outputs:
%    output1 - Description
%    output2 - Description
%
% Example: 
%    Line 1 of example
%    Line 2 of example
%    Line 3 of example
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2

    if (nargin <= 2)
        mode = 1;
    end
    
    for t = 1:length(Struct)
  
        [ ~, bulkCells ] = fitDual.returnGraph( Struct(t), 1 );
        [ dC, dV, ~, iVerts, iCells ] = fitDual.ATN.computeDiffOperators( Struct(t), 1 );
        
        P = zeros(length(iCells),1);
        T = zeros(size(dC,1),1);

        % Extract the pressure values from Struct into an array
        for c = 1:length(iCells)
            P(c) = Struct(t).Cdat(iCells(c)).pressure; 
        end

        bVerts = zeros(length(Struct(t).Bdat),2);
        for b = 1:length(Struct(t).Bdat)
           bVerts(b,1:length(Struct(t).Bdat(b).verts)) = Struct(t).Bdat(b).verts;
        end

        i1 = zeros(length(T),1);
        for b = 1:length(T)
            verts = iVerts(dV(b,:)~=0);
            
            if (length(verts) == 2)
                ind = find(ismember(bVerts,verts,'rows'));
            else
                ind = [];
            end
            
            if (~isempty(ind) && ~isempty(Struct(t).Bdat(ind).tension))
                T(b) = Struct(t).Bdat(ind).tension;
                i1(b) = ind;
            else
                T(b) = 1;
            end
            
        end

        rv = zeros(length(iVerts),2);
        for ii = 1:length(iVerts)
            rv(ii,1) = double(Struct(t).Vdat(iVerts(ii)).vertxcoord);
            rv(ii,2) = double(Struct(t).Vdat(iVerts(ii)).vertycoord);
        end
        
        if (mode == 1)
            sigma = zeros(length(bulkCells),3);
    
            for c = 1:length(bulkCells)
                cVerts = Struct(t).Cdat(bulkCells(c)).nverts;
                for v = 1:length(cVerts)
                   R0 = rv(iVerts==cVerts(v),:);
                   extVerts = Struct(t).Vdat(cVerts(v)).nverts;
                   extVerts = extVerts(~ismember(extVerts,cVerts));

                   bondID = ( abs(dV(:,iVerts==cVerts(v))) & abs(dV(:,iVerts==extVerts)) );
                   rB = rv(iVerts==extVerts,:) - R0;
                   rB = rB / sqrt(sum(rB.^2));
  
                   % 'Bent' Tension contribution
                   if (std(P) > 1e-6 && i1(bondID) > 0 && ~isinf(Struct(t).Bdat(i1(bondID)).radius) ) %&& ~isinf(Struct(t).Bdat(i1(bondID)).radius))
                       if (size(Struct(t).Bdat(i1(bondID)).rBar,2) == 2)
                           rBar = Struct(t).Bdat(i1(bondID)).rBar;
                       else
                           rBar = Struct(t).Bdat(i1(bondID)).rBar';
                       end
                       tB = R0 - rBar;
                       tB = tB / sqrt(sum(tB.^2));
                       tB = [0,-1;1,0] * tB';
                       if (sign(dot(tB,rB)) < 0)
                           tB = -tB;
                       end
                   else
                       tB = rB;
                   end
                   F = T(bondID) * tB;

                   sigma(c,1) = sigma(c,1) + R0(1)*F(1);
                   sigma(c,2) = sigma(c,2) + .5*(R0(2)*F(1)+R0(1)*F(2));
                   sigma(c,3) = sigma(c,3) + R0(2)*F(2);
                end
            end
        else
            
            rb = dV * rv;
            D = sqrt(sum(rb.^2,2));
            D(D==0) = 1;
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
            
        end
        
        if (nargin > 1)
            S = regionprops(L(:,:,t),'Area');
        end
        
        if (mode == 1)
            if (nargin > 1)
                A = [S(bulkCells).Area];
                sigma = bsxfun(@rdivide,sigma,A');
            end
            for c = 1:length(bulkCells)
                Struct(t).Cdat(bulkCells(c)).stress = [sigma(c,1),sigma(c,2);sigma(c,2),sigma(c,3)]; 
            end
        else
            if (nargin > 1)
                A = [S(iCells).Area];
                sigma = bsxfun(@rdivide,sigma,A');
            end
            for c = 1:length(iCells)
                Struct(t).Cdat(iCells(c)).stress = [sigma(c,1),sigma(c,2);sigma(c,2),sigma(c,3)]; 
            end
        end
        
    end
    
end

