function [ Struct ] = estimateMaterial( Struct, n )
    % ESTIMATE MATERIAL 

    [ i1 ] = generate.bondMap( Struct );
    for t = 1:length(Struct)
        [ d0, ~, ~, ~, i0 ] = fitDual.ATN.computeDiffOperators( Struct(t), 1 );

        T = zeros(length(i1{t}),1); % Tension
        R = zeros(length(i1{t}),1); % Length
        rB = zeros(length(i1{t}),2); % Bond direction
        pos = zeros(length(i1{t}),2); % Bond position.

        % Store relevant bond data
        for b = 1:length(i1{t})
            T(b) = Struct(t).Bdat(i1{t}(b)).tension;
            R(b) = Struct(t).Bdat(i1{t}(b)).length;
            rB(b,:) = [Struct(t).Vdat(Struct(t).Bdat(i1{t}(b)).verts(1)).vertxcoord,   ...
                       Struct(t).Vdat(Struct(t).Bdat(i1{t}(b)).verts(1)).vertycoord] - ...
                      [Struct(t).Vdat(Struct(t).Bdat(i1{t}(b)).verts(2)).vertxcoord,   ...
                       Struct(t).Vdat(Struct(t).Bdat(i1{t}(b)).verts(2)).vertycoord];
            rB(b,:) = rB(b,:) / sqrt(sum(rB(b,:).^2));
            pos(b,:) = .5*([Struct(t).Vdat(Struct(t).Bdat(i1{t}(b)).verts(1)).vertxcoord,   ...
                       Struct(t).Vdat(Struct(t).Bdat(i1{t}(b)).verts(1)).vertycoord] + ...
                      [Struct(t).Vdat(Struct(t).Bdat(i1{t}(b)).verts(2)).vertxcoord,   ...
                       Struct(t).Vdat(Struct(t).Bdat(i1{t}(b)).verts(2)).vertycoord]);
        end

        % Build constutitive matrix to fit material properties.
        nParams = ( (n+1) * (n+2) / 2 ); % Number of coefficients needed to describe polynomial in 2d of order n
        conRel = zeros(length(i1{t}),6*nParams);
        col = 1;

        for px = 0:n % Build matrix column by column.
            if (px < n)
                for py = 0:(n-px)
                    % Stiffness Tensor
                    conRel(:,col) = R .* rB(:,1) .* rB(:,1) .* pos(:,1).^px .* pos(:,2).^py;
                    conRel(:,col+nParams) = 2 * R.* rB(:,2) .* rB(:,1) .* pos(:,1).^px .* pos(:,2).^py;
                    conRel(:,col+2*nParams) =  R .* rB(:,2) .* rB(:,2) .* pos(:,1).^px .* pos(:,2).^py;

                    % Intrinsic Tensor
                    conRel(:,col+3*nParams) = -rB(:,1) .* rB(:,1) .* pos(:,1).^px .* pos(:,2).^py;
                    conRel(:,col+4*nParams) = -2 * rB(:,2) .* rB(:,1) .* pos(:,1).^px .* pos(:,2).^py;
                    conRel(:,col+5*nParams) = - rB(:,2) .* rB(:,2) .* pos(:,1).^px .* pos(:,2).^py;
                    col = col + 1;
                end
            else
                % Stiffness Tensor
                conRel(:,col) = R .* rB(:,1) .* rB(:,1) .* pos(:,1).^n;
                conRel(:,col+nParams) = 2 * R.* rB(:,2) .* rB(:,1) .* pos(:,1).^n;
                conRel(:,col+2*nParams) =  R .* rB(:,2) .* rB(:,2) .* pos(:,1).^n;

                % Intrinsic Tensor
                conRel(:,col+3*nParams) = -rB(:,1) .* rB(:,1) .* pos(:,1).^n;
                conRel(:,col+4*nParams) = -2 * rB(:,2) .* rB(:,1) .* pos(:,1).^n;
                conRel(:,col+5*nParams) = - rB(:,2) .* rB(:,2) .* pos(:,1).^n;
            end
        end

        Stiffness = conRel(:,1:3*nParams);
        size(Stiffness)
        N = null(Stiffness);
        nParams
        plot(N)
        pause
        rank(Stiffness)
        kappa =  pinv(Stiffness) * T;
        
        KappaTensor = cell(3,1);
        [ X, Y ] = meshgrid(linspace(1,1024,1024),linspace(1,1024,1024));
        KappaTensor{1} = zeros(size(X));
        KappaTensor{2} = zeros(size(X));
        KappaTensor{3} = zeros(size(X));

        ii = 1;
        for px = 0:n % Build matrix column by column.
            if (px < n)
                for py = 0:(n-px)
                    KappaTensor{1} = KappaTensor{1} + kappa(ii) .* X.^px .* Y.^py;
                    KappaTensor{2} = KappaTensor{2} + kappa(ii+nParams) .* X.^px .* Y.^py;
                    KappaTensor{3} = KappaTensor{3} + kappa(ii+2*nParams) .* X.^px .* Y.^py;
                    ii = ii + 1;
                end
            else
                KappaTensor{1} = KappaTensor{1} + kappa(ii) .* X.^n;
                KappaTensor{2} = KappaTensor{2} + kappa(ii+nParams) .* X.^n;
                KappaTensor{3} = KappaTensor{3} + kappa(ii+2*nParams) .* X.^n;
            end
        end
    end
    
    imagesc(KappaTensor{1})
    pause
    imagesc(KappaTensor{2})
    pause
    imagesc(KappaTensor{3})
    pause
end

