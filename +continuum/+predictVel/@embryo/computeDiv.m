function [ f ] = computeDiv( this, sigma, Z, Phi )
%COMPUTEDIV 

    Z = (Z - min(Z(:)))/(max(Z(:))-min(Z(:)));
    Phi = 2*pi*(((Phi - min(Phi(:)))/(max(Phi(:)) - min(Phi(:))))) - pi;
    
    [z,ph] = this.emb.surf.computeSurfaceCoords( this.emb.Mesh.X );
    ph = -ph;
    ph = mod(ph+pi+2.892,2*pi) - pi;
    z = (z - min(z))/(max(z)-min(z));
    
    % Interpolate onto mesh
    for ii = 1:3
        sigmaTmp = sigma(:,:,ii);
        sigmaTmp(isnan(sigmaTmp)) = 0;
        T(:,ii) = interp2(Z,Phi,sigmaTmp,z,ph);
    end
    
    e1 = this.emb.e2;
    e2 = this.emb.e1;

    Temb = zeros(size(T,1),3,3);
    for ii = 1:3
        for jj = 1:3
            Temb(:,ii,jj) = (e1(:,ii) .* T(:,1) .* e1(:,jj)) + 2 * (e1(:,ii) .* T(:,2) .* e2(:,jj)) + (e2(:,ii) .* T(:,3) .* e2(:,jj));
        end
    end
    
    [ f ] = continuum.predictVel.calculateDivTensor( this.emb.Mesh, Temb );
    [ f ] = continuum.predictVel.interpolateFaceOntoVerts(  this.emb.Mesh, f' );

end

