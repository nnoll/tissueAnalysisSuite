function [ ] = measureVelocity( this, dat, n, time_series )
%MEASUREVELOCITY 

    [ v, xM, yM ] = continuum.PIV( dat, n, n );
    [ ~, ~, v ] = continuum.refineGrid( xM, yM, v );
    
    Z = time_series(1).charts{1};
    Phi = time_series(1).charts{2};

    Z = (Z - min(Z(:)))/(max(Z(:))-min(Z(:)));
    Phi = 2*pi*(((Phi - min(Phi(:)))/(max(Phi(:)) - min(Phi(:))))) - pi;
    
    [z,ph] = this.emb.surf.computeSurfaceCoords( this.emb.Mesh.X );
    ph = -ph;
    ph = mod(ph+pi+2.892,2*pi) - pi;
    z = (z - min(z))/(max(z)-min(z));
    
    for t = 1:length(v)
       % Interpolate onto mesh
        this.v{t} = bsxfun(@times,interp2(Z,Phi,v{t}(:,:,1),z,ph),this.emb.e1) + bsxfun(@times,interp2(Z,Phi,v{t}(:,:,2),z,ph),this.emb.e2);        
    end
end

