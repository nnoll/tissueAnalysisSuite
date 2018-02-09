function [ stress, xVec, yVec ] = takeContinuumLimit( Struct, xG, yG, sigma )
    % TAKE CONTINUUM LIMIT 
    
    xVec = linspace(xG(1),xG(2),218);
    yVec = linspace(yG(1),yG(2),257);
    [X,Y] = meshgrid(xVec,yVec);
    stress = zeros(size(X,1),size(X,2),3);
    for b = 1:length(Struct.Bdat)
        if (~isempty(Struct.Bdat(b).tension))
            v = Struct.Bdat(b).verts;
            rB = [Struct.Vdat(v(1)).vertxcoord,Struct.Vdat(v(1)).vertycoord] ...
               - [Struct.Vdat(v(2)).vertxcoord,Struct.Vdat(v(2)).vertycoord];
            D = sqrt(sum(rB.^2));
            
            R = .5*([Struct.Vdat(v(1)).vertxcoord,Struct.Vdat(v(1)).vertycoord] ...
                  + [Struct.Vdat(v(2)).vertxcoord,Struct.Vdat(v(2)).vertycoord]);
           
            rB = rB / D;
            nB = rB * [0,-1;1,0];
            
            y = nB(1) * (X-R(1)) + nB(2) * (Y-R(2));
            x = rB(1) * (X-R(1)) + rB(2) * (Y-R(2));
            
            I = (exp(-y.^2/(2*sigma^2))/(2*sqrt(2*pi)*sigma)) .* (erf( (D/2 - x) / sqrt(2*sigma^2) ) - erf( (-D/2 - x) / sqrt(2*sigma^2) ) );
            stress(:,:,1) = stress(:,:,1) + rB(1)*rB(1)*Struct.Bdat(b).tension*I;
            stress(:,:,2) = stress(:,:,2) + rB(2)*rB(1)*Struct.Bdat(b).tension*I;
            stress(:,:,3) = stress(:,:,3) + rB(2)*rB(2)*Struct.Bdat(b).tension*I;
        end
    end
end

