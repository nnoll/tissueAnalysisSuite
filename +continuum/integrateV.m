function [ u ] = integrateV( v, x, y )
%INTEGRATEV 

    u = cell(size(v));
    u{1} = zeros(size(v{1}));
    for t = 1:length(v)
        [ vGrad, xG, yG ] = continuum.computeGradV( x, y, u(t) );
        [ vIGrad ] = continuum.reinterpolateTensor( xG, yG, vGrad(1), x, y );
        deltaU = v{t};
        deltaU(1,:) = deltaU(1,:) + (squeeze(v{t}(1,:)) .* squeeze(vIGrad{1}(1,1,:))') ...
                                  + (squeeze(v{t}(2,:)) .* squeeze(vIGrad{1}(1,2,:))');
        deltaU(2,:) = deltaU(2,:) + (squeeze(v{t}(1,:)) .* squeeze(vIGrad{1}(2,1,:))') ...
                                  + (squeeze(v{t}(2,:)) .* squeeze(vIGrad{1}(2,2,:))');
        u{t+1} = u{t} + deltaU;
    end

end

