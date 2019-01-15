function [ PN ] = rescaleGridDual( PN )
    % RESCALE GRID DUAL   
    
    for ii = 1:size(PN,1)
        for jj = 1:size(PN,2)
            if (~isempty(PN{ii,jj}))
                Tscale = mean(sqrt(sum((PN{ii,jj}.d0*PN{ii,jj}.q).^2,2)));
                PN{ii,jj} = PN{ii,jj}.rescaleDual(1/Tscale,[0,0]);
            end
        end
    end
    
end

