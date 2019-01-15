function [ mov ] = plotEstPrimal( PN, mem )
    %PLOT EST PRIMAL
    
    for t = 1:length(PN)
        clf
        imshow(mem(:,:,t))
        hold all
        PN{t}.plotPrimal();
        mov(t) = getframe();
    end
    
end

