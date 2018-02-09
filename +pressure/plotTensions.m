function [ mov ] = plotTensions( PN, mem )
    % PLOT TENSIONS 
    
    for t = 1:length(PN)
        clf
        imshow(mem(:,:,t))
        hold on
        PN{t}.plotTensionMap();
        mov(t) = getframe();
    end

end

