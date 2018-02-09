function [ mov ] = trackCell( series, L, t0 )
%TRACKCELL 

    delta = length(series);
    Ts = t0 - delta + 1;
    
    mov = zeros(size(L,1),size(L,2),3,delta);
    ii = 1;
    for t = Ts:t0
        mov(:,:,1,ii) = ( L(:,:,t) == 0 )+ ( L(:,:,t) == series(ii) );
        mov(:,:,2,ii) = L(:,:,t) == 0;
        mov(:,:,3,ii) = L(:,:,t) == 0;

        ii = ii + 1;
    end

end

