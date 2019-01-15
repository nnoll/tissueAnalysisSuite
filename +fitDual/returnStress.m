function [ sigma, A ] = returnStress( PN, ERes, r0 )

    sigma = zeros(size(PN,1),size(PN,2),3);
    A = zeros(size(PN,1),size(PN,2));
    for ii = 1:(size(PN,1))
        for jj = 1:(size(PN,2))
            if ( ERes(ii,jj) < 50 )
                [sigma(ii,jj,:),A(ii,jj)] = PN{ii,jj}.computeBoundaryStress(r0{ii,jj});
            end
        end
    end

end

