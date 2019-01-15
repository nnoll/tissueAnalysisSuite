function [ sigma ] = symmetrizedStress( sigma )
    % SYMMETRIZED STRESS 
    
    midline = round(size(sigma,1)/2);

    % XX Component
    sigma(1:midline,:,1) = nanmean(cat(3,sigma(1:midline,:,1),sigma(size(sigma,1):-1:midline,:,1)),3);
    sigma(size(sigma,1):-1:midline,:,1) = sigma(1:midline,:,1);
    
    % XY Component
    sigma(1:midline,:,2) = nanmean(cat(3,sigma(1:midline,:,2),-sigma(size(sigma,1):-1:midline,:,2)),3);
    sigma(size(sigma,1):-1:midline,:,2) = -sigma(1:midline,:,2);
    
    % YY Component
    sigma(1:midline,:,3) = nanmean(cat(3,sigma(1:midline,:,3),sigma(size(sigma,1):-1:midline,:,3)),3);
    sigma(size(sigma,1):-1:midline,:,3) = sigma(1:midline,:,3);

end

