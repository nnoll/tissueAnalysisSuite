function [ expMyo ] = embryoMyo( time_series, alpha )
%EMBRYOMYO Summary of this function goes here
%   Detailed explanation goes here

    if (nargin == 1)
        alpha = 0;
    end
    
    for t = 1:54
        if (~isnan(time_series(t).myosin_basal(10,10)))
            expMyo{t}(:,:,3) = time_series(t).myosin_2D{1,1} + .5*time_series(t).myosin_trace + alpha*time_series(t).myosin_basal;
            expMyo{t}(:,:,2) = time_series(t).myosin_2D{1,2};
            expMyo{t}(:,:,1) = time_series(t).myosin_2D{2,2} + .5*time_series(t).myosin_trace + alpha*time_series(t).myosin_basal;
        else
            expMyo{t}(:,:,3) = time_series(t).myosin_2D{1,1} + .5*time_series(t).myosin_trace;
            expMyo{t}(:,:,2) = -time_series(t).myosin_2D{1,2};
            expMyo{t}(:,:,1) = time_series(t).myosin_2D{2,2} + .5*time_series(t).myosin_trace;
        end
        
        expMyo{t}(:,:,1) = expMyo{t}(:,(size(expMyo{t},2):-1:1),1);
        expMyo{t}(:,:,2) = expMyo{t}(:,(size(expMyo{t},2):-1:1),2);
        expMyo{t}(:,:,3) = expMyo{t}(:,(size(expMyo{t},2):-1:1),3);
    end

end

