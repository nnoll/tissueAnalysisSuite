function [ myo ] = myoTensor( time_series )
%MYOTENSOR Summary of this function goes here
%   Detailed explanation goes here

    for t = 1:54
        myo{t}(:,:,1) = time_series(t).myosin_2D{1,1} + .5*time_series(t).myosin_trace;
        myo{t}(:,:,2) = time_series(t).myosin_2D{1,2};
        myo{t}(:,:,3) = time_series(t).myosin_2D{2,2} + .5*time_series(t).myosin_trace;
    end

end

