function [ timeSeries ] = oneCell( L, label, t0 )
    % ONE CELL 
    
    PixelList = L(:,:,t0) == label;
    timeSeries = zeros(t0-1,1);
    for t = (t0-1):-1:1
        Lt = L(:,:,t);
        S = regionprops(Lt,'Area');
        involvedCells = Lt(PixelList);
        involvedCells = involvedCells(involvedCells>1);
        guess = mode(involvedCells);
        if (guess > 0 && sum(involvedCells==guess) >= .75*S(guess).Area)
            timeSeries(t) = guess;
            PixelList = Lt == guess;
        else
            break
        end
    end
    timeSeries(timeSeries==0) = [];
    timeSeries = [timeSeries;label];
    
end

