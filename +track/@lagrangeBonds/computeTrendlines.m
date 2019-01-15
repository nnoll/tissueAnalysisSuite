function [ myo, r, angle ] = computeTrendlines( this, Ti, Tf )
%COMPUTETRENDLINE 

    myo = zeros(length(this.times),1);
    r = zeros(length(this.times),1);
    angle = zeros(length(this.times),1);

    for ii = 1:length(myo)
        indF = find(this.times{ii} == Tf);
        if (~isempty(indF))
            indI = indF - (Tf - Ti);
        else
            indF = length(this.times{ii});
            indI = indF - (Tf - Ti);
        end

        p = polyfit(linspace(0,1,indF-indI+1),this.myo{ii}(indI:indF),1);
        myo(ii) = p(1)/mean(this.myo{ii}(indI:indF));
        
        p = polyfit(linspace(0,1,indF-indI+1),this.r{ii}(indI:indF),1);
        r(ii) = p(1)/mean(this.r{ii}(indI:indF));

        p = polyfit(linspace(0,1,indF-indI+1),mod(this.theta{ii}(indI:indF),pi),1);
        angle(ii) = p(1)/pi;
    end
end

