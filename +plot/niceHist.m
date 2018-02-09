function [ ] = niceHist( data, nbins )
    % NICE HIST 
    
    if (length(data) >=4)
        cmap = jet(length(data));
    elseif (length(data) == 3)
        cmap = [1,0,0;0,1,0;0,0,1];
    elseif (length(data) == 2)
        cmap = [1,0,0;0,0,1];
    else
        cmap = [1,0,0];
    end
    
    for ii = 1:length(data)
        plot.histf(data{ii},nbins,'BarWidth',.75,'FaceColor',cmap(ii,:),'LineWidth',1,'EdgeColor','k','Alpha',.05);
        [~,xi] = hist(data{ii},nbins);
        deltax = xi(2) - xi(1);
        hold on
        [f,xi] = ksdensity(data{ii});
        plot(xi,f*deltax,'LineWidth',3,'Color',cmap(ii,:));
    end

end

