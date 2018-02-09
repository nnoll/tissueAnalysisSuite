function [  ] = plotWeightedTri( this, color )
%PLOTWEIGHTEDTRI Summary of this function goes here
%   Detailed explanation goes here

   if (nargin == 1)
        color = 'b';
   end

    patch('Vertices',this.q,'Faces',this.tri,'FaceColor','none','EdgeColor',color,'LineWidth',2)
    
%     viscircles(this.q,20*sqrt(this.theta),'Color',color,'EnhanceVisibility','off');
    hold on
    [cMap] = plot.brewermap(256,'Blues');
    P = this.p;
    if (std(P) == 0)
        Pcolor = ones(size(P,1),1)*cMap(125,:);
    else
        P = P / mean(P);
        P = (P - min(P)) / (max(P) - min(P));
        xR = linspace(0,1,256);
        Pcolor(:,1) = interp1(xR,cMap(:,1),P);
        Pcolor(:,2) = interp1(xR,cMap(:,2),P);
        Pcolor(:,3) = interp1(xR,cMap(:,3),P);
    end
    
    this.theta = this.theta - min(this.theta);
    R = .5*sqrt(this.theta - min(this.theta)) + .01;

%     scatter(this.q(:,1),this.q(:,2),R,'MarkerFaceColor',color,'MarkerEdgeColor','none')
    for c = 1:size(this.q,1)
        rectangle('Position',[this.q(c,:)-.5*[R(c),R(c)],R(c),R(c)],'Curvature',1,'EdgeColor','b','LineWidth',2,'FaceColor',Pcolor(c,:));
    end
    
end

