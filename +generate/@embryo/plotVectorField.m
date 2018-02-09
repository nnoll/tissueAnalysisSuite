function [ ] = plotVectorField( this, VF )

    this.plotPrimal();
    hold on
    quiver3(this.Mesh.X(:,1),this.Mesh.X(:,2),this.Mesh.X(:,3),VF(:,1),VF(:,2),VF(:,3));

    set(gca,'XLim',1.05*[-1,1],'YLim',1.05*[-1,1],'ZLim',1.05*[-1,1]);
    
end

