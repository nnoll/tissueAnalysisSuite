function [ r, Tri, b0 ] = meshGrid( N, alpha )
    %   GENERATEMESHGRID 
    % N should be odd for proper function.
    
    x1 = linspace(-1,1,N);
    deltax = x1(2)-x1(1);
    x2 = linspace(-1+(deltax/2),1-(deltax/2),N-1);
    
    y2 = linspace(-1,1,(N+1)/2);
    deltay = y2(2) - y2(1);
    y1 = linspace(-1+(deltay/2),1-(deltay/2),(N-1)/2);
    
    [X1,Y1] = meshgrid(x1,y1);
    [X2,Y2] = meshgrid(x2,y2);
    
    yBulk = y2(abs(y2) < 1)';
    
    Y1 = sqrt(3)/2 * Y1;
    Y2 = sqrt(3)/2 * Y2;
    yBulk = sqrt(3)/2 * yBulk;
    xBulk = ones(size(yBulk));
    
    r = [X1(:),Y1(:);X2(:),Y2(:)];
    r = vertcat(r,[xBulk,yBulk;-xBulk,yBulk]);

    Tri = delaunay(r(:,1),r(:,2));
    Mesh = MeshTriangle(Tri,r,'Omega');
    cons = Mesh.freeBoundary;
    
    u = alpha*randn(size(r));
    boundaryVerts =  (abs(abs(r(:,1))-1) <= 1e-3) | (abs(abs(r(:,2))-sqrt(3)/2) <= 1e-3);
%     u(boundaryVerts,:) = 0;
    
    [ DT ] = delaunayTriangulation(Mesh.X + u, cons);
    I0 = DT.isInterior();
    Tri = DT.ConnectivityList(I0,:);
    r = Mesh.X + u;
    Mesh = MeshTriangle(Tri,r,'Omega');
    Mesh = Mesh.Remove_Unused_Vertices;
    
    Tri = Mesh.Triangulation;
    r = Mesh.X;
    b0 = (abs(abs(r(:,1))-1) <= 1e-3) | (abs(abs(r(:,2))-sqrt(3)/2) <= 1e-3);
    r = r/deltax;
    
end

