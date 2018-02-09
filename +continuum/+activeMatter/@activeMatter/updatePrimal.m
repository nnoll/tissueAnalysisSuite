function [ new, T1 ] = updatePrimal( this, u )
    % UPDATE PRIMAL. Takes in the current primal, as well as the
    % displacements, and computes the new (displaced) primal.
    
    if (nargin == 1)
        u = zeros(size(this.Mesh.X));
    end
    
    new = continuum.activeMatter.activeMatter( this.Mesh.X + u, this.Mesh.Triangulation, this.surf );
    [ Tri, t1Edges ] = continuum.activeMatter.returnT1( new.Mesh, new.d1 );
    
    if (~isempty(t1Edges))
        T1 = 1;
    else
        T1 = 0;
    end
    
    while (~isempty(t1Edges))
       new = continuum.activeMatter.activeMatter( new.Mesh.X, Tri, new.surf ); 
       [ Tri, t1Edges ] = continuum.activeMatter.returnT1( new.Mesh, new.d1 );
%        if ((length(t1Edges) == 1 && t1Edges == 12)|| (length(t1Edges) == 2 && all(t1Edges == [18,19])))
%            new.plotPrimal();
%            hold all
%            verts1 = find(new.d0(12,:));
%            verts2 = find(new.d0(18,:));
%            verts3 = find(new.d0(19,:));
%            scatter3(new.Mesh.X(verts1,1),new.Mesh.X(verts1,2),new.Mesh.X(verts1,3),'b')
%            scatter3(new.Mesh.X(verts2,1),new.Mesh.X(verts2,2),new.Mesh.X(verts2,3),'g')
%            scatter3(new.Mesh.X(verts3,1),new.Mesh.X(verts3,2),new.Mesh.X(verts3,3),'r')
%            pause
%        end
    end
    
end

