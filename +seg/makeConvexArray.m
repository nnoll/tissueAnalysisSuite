function [ Struct ] = makeConvexArray( Struct, extCell )
% MAKECONVEXARRAY Move vertices to prevent >180 degree angles
% If concavity in the cells (ie there is an angle > 180) and trying to fit 
% with only positive tensions (as is inferred), there will be a large 
% localized error in the dual (sent to infinity) in order to fit a 
% concave angle. Finds vertices > 180 deg, move vertex in direction
% opposite of angle until < 180 deg. 
% 
%

    if (nargin == 1)
        extCell = 1;
    end

    for t = 1:length(Struct) % Iterate over all time
       three_fold = [Struct(t).Cdat.all_threefold];
       if (isfield(Struct(t).Cdat,'hole'))
           holes = [Struct(t).Cdat.hole];
           if (extCell == 1)
               boundaryCells = [1,find(three_fold==0),find(holes==1)]; 
           else
               boundaryCells = [find(three_fold==0),find(holes==1)]; 
           end
       else
           if (extCell == 1)
               boundaryCells = [1,find(three_fold==0)]; 
           else
               boundaryCells = [find(three_fold==0)]; 
           end
       end
       boundaryVerts = Struct(t).Cdat(boundaryCells);
       boundaryVerts = unique([boundaryVerts(:).nverts]);
       
       for v = 1:length(Struct(t).Vdat) % Iterate over all vertices

           if (~ismember(v,boundaryVerts) && length(Struct(t).Vdat(v).nverts) == 3) % If it is a bulk vertex
               rv = double([Struct(t).Vdat(v).vertxcoord;Struct(t).Vdat(v).vertycoord]);
               nverts = Struct(t).Vdat(v).nverts;
               
               n(:,1) = double([Struct(t).Vdat(nverts(1)).vertxcoord;Struct(t).Vdat(nverts(1)).vertycoord]);
               n(:,2) = double([Struct(t).Vdat(nverts(2)).vertxcoord;Struct(t).Vdat(nverts(2)).vertycoord]);
               n(:,3) = double([Struct(t).Vdat(nverts(3)).vertxcoord;Struct(t).Vdat(nverts(3)).vertycoord]);
               
               N = bsxfun(@minus,n,mean(n,2));
               theta = mod(atan2(N(2,:),N(1,:)),2*pi);
               [~,ind] = sort(theta);
               n = n(:,ind);

               r(:,1) = n(:,1) - rv; r(:,1) = r(:,1)/sqrt(sum(r(:,1).^2));
               r(:,2) = n(:,2) - rv; r(:,2) = r(:,2)/sqrt(sum(r(:,2).^2));
               r(:,3) = n(:,3) - rv; r(:,3) = r(:,3)/sqrt(sum(r(:,3).^2));
               
               z12 = cross([r(:,1);0],[r(:,2);0]);
               z23 = cross([r(:,2);0],[r(:,3);0]);
               z31 = cross([r(:,3);0],[r(:,1);0]);
               theta12 = mod(atan2(z12(3),dot(r(:,1),r(:,2))),2*pi);
               theta23 = mod(atan2(z23(3),dot(r(:,2),r(:,3))),2*pi);
               theta31 = mod(atan2(z31(3),dot(r(:,3),r(:,1))),2*pi);

               if (any([theta12,theta23,theta31] > pi))
                    
                    if (theta12 > pi)
                        deltaR = ((n(:,1)-rv)'*[0,-1;1,0]*(n(:,1)-n(:,2))) / ((n(:,3)-rv)'*[0,-1;1,0]*(n(:,1)-n(:,2)));
                        nrv = rv + 1.5*deltaR*(n(:,3) - rv);
                    elseif (theta23 > pi)
                        deltaR = ((n(:,2)-rv)'*[0,-1;1,0]*(n(:,2)-n(:,3))) / ((n(:,1)-rv)'*[0,-1;1,0]*(n(:,2)-n(:,3)));
                        nrv = rv + 1.5*deltaR*(n(:,1) - rv);
                    else
                        deltaR = ((n(:,3)-rv)'*[0,-1;1,0]*(n(:,3)-n(:,1))) / ((n(:,2)-rv)'*[0,-1;1,0]*(n(:,3)-n(:,1)));
                        nrv = rv + 1.5*deltaR*(n(:,2) - rv);
                    end

                    Struct(t).Vdat(v).vertxcoord = double(nrv(1));
                    Struct(t).Vdat(v).vertycoord = double(nrv(2));

               end
               
               if (any([theta12,theta23,theta31] == pi))
                
                    if (theta12 == pi)
                         nrv = rv + .5*r(:,3);
                    elseif (theta23 == pi)
                         nrv = rv + .5*r(:,1);
                    else
                         nrv = rv + .5*r(:,2);
                    end

                    Struct(t).Vdat(v).vertxcoord = double(nrv(1));
                    Struct(t).Vdat(v).vertycoord = double(nrv(2));

               end
               
           end
       end
       
    end

end

