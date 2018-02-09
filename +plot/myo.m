function [ ] = myo( Struct, smoothSize )
%SKELETONIZE_VERTS Takes the network topology and triangulated position of
%vertices to produce a skeletonized image.

    m = [];
    rB = [];
    nB = [];
    ii = 1;
    for b = 1:length(Struct.Bdat)
        if (~isempty(Struct.Bdat(b).chem))
           m(ii) =Struct.Bdat(b).chem;
           if (nargin  > 1 && smoothSize > 0)
               rB(ii,:) = .5*double([Struct.Vdat(Struct.Bdat(b).verts(1)).vertxcoord,Struct.Vdat(Struct.Bdat(b).verts(1)).vertycoord]...
                                  +   [Struct.Vdat(Struct.Bdat(b).verts(2)).vertxcoord,Struct.Vdat(Struct.Bdat(b).verts(2)).vertycoord]);
               nB(ii,:) = double([Struct.Vdat(Struct.Bdat(b).verts(1)).vertxcoord,Struct.Vdat(Struct.Bdat(b).verts(1)).vertycoord]...
                          -   [Struct.Vdat(Struct.Bdat(b).verts(2)).vertxcoord,Struct.Vdat(Struct.Bdat(b).verts(2)).vertycoord]);
           end
           ii = ii + 1;
        end
    end
    
    if (nargin > 1 && smoothSize > 0)
       nB = bsxfun(@rdivide,nB,sqrt(sum(nB.^2,2)));
               
       D = pdist2(rB,rB);
       SmK = exp(-D.^2/(2*(smoothSize).^2));
       SmK = bsxfun(@rdivide,SmK,sum(SmK,2));

       M = zeros(length(m),3);
       M(:,1) = SmK*bsxfun(@times,m',nB(:,1).*nB(:,1));
       M(:,2) = SmK*bsxfun(@times,m',nB(:,1).*nB(:,2));
       M(:,3) = SmK*bsxfun(@times,m',nB(:,2).*nB(:,2)); 
       m = nB(:,1).*M(:,1).*nB(:,1) + 2*nB(:,1).*M(:,2).*nB(:,2) + nB(:,2).*M(:,3).*nB(:,2);
    end
    
    Tmax = prctile(m,95);
    Tmin = prctile(m,5);

    T = (m-Tmin)/(Tmax-Tmin);
    T(T>1) = 1;
    T(T<0) = 0;
    
    % Convert tension to color.
    cmap = hot(256);
    x = linspace(0,1,256);
    
    Tcolor(:,1) = interp1(x,cmap(:,1),T);
    Tcolor(:,2) = interp1(x,cmap(:,2),T);
    Tcolor(:,3) = interp1(x,cmap(:,3),T);
    
    [ ~, dV, ~, iverts ] = fitDual.ATN.computeDiffOperators( Struct, 1 );  
    Vlist = Struct.Vdat;
    hold on
    for b = 1:size(dV,1)
        v = iverts(dV(b,:)==1);
        nv = iverts(dV(b,:)==-1);
        plot([Vlist(v).vertxcoord,Vlist(nv).vertxcoord],[Vlist(v).vertycoord,Vlist(nv).vertycoord],'Color',Tcolor(b,:),'LineWidth',2);
    end
    
end

