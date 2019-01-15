function [ c ] = myoCorrelation( Struct, Tt, i1, smoothSize )
%TENSIONMYOCORRELATION 

    c = zeros(length(Struct),1);
    
    for t = 1:length(Struct)
        T = Tt{t};
        T = T/mean(T);
        
        m = zeros(size(T));
        badBonds = [];
        for b = 1:length(i1{t})
            if (~isempty(Struct(t).Bdat(i1{t}(b)).chem))
                m(b) = Struct(t).Bdat(i1{t}(b)).chem(1);
            else 
                badBonds = [badBonds,b];
            end
        end
        m = double(m);
        m = m/mean(m);
        
        m(badBonds) = [];
        T(badBonds) = [];
        goodBonds = i1{t};
        goodBonds(badBonds) = [];
        if (nargin > 2 && smoothSize > 0)
           rB = zeros(length(goodBonds),2);
           nB = zeros(length(goodBonds),2);
           ii = 1;

           for b = goodBonds
               rB(ii,:) = .5*double([Struct(t).Vdat(Struct(t).Bdat(b).verts(1)).vertxcoord,Struct(t).Vdat(Struct(t).Bdat(b).verts(1)).vertycoord]...
                          +   [Struct(t).Vdat(Struct(t).Bdat(b).verts(2)).vertxcoord,Struct(t).Vdat(Struct(t).Bdat(b).verts(2)).vertycoord]);
               nB(ii,:) = double([Struct(t).Vdat(Struct(t).Bdat(b).verts(1)).vertxcoord,Struct(t).Vdat(Struct(t).Bdat(b).verts(1)).vertycoord]...
                          -   [Struct(t).Vdat(Struct(t).Bdat(b).verts(2)).vertxcoord,Struct(t).Vdat(Struct(t).Bdat(b).verts(2)).vertycoord]);
               ii = ii + 1;
           end
           nB = bsxfun(@rdivide,nB,sqrt(sum(nB.^2,2)));

           D = pdist2(rB,rB);
           SmK = exp(-D.^2/(2*(smoothSize).^2));

           M = zeros(length(m),3);
           M(:,1) = SmK*bsxfun(@times,m,nB(:,1).*nB(:,1));
           M(:,2) = SmK*bsxfun(@times,m,nB(:,1).*nB(:,2));
           M(:,3) = SmK*bsxfun(@times,m,nB(:,2).*nB(:,2));

           m = nB(:,1).*M(:,1).*nB(:,1) + 2*nB(:,1).*M(:,2).*nB(:,2) + nB(:,2).*M(:,3).*nB(:,2);
        end
        
        ind = m > prctile(m,5) & T < 3;
        if (length(ind) >= 2 && sum(~isnan(m)) >= 1)
            c(t) = corr(m(ind),T(ind));
        end
%         scatter(m,T)
%         pause
    end

end

