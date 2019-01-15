function [ c, nbonds ] = tensionMyoCorrelation( Struct, smoothSize )
%TENSIONMYOCORRELATION 

    c = zeros(length(Struct),1);
    
    for t = 1:length(Struct)
        [ i1 ] = generate.bondMap( Struct(t) );
        i1 = i1{1};
        i1 = i1(i1>0);
        T = zeros(size(i1));
        m = zeros(size(i1));
        badBonds = [];
        for b = 1:length(i1)
            if (~isempty(Struct(t).Bdat(i1(b)).tension) && ~isempty(Struct(t).Bdat(i1(b)).chem))
                T(b) = Struct(t).Bdat(i1(b)).tension;
                m(b) = Struct(t).Bdat(i1(b)).chem;
            else
                badBonds = [badBonds,b];
            end
        end
        T(badBonds) = [];
        m(badBonds) = [];
        m = double(m);
        
        T = T/mean(T);
        m = m/mean(m);

        if (nargin > 1 && smoothSize > 0)
           rB = zeros(length(i1),2);
           nB = zeros(length(i1),2);
           ii = 1;

           for b = i1'
               rB(ii,:) = .5*double([Struct(t).Vdat(Struct(t).Bdat(b).verts(1)).vertxcoord,Struct(t).Vdat(Struct(t).Bdat(b).verts(1)).vertycoord]...
                          +   [Struct(t).Vdat(Struct(t).Bdat(b).verts(2)).vertxcoord,Struct(t).Vdat(Struct(t).Bdat(b).verts(2)).vertycoord]);
               nB(ii,:) = double([Struct(t).Vdat(Struct(t).Bdat(b).verts(1)).vertxcoord,Struct(t).Vdat(Struct(t).Bdat(b).verts(1)).vertycoord]...
                          -   [Struct(t).Vdat(Struct(t).Bdat(b).verts(2)).vertxcoord,Struct(t).Vdat(Struct(t).Bdat(b).verts(2)).vertycoord]);
               ii = ii + 1;
           end
           nB = bsxfun(@rdivide,nB,sqrt(sum(nB.^2,2)));
           nB(badBonds,:) = [];
           rB(badBonds,:) = [];

           D = pdist2(rB,rB);
           SmK = exp(-D.^2/(2*(smoothSize).^2));
           SmK = bsxfun(@rdivide,SmK,sum(SmK,2));

           M = zeros(length(m),3);
           M(:,1) = SmK*bsxfun(@times,m,nB(:,1).*nB(:,1));
           M(:,2) = SmK*bsxfun(@times,m,nB(:,1).*nB(:,2));
           M(:,3) = SmK*bsxfun(@times,m,nB(:,2).*nB(:,2));

           m = nB(:,1).*M(:,1).*nB(:,1) + 2*nB(:,1).*M(:,2).*nB(:,2) + nB(:,2).*M(:,3).*nB(:,2);
        end

%         ind = m > .1 & T < 2;
%         if (sum(ind) >= 1)
        c(t) = corr(m,T);
%         scatter(m,T)
%         pause
%         end

        nbonds(t) = length(i1);
    end

end

