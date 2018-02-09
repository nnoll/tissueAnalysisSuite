function [ c, nbonds ] = tensionMyoCorrelation( Struct, PN, smoothSize )
%TENSIONMYOCORRELATION 

    c = zeros(length(Struct),1);
    if ( nargin == 1 )
        [ Struct2, i1 ] = MI.storeMech( Struct, 1, 1 );
        for t = 1:length(c)
            T = [Struct2(t).Bdat(i1{t}).oldTension];
            m = [Struct2(t).Bdat(i1{t}).chem];
%             m = m(1:2:end);
            nbonds(t) = length(find(i1{t}));
            if (length(T) >= length(m))
                c(t) = corr(m',T(1:length(m))');
            end
        end
    else
        for t = 1:length(Struct)
            [ Struct2, i1 ] = PN{t}.uploadMechanics(Struct(t));
            T = [Struct2.Bdat(i1).tension];
            T = T/mean(T);
            m = [Struct2.Bdat(i1).chem];
            m = double(m);
            m = m/mean(m);
            
            if (nargin > 2 && smoothSize > 0)
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
               
               D = pdist2(rB,rB);
               SmK = exp(-D.^2/(2*(smoothSize).^2));
               SmK = bsxfun(@rdivide,SmK,sum(SmK,2));
%                m = SmK*bsxfun(@times,m',nB); m = dot(nB,m,2)';
%                T = SmK*bsxfun(@times,T',nB); T = dot(nB,T,2)';
               M = zeros(length(m),3);
               M(:,1) = SmK*bsxfun(@times,m',nB(:,1).*nB(:,1));
               M(:,2) = SmK*bsxfun(@times,m',nB(:,1).*nB(:,2));
               M(:,3) = SmK*bsxfun(@times,m',nB(:,2).*nB(:,2));
               
               S = zeros(length(T),3);
               S(:,1) = SmK*bsxfun(@times,T',nB(:,1).*nB(:,1));
               S(:,2) = SmK*bsxfun(@times,T',nB(:,1).*nB(:,2));
               S(:,3) = SmK*bsxfun(@times,T',nB(:,2).*nB(:,2));
               
               m = nB(:,1).*M(:,1).*nB(:,1) + 2*nB(:,1).*M(:,2).*nB(:,2) + nB(:,2).*M(:,3).*nB(:,2);
               m = m';
%                T = nB(:,1).*S(:,1).*nB(:,1) + 2*nB(:,1).*S(:,2).*nB(:,2) + nB(:,2).*S(:,3).*nB(:,2);
%                T = T';
            end
            
            ind = m > prctile(m,2) & T < 2.5;

            if (sum(ind) >= 1)
                c(t) = corr(m(ind)',T(ind)');
            end
%             clf
%             scatter(m,T,'bo')
%             hold on
%             scatter(m(ind),T(ind),'r.')
%             pause
            nbonds(t) = length(i1);
        end
    end

end

