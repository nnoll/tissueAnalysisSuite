function [Vdat]=find_bonds(Vdat,L)

edgenumber=0;

for ii=1:length(Vdat)
    if (length(Vdat(ii).nverts)==4)
        Vdat(ii).fourfold=1;
    end
    nbd=Vdat(ii).nverts;
    N=length(nbd);

    if (N>2)
        clear edge
        edge(N).pix=[];
        for jj=1:N
            Vdat(ii).nbondspix=edge;
            vertindex=Vdat(ii).nverts(jj);
            common=intersect(Vdat(ii).ncells,Vdat(vertindex).ncells);

            if (length(common)==2)
                 image1 = (L==common(1));
                 image2 = (L==common(2));
            else
                break
            end

            image3=zeros(size(L));
            image3(image1==1)=1;
            image3(image2==1)=1;

    %         x1=Vdat(ii).vertxcoord;
    %         y1=Vdat(ii).vertycoord;
    %         x2=Vdat(vertindex).vertxcoord;
    %         y2=Vdat(vertindex).vertycoord;
        %     idx1 = sub2ind(size(L), y1, x1);
        %     
        %     idx2 = sub2ind(size(L), y2, x2);
            image4 = logical(image3);
            image5 = bwmorph(image4,'bridge');
            image6 = image5-image4;
            edge(jj).pix = find(image6==1)';

            Vdat(ii).nbondspix(jj)=edge(jj);


        %       a=Vdat(ii).nbondspix(jj).pix;
        %                     distance=0;
        %                     [y,x]=ind2sub(size(L),a);
        %                     for mm=1:length(a)-1
        %                         distance=distance+sqrt((y(mm+1)-y(mm))^2+(x(mm+1)-x(mm))^2);
        %                     end

            Vdat(ii).nbondslength(jj)=length(Vdat(ii).nbondspix(jj).pix);

            edgenumber=edgenumber+1;

            Vdat(ii).nbonds(jj)=edgenumber;
        end
    end
end

end






