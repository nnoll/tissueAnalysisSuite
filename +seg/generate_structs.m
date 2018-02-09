function [Ltemp,Struct] = generate_structs(Ltrack,bond,clear_border,threefold)

    Ltemp=Ltrack;
    
    if (nargin == 3)
        threefold = 0;
    end
    
    if (clear_border)
        for ii=1:length(Ltemp(1,1,:))
            L=Ltemp(:,:,ii);

            L1=imclearborder(L);
            L2=imsubtract(L,L1)>0;
            L3=imdilate(L2,strel('disk',1));
            L4=imerode(L3,strel('disk',1));
            L1=L1+1;L1(L1==1)=0;L1(L4>0)=1;

            L1(1:end,1)=1;L1(1:end,end)=1;L1(1,1:end)=1;L1(end,1:end,1)=1;

            %Remove spurs on exterior edge
            img = L1==0; 
            img2 = bwmorph(img,'spur');
            img3 = imsubtract(img,img2);
            L1(find(img3)) = 1;
            Ltemp(:,:,ii)=L1;

            clear L1 L2 L3 L4 
        end  
    end

    for ii=1:length(Ltemp(1,1,:))
        ii
        L = Ltemp(:,:,ii);
        Struct(ii) = seg.create_Cdat_Vdat_initial(L,bond,threefold);
        for v = 1:length(Struct(ii).Vdat)
            Struct(ii).Vdat(v).vertxcoord = double(Struct(ii).Vdat(v).vertxcoord);
            Struct(ii).Vdat(v).vertycoord = double(Struct(ii).Vdat(v).vertycoord);
        end
    end 

end
    
