function [Ltemp] = border_cleanup(Ltrack)

Ltemp=Ltrack;

for ii=1:length(Ltemp(1,1,:))
    L=Ltemp(:,:,ii);
    
    L1=imclearborder(L);
    L2=imsubtract(L,L1)>0;
    L3=imdilate(L2,strel('disk',1));
    L4=imerode(L3,strel('disk',1));
    L1=L1+1;L1(L1==1)=0;L1(L4>0)=1;

    L1(1:end,1)=1;
    L1(1:end,end)=1;
    L1(1,1:end)=1;
    L1(end,1:end,1)=1;
   
    Ltemp(:,:,ii)=L1;

    clear L1 L2 L3 L4 
end

end