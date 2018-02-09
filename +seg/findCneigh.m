function [ Struct ] = findCneigh( Struct )

for ii=1:length(Struct)
    for jj=1:length(Struct(ii).Cdat)
        list1=setdiff([unique([Struct(ii).Vdat([Struct(ii).Cdat(jj).nverts]).ncells])],jj);
        
        Struct(ii).Cdat(jj).ncells = list1; %Makes sure list of neighboring cells is unique and non-self containing.
        clear list1;
        
        if (jj > 1) %Not external cell.
            if (~isempty(Struct(ii).Cdat(jj).nverts))
                list2 = [Struct(ii).Cdat(jj).nverts];

                x = [Struct(ii).Vdat(list2).vertxcoord];
                y = [Struct(ii).Vdat(list2).vertycoord];
                xc = [Struct(ii).Cdat(jj).centroid.coord(1)];
                yc = [Struct(ii).Cdat(jj).centroid.coord(2)];
                xtemp = x-xc;
                ytemp = y-yc;

                thetatemp = atan2(ytemp,xtemp);
                [~,thetatemp2] = sort(thetatemp); %Sort vertices based on angle.

                list3 = list2(thetatemp2);
                L3 = length(list3);

                if (L3 > 1)
                    list3(L3+1) = list3(1);
                    Struct(ii).Cdat(jj).orderednverts = list3;
                    Struct(ii).Cdat(jj).inview = 1;
                    for kk=1:length(list3)-1
                        v1=list3(kk);
                        v2=list3(kk+1);
                        shared_cells = Struct(ii).Vdat(v1).ncells(ismembc(Struct(ii).Vdat(v1).ncells,Struct(ii).Vdat(v2).ncells));
                        diff = shared_cells(~ismember(shared_cells,jj));
                        if (length(diff) == 1)
                            Struct(ii).Cdat(jj).orderedncells(kk) = diff;
                            Struct(ii).Cdat(jj).cell_good = 1;
                        else
                            Struct(ii).Cdat(jj).cell_good = 0;
                        end
%                         if (length(setdiff([intersect([Struct(ii).Vdat(v1).ncells],[Struct(ii).Vdat(v2).ncells])],jj))==1)
%                             Struct(ii).Cdat(jj).orderedncells(kk)=setdiff([intersect([Struct(ii).Vdat(v1).ncells],[Struct(ii).Vdat(v2).ncells])],jj);
%                             Struct(ii).Cdat(jj).cell_good=1;
%                         else
%                             Struct(ii).Cdat(jj).cell_good=0;  
%                         end
                    end
                else
                    Struct(ii).Cdat(jj).orderedncells = Struct(ii).Cdat(jj).ncells;
                end
            else
                Struct(ii).Cdat(jj).inview=0;
            end
        end
    end
    
%     vonelist = Struct(ii).Cdat(1).nverts;
%     conelist = Struct(ii).Cdat(1).ncells;
% 
%     nextc = conelist(1);
% 
%     count = 0;
% 
%     vorderedlist = [];
%     corderedlist = [nextc];
% 
%     while (count <= length(conelist))
%         [vtemplist] = Struct(ii).Cdat(nextc).orderednverts;
%         vorderedlist = [vorderedlist,setdiff(intersect(vtemplist,vonelist),vorderedlist)];
%         nextcc = setdiff([Struct(ii).Vdat(vorderedlist(end)).ncells],[1,nextc]);
%         clear vtemplist  nextc
%         nextc = nextcc;
%         corderedlist = [corderedlist,nextc];
%         clear nextcc
%         count = count+1;
%     end
% 
%     Struct(ii).Cdat(1).orderedncells=corderedlist(1:end-2);
%     Struct(ii).Cdat(1).orderednverts=vorderedlist;

    clear vorderedlist corderedlist count vonelist conelist nextc

end

end
