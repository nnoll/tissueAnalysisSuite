function [ bStruct ] = get_bCells( bStruct )
%GET_BCELLS

for t = 1:length(bStruct)
    for b = 1:length(bStruct(t).Bdat)
        bverts = bStruct(t).Bdat(b).verts;
        if (sum(bverts) > 0)
            bStruct(t).Bdat(b).cells = bStruct(t).Vdat(bverts(1)).ncells(ismember(bStruct(t).Vdat(bverts(1)).ncells,bStruct(t).Vdat(bverts(2)).ncells))';
%             bStruct(t).Bdat(b).cells(bStruct(t).Bdat(b).cells==1) = [];
        else
            bStruct(t).Bdat(b).cells = zeros(2,1);
        end
    end
end

end

