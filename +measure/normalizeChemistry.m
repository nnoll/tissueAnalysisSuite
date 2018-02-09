function [ Struct ] = normalizeChemistry( Struct, id )
%NORMALIZECHEMISTRY Summary of this function goes here
%   Detailed explanation goes here


    if (nargin == 1)
        id = 1;
    end
    
    for t = 1:length(Struct)
        chem = [];
        for b = 1:length(Struct(t).Bdat)
           if (~isempty(Struct(t).Bdat(b).chem(id)) &&  Struct(t).Bdat(b).chem(id) > 0 )
              chem = [chem,Struct(t).Bdat(b).chem(id)];
           end           
        end
        meanChem = mean(chem);
        for b = 1:length(Struct(t).Bdat)
           if (~isempty(Struct(t).Bdat(b).chem(id)) &&  Struct(t).Bdat(b).chem(id) > 0 )
              Struct(t).Bdat(b).chem(id) = Struct(t).Bdat(b).chem(id) / meanChem;
           end           
        end
    end
    
end

