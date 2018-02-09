function [ extVerts ] = returnExtVerts( Struct, bulkVerts )
%RETURNEXTVERTS Summary of this function goes here
%   Detailed explanation goes here

    involvedVerts = [];
    for ii = 1:length(bulkVerts)
       involvedVerts = [involvedVerts,Struct.Vdat(bulkVerts(ii)).nverts]; 
    end
    
    extVerts = unique(involvedVerts(~ismember(involvedVerts,bulkVerts)));

end

