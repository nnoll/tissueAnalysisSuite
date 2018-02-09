function [ raw ] = tif( Folder )
%READTIFF Summary of this function goes here
%   Detailed explanation goes here

    Directory = dir(Folder);
    raw = [];
    for t = 3:length(Directory)
        fName = Directory(t).name;
        if ~isempty(strfind(fName,'.tif'))
            pred = load.loadTiff([Folder,fName]);
            raw = cat(3,raw,squeeze(pred));
        end
    end

end

