function [ ] = tif2h5( Folder, mode )
%READTIFF Summary of this function goes here
%   Detailed explanation goes here

    Directory = dir(Folder);
    raw = [];
    for t = 3:length(Directory)
        fName = Directory(t).name;
        if ~isempty(strfind(fName,'.tif'))
            pred = load.loadTiff([Folder,fName]);
            raw = cat(4,raw,squeeze(pred));
        end
    end
    
    hdf5write([Folder,'allFiles.h5'], '/data', uint8(raw));
    
    if (mode==1)
        for t = 1:size(raw,4)
            hdf5write([Folder,'t',num2str(t),'.h5'], '/data', uint8(raw(:,:,:,t)));
        end
    end

end

