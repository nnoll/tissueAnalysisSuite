function [ mem ] = ilastikh5( Folder, mode )
    %READ_ILASTIKH5 
    
    Directory = dir(Folder);
    mem = [];
    for t = 3:length(Directory)
        fName = Directory(t).name;
        if ~isempty(strfind(fName,'.h5'))
            if (nargin == 1 || mode == 0)
                pred = hdf5read([Folder,fName],'exported_data');
            else
                pred = hdf5read([Folder,fName],'volume/prediction');
            end
            mem = cat(3,mem,squeeze(pred(1,:,:,:)));
        end
    end
    
    mem = permute(mem,[2,1,3]);

end

