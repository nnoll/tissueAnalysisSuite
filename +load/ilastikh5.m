function [ mem ] = ilastikh5( Folder, channel)
    % READ_ILASTIKH5 
    % Processes and concatenates all h5 files in Folder, 
    % 
    % Parameters
    % ----------
    % Folder : str
    %   The path containing hdf5 files to process. 
    % channel : int
    %   the index of the slice of the h5 file containing the prediction
    %   1 is foreground background for the two classes, whereas 2 assumes
    %   background, foreground ordering
    %
    % Returns
    % -------
    % mem : Y x X x T float matrix
    %   the image sequence (in time) of the projected "membrane" data
    
    if ~exist('in2','var')
        % second parameter (channel) does not exist, so default is 1
        channel = 1 ;
    end
    Directory = dir(Folder);
    mem = [];
    for t = 3:length(Directory)
        fName = Directory(t).name;
        if contains(fName,'.h5')
            pred = h5read([Folder,fName], '/exported_data');
            % hdf5 changed format of saved data. Old version below.
            % if (nargin == 1 || mode == 0)
            % else
            % pred = hdf5read([Folder,fName],'volume/prediction');
            % end
            mem = cat(3,mem,squeeze(pred(channel,:,:,:)));            
        end
    end
    
    % Change the dimensionality order --> transpose 
    mem = permute(mem,[2,1,3]);

end

