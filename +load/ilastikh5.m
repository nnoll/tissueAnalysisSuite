function [ mem ] = ilastikh5( Folder, channel)
    % READ_ILASTIKH5 
    % Processes and concatenates all h5 files in Folder, 
    % 
    % Parameters
    % ----------
    % Folder : str
    %   The path containing hdf5 files to process, may include wildcard
    %   specifier (example: '/path/to/data/*.h5' or './path/to/data/')
    % channel : int
    %   the index of the slice of the h5 file containing the prediction
    %   1 is foreground background for the two classes, whereas 2 assumes
    %   background, foreground ordering
    %
    % Returns
    % -------
    % mem : Y x X x T float matrix
    %   the image sequence (in time) of the projected "membrane" data
    %
    % Written Nick Noll 2017-2019, Annotated NPMitchell 2019
    
    if ~exist('in2','var')
        % second parameter (channel) does not exist, so default is 1
        channel = 1 ;
    end
    Directory = dir(Folder);
    mem = [];
    for t = 1:length(Directory)
        fName = Directory(t).name;
        if contains(fName,'.h5')
            h5fn = [Folder,fName] ;
            disp(['Reading ' h5fn]) 
            pred = h5read(h5fn, '/exported_data');
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

