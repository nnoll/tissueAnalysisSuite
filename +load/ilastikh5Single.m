function [ mem ] = ilastikh5Single( h5fn, channel)
    % READ_ILASTIKH5 
    % Processes and concatenates a single h5 files in Folder, leaving the
    % third dimension as length 1.
    % 
    % Parameters
    % ----------
    % h5fn : str
    %   The path of a single hdf5 file to process
    % channel : int
    %   the index of the slice of the h5 file containing the prediction
    %   
    % Returns
    % -------
    % mem : Y x X x 1 float matrix
    %   the image of the projected "membrane" data (from iLastik output of 
    %   requested channel)
    %
    % Nick Noll & NPMitchell 2019
    
    if ~exist('in2','var')
        % second parameter (channel) does not exist, so default is 1
        channel = 1 ;
    end
    disp(['Reading ' h5fn]) 
    pred = h5read(h5fn, '/exported_data');
    mem = [] ;
    mem = cat(3,mem,squeeze(pred(channel,:,:,:)));            

    % Change the dimensionality order --> transpose 
    mem = permute(mem,[2,1,3]);

end

