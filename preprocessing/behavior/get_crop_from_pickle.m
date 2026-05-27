function crop_table = get_crop_from_pickle(pickle_file)
%GET_CROP_FROM_PICKLE Load cropping parameters from a DeepLabCut pickle file.
%
% Parameters
% ----------
% pickle_file : char or string
%     Full path to the pickle (.pickle) file containing DeepLabCut
%     metadata and cropping parameters.
%
% Returns
% -------
% crop_table : table
%     Table containing the DeepLabCut model name and cropping
%     parameters with columns:
%
%         model
%         x_min
%         x_max
%         y_min
%         y_max
%
% Notes
% -----
% Requires MATLAB configured with a supported CPython installation.
% The pickle file is loaded using MATLAB's Python interface.
%
% Examples
% --------
% crop_table = get_crop_from_pickle( ...
%     'Y:\data\session1\videoDLC_modelshuffle1.pickle');
%
% disp(crop_table)
% 
% Laura Berkowitz 2026

    % Import Python modules
    pickle = py.importlib.import_module('pickle');
    builtins = py.importlib.import_module('builtins');

    % Open pickle file in binary mode
    fid = builtins.open(pickle_file, 'rb');

    % Load pickle contents
    results = pickle.load(fid);

    % Close file
    fid.close();

    % Extract cropping parameters
    crop_params = double(results{"data"}{"cropping_parameters"});

    % Get filename only
    [~, fname, ~] = fileparts(pickle_file);

    % Extract DLC model name between 'DLC_' and 'shuffle'
    tokens = regexp(fname, 'DLC_(.*?)shuffle', 'tokens');

    if isempty(tokens)
        model = "";
    else
        model = string(tokens{1}{1});
    end

    % Replace zeros with ones
    crop_params(crop_params == 0) = 1;
    
    % Create output table
    crop_table = table( ...
        model, ...
        crop_params(1), ...
        crop_params(2), ...
        crop_params(3), ...
        crop_params(4), ...
        'VariableNames', ...
        {'model', 'x_min', 'x_max', 'y_min', 'y_max'} ...
    );

end