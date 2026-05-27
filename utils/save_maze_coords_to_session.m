function save_maze_coords_to_session(basepath, overwrite)
% SAVE_MAZE_COORDS_TO_BEHAVIOR
%
% Loads maze coordinate CSV files from basepath and inserts them into
% session.behavioralTracking{file}.maze_coords.
%
% Assumes maze coord files already exist in basepath with naming:
%   <video_name>_maze_coords.csv
%
% Example:
%   save_maze_coords_to_behavior(pwd)
%   save_maze_coords_to_behavior(basepath, true)
%
% LB + ChatGPT 2026

if nargin < 1 || isempty(basepath)
    basepath = pwd;
end

if nargin < 2
    overwrite = true;
end

basename = basenameFromBasepath(basepath);

% load session
load(fullfile(basepath, [basename '.session.mat']), 'session');

if ~isfield(session, 'behavioralTracking')
    warning('No behavioralTracking field found.');
    return
end

for file = 1:length(session.behavioralTracking)

    bt = session.behavioralTracking{1,file};

    % skip if maze_coords already exist
    if isfield(bt, 'maze_coords') && ~overwrite
        fprintf('maze_coords already exist for entry %d\n', file);
        continue
    end

    % build expected csv filename
    [~, vidname, ~] = fileparts(bt.notes);

    coord_file = fullfile(basepath, ...
        [vidname '_maze_coords.csv']);

    if ~exist(coord_file, 'file')
        warning('Missing maze coord file:\n%s', coord_file);
        continue
    end

    % load coordinates
    coords_table = readtable(coord_file);

    % optionally rescale if raw coords only
    if ~ismember('x_scaled', coords_table.Properties.VariableNames)

        pixel_distance = bt.pixel_distance;
        pixel_dist_reference = bt.pixel_dist_reference;

        coords_table.x_scaled = ...
            coords_table.x .* ...
            (pixel_dist_reference / pixel_distance);

        coords_table.y_scaled = ...
            coords_table.y .* ...
            (pixel_dist_reference / pixel_distance);
    end

    % save into behavioralTracking
    session.behavioralTracking{1,file}.maze_coords = coords_table;

    fprintf('Loaded maze coords for: %s\n', vidname);

end

% save updated session
save(fullfile(basepath, [basename '.session.mat']), 'session');

fprintf('\nFinished updating behavioralTracking.maze_coords\n');

end