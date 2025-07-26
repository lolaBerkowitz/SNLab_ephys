function plot_paths_for_dataset(dataset_path, save_path, varargin)
%PLOT_PATHS_FOR_DATASET Generate and save behavioral path plots for each session.
%
%   plot_paths_for_dataset(df_path, save_path, exclude_task_string) loads session
%   metadata from the given df_path, filters out sessions containing the specified
%   string(s) in their task_name field, and generates occupancy path plots
%   (by epoch and by trial) for each behavioral session.
%
%   Figures are saved in the specified save_path directory under:
%       - occupancy_by_epoch/
%       - occupancy_by_trial/
%
%   INPUTS:
%       dataset_path             - (char or string) Path to the dataset file
%       save_path           - (char or string) Directory to save the output figures
%       exclude_task_string - (optional, char, string, or cellstr) Task name(s) to exclude 
%                             (e.g., 'home_cage_sleep' or {'taskA','taskB'})
%
%   This function avoids reprocessing if both expected output figures for a
%   session already exist.
%
%   Dependencies: 
%       - SNLab_ephys functions (compile_sessions, plot_path, plot_path_by_trial)
%       - CellExplorer (loadSession)
%

    % Input parser
    p = inputParser;
    addRequired(p, 'dataset_path', @(x) ischar(x) || isstring(x));
    addRequired(p, 'save_path', @(x) ischar(x) || isstring(x));
    addOptional(p, 'exclude_task_string', {'home_cage_sleep'}, ...
        @(x) ischar(x) || isstring(x) || iscellstr(x));
    parse(p, dataset_path, save_path, varargin{:});

    dataset_path = char(p.Results.dataset_path);
    save_path = char(p.Results.save_path);
    exclude_task_string = p.Results.exclude_task_string;

    if ischar(exclude_task_string) || isstring(exclude_task_string)
        exclude_task_string = cellstr(exclude_task_string);
    end

    % Compile sessions
    df = compile_sessions(dataset_path);

    % Ensure output directories exist
    save_path_epoch = fullfile(save_path, 'occupancy_by_epoch');
    if ~exist(save_path_epoch, 'dir')
        mkdir(save_path_epoch);
    end

    save_path_trial = fullfile(save_path, 'occupancy_by_trial');
    if ~exist(save_path_trial, 'dir')
        mkdir(save_path_trial);
    end

    % Filter sessions
    exclude_mask = false(height(df), 1);
    for i = 1:length(exclude_task_string)
        exclude_mask = exclude_mask | contains(df.task_name, exclude_task_string{i});
    end
    behav_sessions = df(~exclude_mask, :);

    % Loop through and plot paths by epoch and trial
    for i = 1:length(behav_sessions.basepath)
        basepath = behav_sessions.basepath{i};
        [~, basename] = fileparts(basepath);
        session = loadSession(basepath, basename);

        % Define expected output file paths
        epoch_file = fullfile(save_path_epoch, [basename, '_behavior_epochs.png']);
        trial_file = fullfile(save_path_trial, [basename, '_behavior_trials.png']);

        % Skip iteration if both figures already exist
        if exist(epoch_file, 'file') && exist(trial_file, 'file')
            continue
        end

        % Load CellExplorer animal behavior file
        load(fullfile(basepath, [basename, '.animal.behavior.mat']));

        % Plot and save epoch figure if missing
        if ~exist(epoch_file, 'file')
            fig = plot_path(behavior, session);
            saveas(fig, epoch_file);
            close(fig)
        end

        % Plot and save trial figure if missing
        if ~exist(trial_file, 'file')
            fig = plot_path_by_trial(behavior);
            saveas(fig, trial_file);
            close(fig)
        end
    end
end
