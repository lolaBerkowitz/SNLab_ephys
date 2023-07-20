all_files = dir(data_path);
all_files = dir(data_path);

% for a given basepath, keep the most recently edited events file and
% resave as interictal_spikes.events.mat 

basenames = unique({all_files.name}');
basenames(contains(basenames,'.')) = [];

%% save basepaths with matching file 
for idx = 1:length(basenames)
    basepath = fullfile(data_path,basenames{idx});
    % check if basepath contains ieds
    matching_files = dir(fullfile(basepath, '*interictal_spikes*.mat'));
    if length(matching_files) == 1
        missing_file(idx) = false;
        continue
    else
       missing_file(idx) = true;
    end
end

basenames(missing_file);


%%
for idx = 1:length(basenames)
    basepath = fullfile(data_path,basenames{idx});
    % check if basepath contains ieds
    matching_files = dir(fullfile(basepath, '*interictal_spikes*.mat'));
    if length(matching_files) > 1
        % find oldest file 
        [~, old] = min([matching_files.datenum]);
        [~, new] = max([matching_files.datenum]);
        % rename oldest file 
        movefile(fullfile(basepath, matching_files(old).name), fullfile(basepath, 'interictalspikes_old.events.mat'));
        % rename newest file 
        movefile(fullfile(basepath, matching_files(new).name), fullfile(basepath,[basenames{idx}, '.interictal_spikes.events.mat']));
    elseif length(matching_files) == 1
        if ismember({matching_files.name},[basenames{idx}, '.interictal_spikes.events.mat'])
            continue
        end
        % rename the file
        movefile(fullfile(basepath, matching_files.name), fullfile(basepath,[basenames{idx}, '.interictal_spikes.events.mat']));
    elseif isempty(matching_files)
        continue
    end
end