function df = compile_sessions(data_folder,varargin)
% compile_sessions finds subdirectories in folder and outputs a table with
% each basepath
%
% input: 
%  data_folder: path to parent directory containing sessions or csv to
%  already compiled sessions.
%  save_path (optional): path to save df as csv (default, none)
% ouput: 
%  table of basepaths,basename from data_folder
% 
% LB 2022

% handle inputs
p = inputParser; 
p.addParameter('save_path',[],@ischar)

parse(p,varargin{:});
save_path = p.Results.save_path;

% if input to csv, load table, else create table from folder contents.
if contains(data_folder,'.csv')
    df = readtable(data_folder,'Delimiter','comma');
else
    
    % use dir to find subdirs and keep only directories
    folders = dir(fullfile(data_folder, '**\*.*'));
    folders = folders(~ismember({folders.name},{'.','..'}),:); % remove ., ..
    folders = folders([folders.isdir],:);
    % create table and save basepaths for each subdirectory
    df = table;
    df.basepath = cellfun(@(S) fullfile(unique({folders.folder}), S), {folders.name}, 'Uniform', 0)';
    df.basename = {folders.name}';
end


% save basepaths if save_path is input
if ~isempty(save_path)
    writetable(df,fullfile(save_path,'sessions_list.csv'));
end

end
