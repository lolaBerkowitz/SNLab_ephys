function batch_split_dat(metadata_csv_path)
% pulls path and port info from metata csv and runs split_dat.m 
%
% input: 
% metadata_csv_path: excel file with column variables: 
%       data_path: path to intan data to split
%       save_path: path to project folder where data should be saved
%       Port_A: contains animal id recorded at that port
%       Port_B: contains animal id recorded at that port
%       Port_C: contains animal id recorded at that port
%       Port_D: contains animal id recorded at that port
%       split: empty or contains 'done'. 
% ouptput: 
%   runs split_dat.m on data paths that have not been indicated as 'done'. 
%   metadata csv is updated. 
% 
%   LBerkowitz 03/2022

% load csv as table
metadata = readtable(metadata_csv_path);

% loop through paths 

run_idx = find(~contains(metadata.split,'done')); % run files that have not been run

for file = run_idx'
    data_path = metadata.data_path{file};
    save_path = metadata.save_path{file};
    subject_order = {metadata.Port_A{file},metadata.Port_B{file}...
                    ,metadata.Port_C{file},metadata.Port_D{file}};
    
    disp(['Running split dat on data folder :', data_path])
    split_dat(data_path,save_path, subject_order)
    metadata.split{file} = 'done';
end
% Save csv 
writetable(metadata,metadata_csv_path);
end