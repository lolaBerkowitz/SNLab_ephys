function batch_analysis(csv_file_path, function_name)
% Reads in a CSV file containing basepaths, and a function name, and runs each
% basepath in the function.
%
% Inputs:
% - csv_file_path: The path to the CSV file containing basepaths.
% - function_name: The name of the function to run each basepath in.

% Read the CSV file.
basepaths_table = readtable(csv_file_path);

% Convert the basepaths to a cell array.
basepaths = cellstr(basepaths_table.basepath);

% Loop over each basepath and run the specified function.
for i = 1:length(basepaths)
    feval(function_name, basepaths{i});
end
