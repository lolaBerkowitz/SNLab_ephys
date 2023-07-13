function t = load_time(basepath)

% Load info.rhd for recording parameters
[~, ~, ~, ~,...
    ~, ~, frequency_parameters,~ ] = ...
    read_Intan_RHD2000_file_snlab(basepath);

fileinfo = dir(fullfile(basepath,'time.dat'));
num_samples = fileinfo.bytes/4; % int32 = 4 bytes
fid = fopen(fullfile(basepath,'time.dat'), 'r');
t = fread(fid, num_samples, 'int32');
fclose(fid);

% convert to seconds
t = t / frequency_parameters.amplifier_sample_rate;

end