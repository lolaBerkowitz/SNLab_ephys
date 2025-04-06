

function v = load_analogin(basepath)
% Load info.rhd for recording parameters
[amplifier_channels, ~, aux_input_channels, ~,...
    ~, ~, frequency_parameters,board_adc_channels ] = ...
    read_Intan_RHD2000_file_snlab(basepath);


num_channels = length(board_adc_channels); % ADC input info from header file
filename = fullfile(basepath,'analogin.dat');
fileinfo = dir(fullfile(basepath,'analogin.dat'));
num_samples = fileinfo.bytes/(num_channels * 2); % uint16 = 2 bytes

fid = fopen(fullfile(basepath,'analogin.dat'), 'r');
v = fread(fid, [num_channels, num_samples], 'uint16');
fclose(fid);

% 
end