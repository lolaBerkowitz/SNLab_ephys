function trim_dat(basepath,varargin)
% trim_dat removes parts of the dat file that does not contain neural data, 
% either from unplugs or for multianimal recordings where subjects have
% different start times. 
% Dependenies: one signal per channel dat files from intan RHD recording
% system. 
%
% inputs: 
%   basepath: path to SNLab data directory for a given session. 
% output: 
% 
% 
p = inputParser;
addParameter(p,'dat_name','digitalin.dat',@ischar) %snlab rig wiring for events
parse(p,varargin{:});

dat_name = p.Results.dat_name;

% load rhd to obtain recording parameters 
% Load info.rhd for recording parameters
[amplifier_channels, ~, aux_input_channels, ~,...
    ~, ~, frequency_parameters,~ ] = ...
    read_Intan_RHD2000_file_snlab(basepath);

% load amp, digitalin,
amp = load_amplifier(basepath,1:64);
digitalIn = process_digitalin(basepath,fs);

temp_samples = trim_dat_epoch(port,2) - trim_dat_epoch(port,1);
batch{port} = ceil(linspace(trim_dat_epoch(port,1),temp_samples,...
    ceil(temp_samples/frequency_parameters.amplifier_sample_rate/4)));

amp_file = fopen([basepath,filesep,'amplifier_.dat'],'w');
tic
fwrite(amp_file, amp.Data.mapped(:,1:idx), 'uint16');
fclose(amp_file);
toc


temp_samples = trim_dat_epoch(port,2) - trim_dat_epoch(port,1);
batch{port} = ceil(linspace(trim_dat_epoch(port,1),temp_samples,...
    ceil(trim_dat_epoch(port,2)/frequency_parameters.amplifier_sample_rate/4)));
end

function write_amp(amp,basepath,batch,idx)
amp_file = fopen(fullfile(basepath,'amplifier.dat'),'w');
% loop though batches
for ii = 1:length(batch)-1
    disp(['batch ',num2str(batch(ii)+1),' to ',num2str(batch(ii+1)),...
        '   ',num2str(ii),' of ',num2str(length(batch)-1)])
    % write to disk
    fwrite(amp_file,amp.Data.mapped(idx,batch(ii)+1:batch(ii+1)), 'int16');
end


end