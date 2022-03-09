function amp = load_amplifier(basepath,amplifier_channels)
% load intan amplifier file mapped 
basename = basenameFromBasepath(basepath);

% Load the file
n_channels = size(amplifier_channels,2);

if exist([basepath,filesep,'amplifier.dat'],'file')
    contFile = fullfile(basepath,'amplifier.dat');
elseif exist([basepath,filesep,basename,'.dat'],'file')
    contFile = fullfile(basepath,[basename,'.dat']);
else
    disp('no amplifier file found')
    return
end

file = dir(contFile);
samples = file.bytes/(n_channels * 2); %int16 = 2 bytes
amp = memmapfile(contFile,'Format',{'int16' [n_channels, samples] 'mapped'});
end