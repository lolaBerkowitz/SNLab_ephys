function aux = load_auxiliary(basepath,n_channels)
% load intan auxiliary file mapped 
% Load the file
contFile = fullfile(basepath,'auxiliary.dat');
file = dir(contFile);
samples = file.bytes/(n_channels * 2); %int16 = 2 bytes
aux = memmapfile(contFile,'Format',{'uint16' [n_channels samples] 'mapped'});

end