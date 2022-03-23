function plot_impedance(file,ref_file)
% plots intan impedance measurement from two timepoints. Assumes file is
% stored in folder named after probe. 

[basepath,file_name] = fileparts(file);
[~,ref_name] = fileparts(ref_file);

% basename should be name of probe
basename = basenameFromBasepath(basepath);

df2 = readtable(ref_file); % reference
df = readtable(file); % 

% generate plot
fig = figure;
fig.Color = [1 1 1];
histogram(df2.ImpedanceMagnitudeAt1000Hz_ohms_/1000000,75);
hold on;
histogram(df.ImpedanceMagnitudeAt1000Hz_ohms_/1000000,75);
title(basename)
xlabel('Impedance (MOhm)')
ylabel('Count')
legend({ref_name,file_name})

saveas(fig,fullfile(basepath,[basename,'_explant_histogram.tiff']),'tiff')
end