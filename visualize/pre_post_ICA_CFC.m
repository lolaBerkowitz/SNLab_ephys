function pre_post_comod = pre_post_ICA_CFC(ica,pyrCh,basepath,save_path,varargin)
% Plot the comodulagram for theta from CA1pyr and IC
%
% Currently chooses theta from REMstate as epochs tend to be shorter and theta is clean. 
%
% To-do: 
%   - create input to choose awake Theta state LB 10/2022
%
% inputs:
%   ica - output of ICA analysis
%   channel - channel of CA1 pyramidal layer
%   basepath - path to session that contains SleepStates.states.mat
%
%  (optional)
%   phaserange: vector of frequency range of interest
%               lowerbound:step:upperbound (default 5:0.5:12)
%   amprange:   vecto of amplitude range (default 30:5:200)
%
%   outputs:
%       saves figure of coupling across ICA components

p = inputParser;
addParameter(p,'phaserange',5:0.5:12,@isnumeric);
addParameter(p,'amprange', 30:5:200,@isnumeric);
addParameter(p,'state_tag','REMstate',@ischar);

parse(p,varargin{:});
phaserange = p.Results.phaserange;
amprange = p.Results.amprange;
state_tag = p.Results.state_tag; 

basename = basenameFromBasepath(basepath);

% load ica from basepath
% TO-DO

% [lfpTheta,lfpICA] = getLFP(pyrCh,'basepath',basepath,'interval',interval); % Only take 1000 seconds worth of data to keep the computation quick
% load lfp (all lfp)
lfp = load_lfp(basepath,ica.channels);

ica_array = analogSignalArray(...
    'data',[double(lfp.data(:,ica.channels == pyrCh))';  ica.unmixing * double(lfp.data')],...
    'timestamps',lfp.timestamps',...
    'sampling_rate',lfp.samplingRate);
    
% Choose intervals 
% behavior epochs 
epoch_df = load_epoch('basepath',basepath);
epoch = IntervalArray([epoch_df.startTime, epoch_df.stopTime]);

% find pre/post idx 
pre_task_post_mat = find_multitask_pre_post(epoch_df.environment,'pre_sleep_common',true); 

% load state of interest 
load(fullfile(basepath,[basename,'.SleepState.states.mat']))
state_epoch = IntervalArray(SleepState.ints.(state_tag));

pre_post_comod = {};
for task_n = 1:size(pre_task_post_mat,1)

    temp_idx = pre_task_post_mat(task_n,:); 
    
    % get all pre-task state ep 
    i = 1;
    pre_post_comod_temp = {};
    for pre_post  = [1,3]
        
        temp_ = ica_array(epoch(temp_idx(pre_post)) & state_epoch);

        lfp_.data = temp_.data;
        lfp_.timestamps = temp_.timestamps;
        lfp_.samplingRate = temp_.sampling_rate;
        lfp_.channels = 1:size(temp_.data,2);

        % Create the plot
        [comod_] = CFCPhaseAmp(lfp_,phaserange,amprange,...
            'phaseCh',1,...
            'ampCh',2:(size(ica.weights,1)+1),...
            'perm_test',false,...
            'units','MI',...
            'makePlot',true);
        pre_post_comod_temp{i} = comod_.comod;
        i = i+1;
    end
    
    pre_post_comod{task_n} = pre_post_comod_temp;
end

save(fullfile(save_path,[basename,'pre_post_ICA_CFC.mat']),'pre_post_comod')

end


function lfp = load_lfp(basepath,channelOrder)

basename = basenameFromBasepath(basepath);
session = loadSession(basepath,basename);

% find lfp file
lfp_path = dir(fullfile(basepath,'*.lfp'));
lfp_path = fullfile(lfp_path.folder,lfp_path.name);

% Load lfp
lfp.data = loadBinary(lfp_path,...
    'frequency',session.extracellular.srLfp,...
    'nchannels',session.extracellular.nChannels,...
    'start',0,...
    'duration',Inf,...
    'channels',channelOrder);

% lfp sample rate
Fs = session.extracellular.srLfp;

% prepare structure for below fuctions
lfp.samplingRate = Fs;
lfp.timestamps = (0:length(lfp.data)-1)/Fs;
lfp.channels = channelOrder; 

end