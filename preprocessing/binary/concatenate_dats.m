function concatenate_dats(basepath,varargin)
%
%
% [concatenate_dats - Concatenate raw .dat files found in a session folder
% for intan type recordings] 
% 
% [ALGORITHM OUTLINE: looks for .dat files in a folder (or in subfolders) to
% concatenate together.  The concatenation happens via system commands 
% ("cat" command for linux/mac, "copy" command if windows/pc).  Uses
% different assumptions to find and recognize relevant .dats depending on
% the acquisition system]  
% 
% [REQUIREMENTS: Assumes you are in or pointed to a directory containing 
% subdirectories for various recording files from a single session. *It is 
% assumed that an earlier-acquired data file/folder will have a name that
% is sorted alphanumerically earlier.  Alphanumeric sorting order is
% assumed to be the recording temporal sequence.
% Works with acquisition systems: Intan  - 
%   1) intan: wherein subfolders are inside the session folder.  Each
%   subfolder contains simultaneously-recorded .dat files recorded for a
%   continuous period of time.  Start/stop recording commands each create a
%   new folder.  *It is assumed that the alphanumeric sorting of these 
%   folders corresponds with their sequence in acquisiton time.*  
%   These folders contain
%       - info.rhd files with metadata about the recording. 
%       - amplifier.dat - int16 file with usually neural data from the
%           headstage
%       - auxiliary.dat (optional) - uint16 file from auxiliary channels on
%           the headstage - often accelerometer
%       - analogin.dat (optional) - uint16 file from analogin channels on 
%           main board 
%       - digitalin.dat (optional) - uint16 file recording all 16 digital 
%           channels on main board 
%       - time.dat - int32 file giving recording sample indicies (e.g. 
%           0,1,2,3...) for each sample recorded in other channels
%       - supply.dat - uint16 file showing voltage supplied to (?preamp?)]
%   
%
%  USAGE
%
%    [concatenateDats(basepath,sortFiles)]
%
%  INPUTS
%
%    [basepath]        [computer path to session folder.  Defaults to
%                      current folder if no input given]
%    [sortFiles]               [boolean denoting whether to sort files according 
%                              to time of recording (1) or
%                              not (0) and thus sort them alphabetically 
%                              Default = 1]
%    [legacy]                  [there should only be 2 inputs, so this third
%                              input is here temporarily for legacy reasons.
%                              If for some reason you are using an old script
%                              which expects the "sortFiles" input to be 3rd,
%                              this will temporarily override your second input]
%                              
%                              
%
%  OUTPUT
%     [Operates on files in specified folder.  No output variable]
%
%  EXAMPLES
%      [Can be called directly or via bz_PreprocessExtracellEphysSession.m]
%
%  SEE ALSO
%
% [Brendon Watson, Antonio FR, Ralitsa Todorova, Laura Berkowitz] [2018-2024]
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.


%% Handling inputs
p = inputParser;
addParameter(p,'sortFiles',true,@islogical) 
addParameter(p, 'overwrite',false,@islogical)

parse(p,varargin{:});

sortFiles = p.Results.sortFiles;
overwrite = p.Results.overwrite;


[~, basename] = fileparts(basepath);


%% If the dats are already merged quit
if exist(fullfile(basepath,[basename,'.dat']),'file') & ~overwrite
    disp('.dat already exists in session directory, not merging subdats')
    return
end

%% Find all .dat paths in subfolders 

otherdattypes = {'analogin';'digitalin';'auxiliary';'time';'supply'};

bad_otherdattypes = [];

for odidx = 1:length(otherdattypes)
    newpaths.(otherdattypes{odidx}) = fullfile(basepath,[otherdattypes{odidx}, '.dat']);
end

% use dir to find folder contents
d = dir(basepath);
d(contains({d.name},{'.','..'})) = [];
datpaths = {};
datsizes.amplifier = [];
recordingnames = {};
rcount = 0; %Count of good subfolders
fun_existsdate = @(x) any(strfind(ismember(x,'1234567890'),[ones(1,6) 0 ones(1,6)])); % Intan's built in date format is 6 numbers followed by 6 numbers: YYMMDD_HHMMss

sr = []; % sampling rate 
for a = 1:length(d)
    
    %look in each subfolder
    if d(a).isdir && any(~ismember(d(a).name,'.')) && fun_existsdate(d(a).name) % d(a).name needs to be a directory without dots and contain a date to count as a session subfolder
       % grab amplifier sampling rate 
        [~, ~, ~, ~,...
            ~, ~, frequency_parameters,~ ] = ...
            read_Intan_RHD2000_file_snlab(fullfile(basepath,d(a).name));

        sr(a) = frequency_parameters.amplifier_sample_rate;
            
        %Check for amplifier.dat or subfolderbaseName.dat
        if exist(fullfile(basepath,d(a).name,[d(a).name,'.dat']),'file')
            ampfile = fullfile(basepath,d(a).name,[d(a).name,'.dat']);

        elseif exist(fullfile(basepath,d(a).name,'amplifier_analogin_auxiliary_int16.dat'),'file')%Luke Sjulson's Modified code to record all 16bit signals in one file
            ampfile = fullfile(basepath,d(a).name,'amplifier_analogin_auxiliary_int16.dat');

        else
            ampfile = fullfile(basepath,d(a).name,'amplifier.dat');
            digitalinfile = fullfile(basepath,d(a).name,'digitalin.dat');

        end

        if ~exist(ampfile,'file') && exist(digitalinfile,'file')
             % create an empty ampfile. This will allow us to preprocess sessions without recorded brain activity (but e.g. digital events)
            if ~exist(ampfile,'file'), fclose(fopen(ampfile, 'w')); end
        end
        if exist(ampfile,'file')
            rcount = rcount+1;
            datpaths.amplifier{rcount} = ampfile;
            t = dir(ampfile);
            datsizes.amplifier(rcount) = t.bytes;
            recordingnames{rcount} = d(a).name;
            for odidx = 1:length(otherdattypes)%loop through other .dat types found here
                % eval([otherdattypes{odidx} 'datpaths{rcount} = fullfile(basepath,recordingnames{rcount},''' otherdattypes{odidx} '.dat'');'])
                datpaths.(otherdattypes{odidx}){rcount} = fullfile(basepath,recordingnames{rcount},[otherdattypes{odidx} '.dat']);
                %eval(['d2 = dir(' otherdattypes{odidx} 'datpaths.amplifier{rcount});'])
                d2 = dir(datpaths.(otherdattypes{odidx}){rcount});
                if isempty(d2)
                    bad_otherdattypes(odidx) = 1;
                else
                    %eval([otherdattypes{odidx} 'datsizes(rcount) = d2(1).bytes;'])
                    datsizes.(otherdattypes{odidx})(rcount) = d2(1).bytes;
                end
            end
        end

    end
end

otherdattypes(find(bad_otherdattypes)) = [];%if there weren't analogin or digitalin in some recording
if isempty(datpaths) || isempty(datpaths.amplifier)
    disp('No .dats found in subfolders.  Exiting bz_ConcatenateDats.')
    return
end

%% Grab sampling rate 

if length(unique(sr)) > 1 
    error('Dat files are different sampling rates')
else
    sr = mode(sr);
end
%% Sort files according to time of recording

if sortFiles
    names2sort = cellfun(@(X) str2num(X(end-5:end)),recordingnames,'UniformOutput',false);
    names2sort = cell2mat(names2sort);
    if ~isempty(names2sort)
        disp('Assuming the last 6 digits reflect recording time.')
    else
        names2sort = 1:length(recordingnames);
        disp('Last 6 digits not numeric... sorting alphabetically')
    end
    [~,I] = sort(names2sort);
    recordingnames = recordingnames(I);
    datpaths.amplifier = datpaths.amplifier(I);
    datsizes.amplifier = datsizes.amplifier(I);
    for odidx = 1:length(otherdattypes)
        datpaths.(otherdattypes{odidx}) = datpaths.(otherdattypes{odidx})(I);
        datsizes.(otherdattypes{odidx}) = datsizes.(otherdattypes{odidx})(I);
    end
end

%% Concatenate

newdatpath = fullfile(basepath,[basename,'.dat']);
if isunix
    cs = strjoin(datpaths.amplifier);
    catstring = ['! cat ', cs, ' > ',newdatpath];
elseif ispc  
    if length(datpaths.amplifier)>1
        for didx = 1:length(datpaths.amplifier)-1
            datpathsplus{didx} = [datpaths.amplifier{didx} ' +'];
        end
        %Last file string shouldn't end with '+'
        datpathsplus{length(datpaths.amplifier)} = [datpaths.amplifier{length(datpaths.amplifier)}];
    else
        datpathsplus = datpaths.amplifier;
    end
    cs = strjoin(datpathsplus);
    catstring = ['! copy /b ', cs, ' ',newdatpath];
end

% action
disp('Concatenating Amplifier Dats... be patient')
eval(catstring)%execute concatention
    
%% Check that size of resultant .dat is equal to the sum of the components
t = dir(newdatpath);
if t.bytes ~= sum(datsizes.amplifier)
    error('New .dat size not right.  Exiting')
    return
else
    sizecheck.amplifier = true;
    disp('Primary .dats concatenated and size checked')
end

%% save times of each individual file concatenated

%% Also concatenate the other .dats
disp('Concatenating Other Dats..... continue to be patient')
for odidx = 1:length(otherdattypes)

    if isunix
        cs = strjoin(datpaths.(otherdattypes{odidx}));
        catstring = ['! cat ', cs, ' > ',newpaths.(otherdattypes{odidx})];
    elseif ispc%As of 4/9/2017 - never tested
        if length(datpaths.(otherdattypes{odidx}))>1
            for didx = 1:length(datpaths.(otherdattypes{odidx}))
                datpathsplus{didx} = [datpaths.(otherdattypes{odidx}){didx} '+'];
            end
          %Last file string shouldn't end with '+'
          datpathsplus{length(datpaths.(otherdattypes{odidx}))} = datpaths.(otherdattypes{odidx}){didx};
        else
            datpathsplus = datpaths.(otherdattypes{odidx});
        end
        cs = strjoin(datpathsplus);
        catstring = ['! copy /b ', cs, ' ',newpaths.(otherdattypes{odidx})];
    end

    
    eval(catstring)%execute concatenation

    % Check that size of resultant .dat is equal to the sum of the components
    t = dir(newpaths.(otherdattypes{odidx}));
    if t.bytes ~= sum(datsizes.(otherdattypes{odidx}))
        error(['New ' otherdattypes{odidx} '.dat size not right.  Exiting after .dats converted.  Not deleting'])
        sizecheck.(otherdattypes{odidx}) = false;
    else
        sizecheck.(otherdattypes{odidx}) = true;
        disp([otherdattypes{odidx} ' concatenated and size checked'])
    end
end


%% Get time points from the time.dat
%Use the timestamps from time.dat to get the sort order
%Number of samples in time.dat. First timepoint, last timepoint
%Convert from number of samples to recording time of start/ends
for ff = 1:length(datpaths.time)
    
	f = fopen(datpaths.time{ff},'r'); 
    % Determine total number of samples in file
    fileStart = ftell(f);
    
    %Read the first time point
    firsttimepoint = fread(f,1,'int32');
    status = fseek(f,-4,'eof'); %int32 = 4 bytes
    lasttimepoint = fread(f,1,'int32');
    fileStop = ftell(f);
    
    firstlasttimepoints(ff,:) = [firsttimepoint lasttimepoint];
    numsamples(ff) = fileStop./4;
    if ff==1
        transitiontimes_samp = firstlasttimepoints(ff,:);
    else
        transitiontimes_samp(ff,:) = firstlasttimepoints(ff,:)+transitiontimes_samp(ff-1,2)+1;
    end
end

disp(['Calculating merge times based on wideband samplingRate of ',num2str(sr),'Hz.'])
transitiontimes_sec = transitiontimes_samp./sr; %convert to seconds


%% Make the events.mat file that saves all the merge information

eventsfilename = fullfile(basepath,[basename,'.MergePoints.events.mat']);

MergePoints.timestamps = transitiontimes_sec;
MergePoints.timestamps_samples = transitiontimes_samp;
MergePoints.firstlasttimpoints_samples = firstlasttimepoints;
MergePoints.foldernames = recordingnames;
MergePoints.filesmerged = datpaths;
MergePoints.filesizes = datsizes;
MergePoints.sizecheck = sizecheck;
MergePoints.detectorinfo.detectorname = 'bz_ConcatenateDats';
MergePoints.detectorinfo.detectiondate = datestr(now,'yyyy-mm-dd');

%Saving SleepStates
save(eventsfilename,'MergePoints');


