% Process DLC for postprocess

function [ts, x, y, angles] = process_DLC_for_ephys(basepath,varargin)
 % Unpacks DLC CSV, saves all position to structured matfile, and outputs
 % target coordinates for saving in CellExplorer format. 
% 
% Laura Berkowitz 2021

% Inputs:
%   basepath: path to raw data and dlc csv
%
% (Optional)
%   nvars: number of bodyparts tracked in DLC. Default is 2. 
%   target: cell array containing bodypart variable name to be saved for CellExplorer. 
%   compute_position: default is true. Will use two body parts to compute 
%        new position coordinate based on center of these two coordinates 
%        (ex/ right_ear & left_ear to compute head position). Must include
%        variable names in 'target'. 


p = inputParser;
p.addParameter('nvars',2,@isnumeric);
p.addParameter('target',{'cap_l','cap_r'},@ischar); %left and right ears
p.addParameter('compute_position',true,@islogical);

p.parse(varargin{:});

% unpack DLC csv and outs header and coordinates (note low liklihood
%% coordinates are converted to nan in this function
[header,tsxy]= open_dlc(basepath,nvars);

%% load digital events

% 
if compute_position 
    
    % pull cordinates for target
    p1 = tsxy(:,ismember(header(1,:),target{1}) & ismember(header(2,:),{'x','y'}));
    p2 = tsxy(:,ismember(header(1,:),target{2}) & ismember(header(2,:),{'x','y'}));

    % median between x and y is center of the head, so lets make those our new
    % x/y values. 
    x = nanmedian([p1(:,1) p2(:,1)],2)';
    y = nanmedian([p1(:,2) p2(:,2)],2)';
    
%      % Compute head angle
%     angles = rad2deg(XYangleLED(p1(:,1),p1(:,2),p2(:,1),p2(:,2)))';
%     if sum(ismember([front{1}],'nan')) == 3
%         angles = wrapTo360(angles-90); % account for markers on side of head instead of front/back
%     end

else
    % get x/y coordinates indicated by target
    p1 = tsxy(:,ismember(header(1,:),target{1}) & ismember(header(2,:),{'x','y'}));
end


% Sync start of spikes with video 
if ~isempty(dir([basepath,filesep,'digitalin.events.mat']))
   load([basepath,filesep,'digitalin.events.mat'])
    
% for neuralynx bitfields     
elseif ~isempty(dir([basepath,filesep,'*.ntt']))
    ts = Nlx2MatVT([basepath,filesep,'VT1.nvt'],[1,0,0,0,0,0],0,1);
    % remove trailing frames (Neuralynx capture has delay in closing video
    % capture, so frames may bigger). 
    ts = ts(1:length(angles));
else 
    % For now, just base timestamps off of DLC frames. These sessions are noted and will likely not be used given the ambiguous offset.  
    ts = [0:1/30:tsxy(end,1)/30]*10^6;
end

end

function [header,tsxy]= open_dlc(basepath,nvars)
% imports output from DLC tracking.  
%dlc_file is cell array of 'dir' output.
% for example,
% file = table2cell(struct2table(dir([path,filesep,'*.csv'])));
% where path is the path to the dlc output

file = table2cell(struct2table(dir([basepath,filesep,'*.csv'])));
dlc_idx = contains(file(:,1),'DLC') & ~contains(file(:,1),'above95th');

if sum(dlc_idx) == 0
    disp('No DLC tracking found')
    return
end

dlc_file = file(dlc_idx,:);

% load header
fileID = fopen(fullfile(dlc_file{2},dlc_file{1}),'r');
formatspec = ['%q',repmat('%q',1,nvars*3),'%[^\n\r]']; % nvars*3 as DLC gives x,y,liklihood
dataArray = textscan(fileID, formatspec,...
    3-2+1, 'Delimiter',',', 'TextType', 'string', 'HeaderLines',...
    2-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
header = [dataArray{1:end-1}];
clearvars fileID formatspec dataArray ans;

% load data
fileID = fopen(fullfile(dlc_file{2},dlc_file{1}),'r');
formatspec = ['%f',repmat('%f',1,nvars*3),'%[^\n\r]'];
dataArray = textscan(fileID, formatspec,...
    'Delimiter', ',','TextType', 'string', 'EmptyValue', NaN,...
    'HeaderLines' ,4-1,'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
tsxy = [dataArray{1:end-1}];

% Turn low liklihood points into NAN and smooth remaining coords
idx = tsxy(:,contains(header(2,:),'likelihood')) < .95;
xloc = find(contains(header(2,:),'x'));
yloc = find(contains(header(2,:),'y'));
for l = 1:size(idx,2)
    tsxy(idx(:,l),xloc(l):yloc(l))=NaN;
end

end