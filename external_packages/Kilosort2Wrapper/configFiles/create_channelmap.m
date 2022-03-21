
function create_channelmap(basepath,varargin)
%  create a channel map file based on the .xml file. Original script from
%  kilosort 1 wrapper. Modified to not use channels that user set as 'skip'
%  in the .xml file. Modified further to accomodate probe geometry. 
%
%  Assumes xml file is named after raw data folder basename in line with 
%  buzcode data formatting standards https://github.com/buzsakilab/buzcode/wiki/Data-Formatting-Standards
%
%  Default map is for neuronexus A1x64-Poly2-6mm-23s-160. otherwise, add
%  xcoords/ycoords,Nshanks for given probe. 
%
%  Modified by Eliezyer de Oliveira, 02/03/2020
%  Modified by Laura Berkowitz, 08/11/2021

p = inputParser; 
p.addParameter('xcoords',[6,19,repmat([0,30],1,31)]');
p.addParameter('ycoords', [linspace(0,1449.0,64)]');
p.parse(varargin{:});

% unpack coordinates for channels 
xcoords = p.Results.xcoords;
ycoords = p.Results.ycoords;

% load basename.xml
[~,basename] = fileparts(basepath);
d   = dir([basepath,filesep,basename,'.xml']);
par = LoadXml(fullfile(basepath,d.name));

% Get nshanks from electrode groups or anatomical groups
if ~isfield(par,'nElecGps')
    Nshanks = length(par.ElecGp);
else
    warning('No Electrode/Spike Groups found in xml.  Using Anatomy Groups instead.')
    Nshanks = length(par.AnatGrps);
end

% Pull Nchannels and sample rate from xml
Nchannels = par.nChannels;
fs = par.SampleRate; 

kcoords = (reshape(repmat(1:Nshanks, Nchannels/Nshanks, 1), Nchannels, 1));

connected   = true(Nchannels, 1); 

%getting channels labeled as a skip to set as disconnected
%with the code below I assume your anatomical group is the same as your
%spike group, so make sure it is if you want to skip channels correctly
aux_skip      = cell2mat({par.AnatGrps.Skip});
aux_ch        = cell2mat({par.AnatGrps.Channels});
disconnect_ch = aux_ch(logical(aux_skip))+1;

chanMap   = aux_ch + 1;
chanMap0ind = chanMap - 1;

connected(disconnect_ch) = false;

save(fullfile(basepath,'chanMap.mat'), ...
    'chanMap','connected', 'xcoords', 'ycoords', 'kcoords', 'chanMap0ind', 'fs')
end