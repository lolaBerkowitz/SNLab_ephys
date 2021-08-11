
function create_channelmap(basepath,varargin)
%  create a channel map file based on the .xml file. Original script from
%  kilosort 1 wrapper. Modified to not use channels that user set as 'skip'
%  in the .xml file. Modified further to accomodate probe geometry. 
%
%  Default map is for neuronexus A1x64-Poly2-6mm-23s-160 
%
%  Modified by Eliezyer de Oliveira, 02/03/2020
%  Modified by Laura Berkowitz, 08/11/2021

p = inputParser; 
p.addParameter('xcoords',[6,19,repmat([0,30],1,31)]');
p.addParameter('ycoords', [linspace(0,1449.0,64)]');
p.addParameter('Nshanks',1)
p.addParameter('fs',30000)
p.parse(varargin{:});
if ~exist('basepath','var')
    basepath = cd;
end

d   = dir('*.xml');
par = LoadXml(fullfile(basepath,d(1).name));

xcoords = p.Results.xcoords;
ycoords = p.Results.ycoords;
if ~isfield(par,'nElecGps')
    warning('No Electrode/Spike Groups found in xml.  Using Anatomy Groups instead.')
    tgroups = par.ElecGp;
    ngroups = length(tgroups);
else
    t = par.AnatGrps;
    ngroups = length(par.AnatGrps);
    for g = 1:ngroups
        tgroups{g} = par.AnatGrps(g).Channels;
    end
end

% Unpack to save
Nchannels = length(xcoords);
Nshanks = p.Results.Nshanks;
fs = p.Results.fs;

kcoords = (reshape(repmat(1:Nshanks, Nchannels/Nshanks, 1), Nchannels, 1));

connected   = true(Nchannels, 1); 
%getting channels labeled as a skip to set as disconnected
%with the code below I assume your anatomical group is the same as your
%spike group, so make sure it is if you want to skip channels correctly
aux_skip      = cell2mat({par.AnatGrps.Skip});
aux_ch        = cell2mat({par.AnatGrps.Channels});
disconnect_ch = aux_ch(logical(aux_skip))+1;

chanMap   = 1:Nchannels;
chanMap0ind = chanMap - 1;

connected(disconnect_ch) = false;

save(fullfile(basepath,'chanMap.mat'), ...
    'chanMap','connected', 'xcoords', 'ycoords', 'kcoords', 'chanMap0ind', 'fs')
end