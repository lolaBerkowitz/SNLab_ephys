function complete_bool = check_lfp(basepath, varargin)
%check_lfp Verify completeness of LFP file relative to the DAT file
%   This function checks whether a downsampled LFP file (.lfp) is complete
%   by comparing it against the corresponding raw wideband data file (.dat).
%   If an LFP file exists, the function reads file sizes, computes the
%   expected downsampling ratio, and verifies that the ratio of samples in
%   the .dat and .lfp files matches the session metadata.
%
%   Syntax:
%   check_lfp(basepath)
%   check_lfp(basepath, 'ParameterName', ParameterValue, ...)
%
%   Inputs:
%   basepath - Path to directory containing the .dat and .lfp files
%              (default: current directory).
%
%   Optional Name-Value Parameters:
%   'datFile'  - Name of the .dat file (default: [basename '.dat'])
%   'inFs'     - Input sampling rate (default: session.extracellular.sr)
%   'outFs'    - LFP sampling rate (default: session.extracellular.srLfp)
%   'localDir' - Directory to search for the .lfp file 
%                (default: basepath)
%   'session'  - Session metadata structure 
%                (default: loaded via getSession.m from basepath)
%
%   Outputs:
%   boolean: If True, lfp file is complete.
%   - Confirms whether the .lfp file is complete.
%   - Warns if the LFP file exists but has fewer samples than expected.
%
%   Verification Procedure:
%   1. Load session metadata (sampling rates, number of channels).
%   2. Compute total samples in the .dat file.
%   3. If an .lfp file exists:
%        a. Compute total samples in the .lfp file.
%        b. Compare observed sample ratio with expected ratio 
%           (inFs/outFs).
%        c. Report whether the .lfp file is complete.
%   4. If no .lfp file exists, the function does nothing.
%
%   Example:
%   % Verify LFP file completeness in current directory
%   check_lfp
%
%   % Specify basepath
%   check_lfp('/data/mouse1/session1')
%
%   Dependencies:
%   - getSession.m (from CellExplorer)
%   - checkFile.m (utility to check file existence and size)
%
%   Note: This function does not generate a new .lfp file, but will
%   indicate if rewriting is necessary.

%% Input handling
if ~exist('basepath', 'var')
    basepath = pwd;
end
p = inputParser;
addParameter(p, 'datFile', [], @ischar);
addParameter(p, 'outFs', [], @isnumeric);
addParameter(p, 'inFs', [], @isnumeric);
addParameter(p, 'localDir', [], @isfolder);
addParameter(p, 'session', [], @isstruct);

parse(p, varargin{:})
datFile = p.Results.datFile;
outFs = p.Results.outFs;
inFs = p.Results.inFs;
localDir = p.Results.localDir;
session = p.Results.session;

if isempty(session)
    try
        session = getSession('basepath', basepath);
    catch ME
        error('LFPfromDat:SessionLoadError', 'Could not load session metadata. Please ensure getSession.m is available and a .session.mat file exists in the basepath. Error: %s', ME.message);
    end
end

[~ , basename] = fileparts(basepath);

if isempty(datFile)
    datFile = [basename, '.dat'];
elseif ~strcmp(datFile(end-3:end), '.dat')
    datFile = [datFile, '.dat'];
end


sizeInBytes = 2; % int16 is 2 bytes
%% Housekeeping
fInfo = checkFile('basepath', basepath, 'filename', datFile, 'searchSubdirs', false);

if isempty(inFs)
    inFs = session.extracellular.sr;
end

nbChan = session.extracellular.nChannels;

if isempty(outFs)
    outFs = session.extracellular.srLfp;
end

sampleRatio = (inFs / outFs);

% Output file
if ~isempty(localDir)
    flfp = fullfile(localDir, [basename, '.lfp']);
else
    flfp = fullfile(basepath, [basename, '.lfp']);
end


totalSamples = fInfo.bytes / (nbChan * sizeInBytes);

if exist(flfp, 'file') || ...
        exist([basepath, filesep, basename, '.eeg'], 'file')
    fprintf('LFP file already exists, checking if writing was complete \n')

    % Check if Lfp samples are equal to dat file 
    lfp_Info = checkFile('basepath', basepath, 'filename', [basename,'.lfp'], 'searchSubdirs', false);
    totalSamples_lfp = lfp_Info.bytes / (nbChan * sizeInBytes);
    
    % sampleRatio should be equal to totalSamples/totalSamples_lfp
    sampleRatio_obs = totalSamples/totalSamples_lfp;
    
    if ismembertol(sampleRatio,sampleRatio_obs,1e-6)
        fprintf('LFP file complete \n')
        complete_bool =  true;
    else
        fprintf('LFP file exists but not complete \n')
        complete_bool = false;
    end


end