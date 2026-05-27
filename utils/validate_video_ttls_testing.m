function report = validate_video_ttls_testing(basepath, varargin)

% VALIDATE_VIDEO_TTLS
%
% Checks:
%   1) TTLs form contiguous video segments
%   2) Epoch channels contain valid start/end pairs
%   3) Each TTL segment falls within ONE epoch
%   4) TTL segments match AVI frame counts
%   5) Detects orphan TTLs (camera on but not recording)
%   6) Detects duplicate/spurious epoch TTLs across channels
%
% OUTPUT:
%   report struct containing all detected issues

%% Parameters
p = inputParser;

addParameter(p, 'minFrames', 50);
addParameter(p, 'gapThresh', 50);
addParameter(p, 'frameTolerance', 2);

parse(p, varargin{:});

minFrames = p.Results.minFrames;
gapThresh = p.Results.gapThresh;
frameTol = p.Results.frameTolerance;

%% Load data

load(fullfile(basepath,'digitalIn.events.mat'), ...
    'digitalIn');

video_ts = digitalIn.timestampsOn{1,1};

% epoch channels
numEpochSets = size(digitalIn.timestampsOn,2) - 1;

epochSets = cell(1, numEpochSets);

for i = 1:numEpochSets
    epochSets{i} = digitalIn.timestampsOn{1,i+1};
end

%% ------------------------------------------------------------------------
% 1) Segment video TTLs
% -------------------------------------------------------------------------

d = diff(video_ts);

if isempty(gapThresh)
    gapThresh = median(d) * 5;
end

breakIdx = find(d > gapThresh);

segmentStarts = [1; breakIdx+1];
segmentEnds   = [breakIdx; length(video_ts)];

segments = struct();

for i = 1:length(segmentStarts)

    idx = segmentStarts(i):segmentEnds(i);

    segments(i).idx = idx;
    segments(i).start_ts = video_ts(idx(1));
    segments(i).end_ts = video_ts(idx(end));
    segments(i).nFramesTTL = length(idx);

end

%% ------------------------------------------------------------------------
% 2) Load AVI files + frame counts
% -------------------------------------------------------------------------

videoFiles = dir(fullfile(basepath, '*.avi'));

videoInfo = struct();

for v = 1:length(videoFiles)

    fname = fullfile(basepath, videoFiles(v).name);

    try

        vr = VideoReader(fname);

        % fast estimate
        nFrames = floor(vr.Duration * vr.FrameRate);

        videoInfo(v).name = videoFiles(v).name;
        videoInfo(v).nFrames = nFrames;

    catch

        warning('Could not read video: %s', fname);

        videoInfo(v).name = videoFiles(v).name;
        videoInfo(v).nFrames = NaN;

    end
end

% sort alphabetically
[~, idx] = sort({videoInfo.name});
videoInfo = videoInfo(idx);

%% ------------------------------------------------------------------------
% 3) Match TTL segments to videos
% -------------------------------------------------------------------------
videoMatchIssues = [];

segment_list = [];
for s = 1:length(segments)
    segment_list = [segment_list segments(s).nFramesTTL];
end

for n_vid = 1:length(videoInfo)

    if any(abs(videoInfo(n_vid).nFrames - segment_list) <= 40)
        idx = find(abs(videoInfo(n_vid).nFrames - segment_list));
        videoMatchIssues(n_vid).frame_diff = abs(videoInfo(n_vid).nFrames - segment_list(idx));

    end
end


%% ------------------------------------------------------------------------
% 4) Check epoch integrity
% -------------------------------------------------------------------------

epochIssues = [];

for m = 1:numEpochSets

    ep = epochSets{m};

    % odd number of timestamps
    if mod(length(ep),2) ~= 0

        epochIssues(end+1).mouse = m; %#ok<AGROW>
        epochIssues(end).issue = ...
            'Odd number of epoch timestamps';

    end
end

%% ------------------------------------------------------------------------
% 5) Detect duplicate/spurious epoch TTLs across channels
% -------------------------------------------------------------------------

spuriousEpochTTLs = [];

for m1 = 1:numEpochSets-1

    ep1 = epochSets{m1};

    for m2 = (m1+1):numEpochSets

        ep2 = epochSets{m2};

        [sharedTS, idx1, idx2] = intersect(ep1, ep2);

        if ~isempty(sharedTS)

            for k = 1:length(sharedTS)

                spuriousEpochTTLs(end+1).timestamp = ...
                    sharedTS(k); %#ok<AGROW>

                spuriousEpochTTLs(end).mouse1 = m1;
                spuriousEpochTTLs(end).mouse2 = m2;

                spuriousEpochTTLs(end).idx1 = idx1(k);
                spuriousEpochTTLs(end).idx2 = idx2(k);

                % determine start/end labels
                if mod(idx1(k),2)==1
                    spuriousEpochTTLs(end).type1 = 'start';
                else
                    spuriousEpochTTLs(end).type1 = 'end';
                end

                if mod(idx2(k),2)==1
                    spuriousEpochTTLs(end).type2 = 'start';
                else
                    spuriousEpochTTLs(end).type2 = 'end';
                end
            end
        end
    end
end

%% ------------------------------------------------------------------------
% 6) Assign segments to epochs
% -------------------------------------------------------------------------

assignmentIssues = [];

for s = 1:length(segments)

    seg = segments(s);

    matched = false;
    matches = [];

    for m = 1:numEpochSets

        ep = epochSets{m};

        % skip malformed epochs
        if mod(length(ep),2) ~= 0
            continue
        end

        for e = 1:2:length(ep)

            start_ep = ep(e);
            end_ep = ep(e+1);

            if seg.start_ts >= start_ep && ...
                    seg.end_ts <= end_ep

                matched = true;

                matches = [matches; m e]; %#ok<AGROW>

            end
        end
    end

    if ~matched

        assignmentIssues(end+1).segment = s; %#ok<AGROW>
        assignmentIssues(end).issue = ...
            'No matching epoch';

    elseif size(matches,1) > 1

        assignmentIssues(end+1).segment = s; %#ok<AGROW>
        assignmentIssues(end).issue = ...
            'Segment matches multiple epochs';

    else

        segments(s).assignedMouse = matches(1,1);
        segments(s).epochIdx = matches(1,2);

    end
end

%% ------------------------------------------------------------------------
% 7) Detect orphan TTL segments
% -------------------------------------------------------------------------

orphanSegments = [];

for s = 1:length(segments)

    if segments(s).nFramesTTL < minFrames

        orphanSegments(end+1) = s; %#ok<AGROW>

    end
end

%% ------------------------------------------------------------------------
% Build report
% -------------------------------------------------------------------------

report = struct();

report.segments = segments;

report.videoInfo = videoInfo;

% report.videoMatchIssues = videoMatchIssues;

report.epochIssues = epochIssues;

report.spuriousEpochTTLs = spuriousEpochTTLs;

report.assignmentIssues = assignmentIssues;

report.orphanSegments = orphanSegments;

%% ------------------------------------------------------------------------
% Print summary
% -------------------------------------------------------------------------

fprintf('\n========================================\n');
fprintf('TTL + VIDEO VALIDATION REPORT\n');
fprintf('========================================\n\n');

% fprintf('Detected TTL segments: %d\n', nSeg);
% fprintf('Detected AVI files:    %d\n\n', nVid);

fprintf('Video matching issues:        %d\n', ...
    length(videoMatchIssues));

fprintf('Epoch integrity issues:       %d\n', ...
    length(epochIssues));

fprintf('Spurious epoch TTL overlaps:  %d\n', ...
    length(spuriousEpochTTLs));

fprintf('Epoch assignment issues:      %d\n', ...
    length(assignmentIssues));

fprintf('Orphan TTL segments:          %d\n\n', ...
    length(orphanSegments));

%% Print detailed issues

if ~isempty(videoMatchIssues)

    fprintf('\n--- VIDEO MATCH ISSUES ---\n');

    disp(videoMatchIssues)

end

if ~isempty(epochIssues)

    fprintf('\n--- EPOCH ISSUES ---\n');

    disp(epochIssues)

end

if ~isempty(spuriousEpochTTLs)

    fprintf('\n--- SPURIOUS EPOCH TTLS ---\n');

    for i = 1:length(spuriousEpochTTLs)

        s = spuriousEpochTTLs(i);

        fprintf(['Timestamp %.6f shared between ' ...
            'Mouse %d (%s idx=%d) and ' ...
            'Mouse %d (%s idx=%d)\n'], ...
            s.timestamp, ...
            s.mouse1, s.type1, s.idx1, ...
            s.mouse2, s.type2, s.idx2);

    end
end

if ~isempty(assignmentIssues)

    fprintf('\n--- ASSIGNMENT ISSUES ---\n');

    disp(assignmentIssues)

end

if ~isempty(orphanSegments)

    fprintf('\n--- ORPHAN TTL SEGMENTS ---\n');

    disp(orphanSegments)

end

fprintf('\nValidation complete.\n\n');

end