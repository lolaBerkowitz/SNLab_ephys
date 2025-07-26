
function plot_object_explore_results(basepath,explore_vectors,results)

% load position 
 basename = basenameFromBasepath(basepath);
 
 % load animal behavior and session files
 session = loadSession(basepath);
 load(fullfile(basepath,[basename,'.animal.behavior.mat']),'behavior')
 
 % load behavior epochs
 behave_ep = behavior_funcs.load_epochs(basepath);
 % load trials epochs
 trial_ep = behavior_funcs.load_trials(basepath);
 
 % get xy as analog signal array for easy epoching
 positions = analogSignalArray(...
     'data',[behavior.position.x;behavior.position.y],...
     'timestamps',behavior.timestamps,...
     'sampling_rate',behavior.sr);
    fig = figure; 
    fig.Color = [1 1 1];
 i = [1,2; 3,4];
   
 for ep = 1:behave_ep.n_intervals 
    % get intervalArray of current epoch
    cur_ep = behave_ep(ep) & trial_ep;   
    pos = positions(cur_ep);
    distance_vector = explore_vectors{1,ep}.object_distances;  
    cue_angle = explore_vectors{1,ep}.cue_angle;  
    
    angle_idx = cue_angle >= -45 & cue_angle <= 45;
    dist_idx = distance_vector <= 6;
 
    explore_idx = angle_idx & dist_idx;
    ax(i(ep,1)) = subplot(3,behave_ep.n_intervals,i(ep,1));
    plot(pos.data(:,1),pos.data(:,2),'Color','k')
    hold on;
    scatter(pos.data(explore_idx(1,:)',1),pos.data(explore_idx(1,:)',2),20,cue_angle(1,explore_idx(1,:))','filled')
    scatter(pos.data(explore_idx(2,:)',1),pos.data(explore_idx(2,:)',2),20,cue_angle(2,explore_idx(2,:))','filled')
    title('Object-exploration')
    c = colorbar;
    caxis([-180 180]);
    axis equal
    axis off
    colormap(ax(i(ep,1)), get_cmap/255)
    c.Label.String = 'Heading relative to cue';

    ax(i(ep,2)) = subplot(3,behave_ep.n_intervals,i(ep,2));
    [occ,~] = behavior_funcs.occ_map(pos.data(:,1),pos.data(:,2),2,60);
    b = imagesc(occ);
    set(b,'AlphaData',~isnan(occ))
    colormap(  ax(i(ep,2)),viridis(255))
    set(  ax(i(ep,2)),'YDir','normal')    
    title('Occupancy')
    axis equal
    axis off
    c = colorbar;
    c.Label.String = 'Occupancy (s)';
 
 end
 
    ax(5) = subplot(3,behave_ep.n_intervals,5);
    b = bar(explore_vectors{1, 1}.bin_explore','stacked');
    if ismember(results.moved_object_id(1),'A')
        legend({'Moved','Fixed'})
        b(1).FaceColor = [1 1 1];
        b(2).FaceColor = [.25 .25 .25];
    else
        legend({'Fixed','Moved'})
        b(2).FaceColor = [1 1 1];
        b(1).FaceColor = [.25 .25 .25];
    end
    xlabel('Time bin')
    ylabel('Exploration time (s)')
    title(['Baseline DR: ',num2str(results.DR_overall(1))])
    axis square
    
    ax(6) = subplot(3,behave_ep.n_intervals,6);
    b = bar(explore_vectors{1, 2}.bin_explore','stacked');
    if ismember(results.moved_object_id(1),'A')
        legend({'Moved','Fixed'})
        b(1).FaceColor = [1 1 1];
        b(2).FaceColor = [.25 .25 .25];
    else
        legend({'Fixed','Moved'})
        b(2).FaceColor = [1 1 1];
        b(1).FaceColor = [.25 .25 .25];
    end

    xlabel('Time bin (5 min/bin)')
    ylabel('Exploration time (s)')
    title(['Test DR: ',num2str(results.DR_overall(2))])
    axis square
end

function cmap = get_cmap()
cmap = [0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
1,0,1;
4,0,5;
7,1,9;
10,1,13;
13,2,17;
16,2,21;
19,2,25;
22,3,29;
26,3,33;
29,3,37;
32,4,41;
35,4,45;
38,5,49;
41,5,53;
44,5,57;
47,6,61;
50,6,65;
54,6,69;
57,7,73;
60,7,77;
63,8,81;
66,8,85;
69,8,89;
72,9,93;
75,9,97;
78,9,101;
81,10,105;
85,10,109;
88,11,113;
91,11,117;
94,11,121;
97,12,125;
100,12,129;
103,12,133;
106,13,137;
109,13,141;
112,14,145;
116,14,149;
119,14,153;
122,15,157;
125,15,161;
128,16,165;
131,16,170;
134,16,174;
137,17,178;
140,17,182;
144,17,186;
147,18,190;
150,18,194;
153,19,198;
156,19,202;
159,19,206;
162,20,210;
165,20,214;
168,20,218;
171,21,222;
175,21,226;
178,22,230;
181,22,234;
184,22,238;
187,23,242;
190,23,246;
193,23,250;
196,24,254;
196,24,254;
193,23,250;
190,23,246;
187,23,242;
184,22,238;
181,22,234;
178,22,230;
175,21,226;
171,21,222;
168,20,218;
165,20,214;
162,20,210;
159,19,206;
156,19,202;
153,19,198;
150,18,194;
147,18,190;
144,17,186;
140,17,182;
137,17,178;
134,16,174;
131,16,170;
128,16,166;
125,15,161;
122,15,157;
119,14,153;
116,14,149;
112,14,145;
109,13,141;
106,13,137;
103,12,133;
100,12,129;
97,12,125;
94,11,121;
91,11,117;
88,11,113;
85,10,109;
81,10,105;
78,9,101;
75,9,97;
72,9,93;
69,8,89;
66,8,85;
63,8,81;
60,7,77;
57,7,73;
54,6,69;
50,6,65;
47,6,61;
44,5,57;
41,5,53;
38,5,49;
35,4,45;
32,4,41;
29,3,37;
26,3,33;
22,3,29;
19,2,25;
16,2,21;
13,2,17;
10,1,13;
7,1,9;
4,0,5;
1,0,1;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0;
0,0,0];
end
%     
%     ax(5) = subplot(2,4,5);
%     scatter(pos.data(:,1),pos.data(:,2),10,distance_vector(2,:)','filled')
%     title('Position by distance')
%     colorbar
%     axis square 
%     axis off
%     colormap(ax(5),'jet')
