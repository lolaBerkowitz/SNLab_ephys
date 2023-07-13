basepath = 'Y:\laura_berkowitz\behavior_validation\appps1_cpp\data\cohort 2\3775N_pretest-07122023121951DLC_resnet50_mouse_cppJun28shuffle1_175000_filtered.csv';
save_path = 'Y:\laura_berkowitz\behavior_validation\appps1_cpp\data\cohort 2';

df = load_dlc_csv(basepath);
fs = 30;

% Choose neck 

df_neck = df(:,contains(fieldnames(df),'neck'));
df_neck.neck_x_3(df_neck.neck_2_likelihood_3 < .95) = nan;
df_neck.neck_1_y_3(df_neck.neck_2_likelihood_3 < .95) = nan;

ts = linspace(0,length(df_neck.neck_x_3)/fs,length(df_neck.neck_x_3));
x = df_neck.neck_x_3;
y = df_neck.neck_1_y_3;

good_idx = manual_trackerjumps(ts,...
    x,...
    y,...
    0,... % start
   length(df_neck.neck_x_3),... % end
    save_path,'darkmode',false);

x(~good_idx) = nan;
y(~good_idx) = nan;


arena_A = x < 435;
arena_B = x > 930;
% convert xy to cm from pixel space 
% get maxmin of x



ts_dff = [diff(ts) min(diff(ts))];

ts_A = sum(ts_dff(arena_A));
ts_B = sum(ts_dff(arena_B));

DR = (ts_A - ts_B) / (ts_A + ts_B);



disp(DR)