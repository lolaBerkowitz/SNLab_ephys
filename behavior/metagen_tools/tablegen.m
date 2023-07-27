
function output = tablegen (basepath)
    % Get all the videos from the directory
    file_struct = dir(fullfile(basepath,'*.avi'));
    %remove hidden files
    file_struct = file_struct(~startsWith({file_struct.name}, '.'));
    % Extract base name
    bstring=char(basepath);
    File_name=bstring(1+max(strfind(basepath,'/')):end); % ***** edit '\' or '/' depending on mac or windows
    % Extract the timestamps from the filenames
    timestamps = regexp({file_struct.name}, '-\d+', 'match');
    timestamps = cellfun(@(x) extractAfter(x,'-'),timestamps,'UniformOutput',false);
    timestamps = cellfun(@(x) str2double(x),timestamps,'UniformOutput',false);
    [~, sortedIndices] = sort( [timestamps{:}]);
    file_struct=file_struct(sortedIndices);
    subid=[];
    basepath3=[];
    genotype=[];
    age=[];
    session_date=[];
    dob=[];
    vidname=[];
    basename=[];
    exposure=[];
    maze_width_cm=[];
    maze_length_cm=[];
    maze_type=[];
    maze_color=[];
    cue_position=[];
    task_name=[];
    condition=[];
    treatment=[];
    objects=[];
    paradigm=[];
    moved_object=[];
    trial_start_1=[];
    trial_stop_1=[];
    trial_start_2=[];
    trial_stop_2=[];
    trial_start_3=[];
    trial_stop_3=[];
    trial_start_4=[];
    trial_stop_4=[];
    % none columns
    none=[];
    for i=1:length(sortedIndices)
        %nones and constant columns
        maze_length_cm=[maze_length_cm;39];
        none=[none;"none"];
        task_name=[task_name;"cpp"];
        %subid
        vid=file_struct(i).name;
        subid=[subid;vid(1:4)];
        %Basepath
        pos=find(sortedIndices==i);
        basepath3=[basepath3;string(basepath)];
        %Session Date
        daytemp=extractBefore(file_struct(pos).date," "); 
        date=datetime(daytemp, 'InputFormat', 'dd-MMM-yyyy');
        date=datestr(date, 'mm/dd/yy');
        session_date=[session_date;string(date)];
        %video name
        vidname=[vidname;string(vid)];
        %base name
        basename=[basename;string(File_name)];
        %exposure and condition
        
        if contains(basename,'pretest')
            
            if contains(vid,'pretest')
                condition=[condition;"pretest"];
                exposure=[exposure;1];
            elseif contains(vid,'pairing')|| contains(vid,'paring')
                exposure=[exposure;2];
                condition=[condition;"pairing"];
            end
         elseif contains(basename,'posttest')
            if contains(vid,'posttest')
                exposure=[exposure;5];
                condition=[condition;"posttest"];

            elseif contains(vid,'pairing')|| contains(vid,'paring')
                exposure=[exposure;4];
                condition=[condition;"pairing"];

            end
        elseif contains(basename,'day02')
            condition=[condition;"pairing"];
            exposure=[exposure;3];
        else
            condition=[condition;'pairinadsfafdsfsdg'];
            exposure=[exposure;9];
        end
        %maze width and type and color
        if contains(vid, 'pairing') || contains(vid, 'paring')
            maze_width_cm=[maze_width_cm;23];
            maze_type=[maze_type;"open_field"];
            if contains(vid, '_A-')
                maze_color=[maze_color;"yellow"];
            elseif contains(vid, '_B-')
                maze_color=[maze_color;"blue"];
            end
        else
            maze_width_cm=[maze_width_cm;70.5];
            maze_type=[maze_type;"dual_open_field"];
            maze_color=[maze_color;"multi"];
        end
    end
    %NaN and none columns
    na=cell(length(none), 1);
    trial_start_1=na;
    trial_stop_1=na;
    trial_start_2=na;
    trial_stop_2=na;
    trial_start_3=na;
    trial_stop_3=na;
    trial_start_4=na;
    trial_stop_4=na;
    genotype=na;
    age=na;
    dob=na;
    objects=none;
    cue_position=none;
    paradigm=none;
    moved_object=none;
    treatment=none;
    notes=na;
    output=table(subid,basepath3,genotype,age,session_date,dob,vidname,basename,exposure,...
        maze_width_cm,maze_length_cm,maze_type,maze_color,cue_position,task_name,condition,...
        treatment,objects,paradigm,moved_object,trial_start_1,trial_stop_1,trial_start_2,...
        trial_stop_2,trial_start_3,trial_stop_3,trial_start_4,trial_stop_4);
    