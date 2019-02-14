
function get_stiminfo(thorimage_dir, thorsync_dir, output_dir,...
                      date, fly_num)
% Can be called with all arguments, fly_num missing, or date and fly_num
% missing.

% Gets only the final directory in the path, whether or not there is a trailing
% filesep character.
[~, thorimage_id, ~] = fileparts(strip(thorimage_dir, filesep));
thorimage_id = strip(thorimage_id, '_');

%thorimage_xml_filename = sprintf('_%s.xml', thorimage_id);
%thorimage_xml_filepath = fullfile(folder, 'xml_files', thorimage_xml_filename);

thorimage_xml_filepath = fullfile(thorimage_dir, 'Experiment.xml');

temp = xml2struct(thorimage_xml_filepath);
ThorImageExperiment = temp.ThorImageExperiment;
%%
% frame_num = str2double(ThorImageExperiment.Streaming.Attributes.frames);
% pixelX = str2double(ThorImageExperiment.LSM.Attributes.pixelX);
% pixelY = str2double(ThorImageExperiment.LSM.Attributes.pixelY);
% 
% siz = [pixelY pixelX frame_num];
% ti.ts = frame_num;
%% load thorsync

% TODO TODO strip /, split on /, and take last part, in case full path passed in
tsync_filename = sprintf('%s.mat', thorsync_dir);

% Probably won't work w/ Windows specific filesep, and filesep probably wouldn't
% work with forward slashed paths in Windows.
[parent_dir, ~, ~] = fileparts(['/' strip(thorimage_dir, '/')]);
tsync_filepath = fullfile(parent_dir, 'thorsync', tsync_filename);
ai = load(tsync_filepath);
% TODO just load from hdf5 as she does elsewhere

% TODO had I commented this? should it be uncommented?
%ai.fpid = smoothdata(ai.pid, 'gaussian', 100);
%%


%%
ti.ts = str2double(ThorImageExperiment.Streaming.Attributes.frames);
ti.averageMode = str2double(ThorImageExperiment.LSM.Attributes.averageMode);
ti.averageNum = str2double(ThorImageExperiment.LSM.Attributes.averageNum);

fr = str2double(ThorImageExperiment.LSM.Attributes.frameRate);
if ti.averageMode == 1
   fr = fr/ti.averageNum; 
end
ti.fr = fr;

[~,ic,fc,~] = pulsewidth(ai.scopePin);
ti.block_ic_idx = round(ic);
ti.block_fc_idx = round(fc);
ti.block_ict = ai.time(ti.block_ic_idx);
ti.block_fct = ai.time(ti.block_fc_idx);
ti.num_trials = length(ti.block_fct);                 % # of blocks

[~,ic,fc,~] = pulsewidth(ai.olfDispPin);
ti.stim_ic_idx = round(ic);
ti.stim_fc_idx = round(fc);
ti.stim_ict = ai.time(ti.stim_ic_idx);
ti.stim_fct = ai.time(ti.stim_fc_idx);
ti.num_stim = length(ti.stim_fct);                  % # of odor pulses

ti.spt = ti.num_stim / ti.num_trials;               % stimuli per trial
ti.fpt = ti.ts/ti.num_trials;                       % frames per trial


%%

% TODO what exactly is Frame_Out?
frameOutLogical = logical(ai.Frame_Out);
frameOutDiff = diff(frameOutLogical);
risingEdge = find(frameOutDiff>0);
fallingEdge = find(frameOutDiff<0);

rLT = ai.time(risingEdge);
fLT = ai.time(fallingEdge);

gap = diff(rLT);
% 0.0169 is the time between consective times in frame out counter.
% TODO may change w/ dimensions or some other (?) imaging settings. compute each
% time?
% Checking if gap is >10% more than expected.
idx = find(gap>(1.1*.0169));
idx = [0 idx' numel(rLT)];

idx_by_stim = cell(ti.num_trials,1);
for i = 1:ti.num_trials
    idx_by_stim{i} = (idx(i)+1):idx(i+1);
end

av_frames_by_stim = cell(ti.num_trials,1);
times_by_stim = cell(ti.num_trials,1);
frame_times = [];
if ti.averageMode == 1
   for i = 1:ti.num_trials
       idx = idx_by_stim{i};
       av_frames_by_stim{i} = idx(ti.averageNum):ti.averageNum:idx(end);
       times_by_stim{i} = fLT(av_frames_by_stim{i});
       frame_times = [frame_times; times_by_stim{i}];
   end
end
ti.frames_by_stim = idx_by_stim;
ti.av_frames_by_stim = av_frames_by_stim;
ti.times_by_stim = times_by_stim;
ti.frame_times = frame_times;

%% calculate trial_start and trial_end frames
ti.trial_end = cumsum(cellfun(@numel,ti.av_frames_by_stim'));
ti.trial_start = [0 ti.trial_end(1:end-1)]+1;

%% calculate frames when olfDispPin is high
start_frame = zeros(1,ti.num_trials);
end_frame = zeros(1,ti.num_trials);
for i = 1:ti.num_stim
    idx = find(ti.frame_times>=ti.stim_ict(i) & ti.frame_times<=ti.stim_fct(i));
    start_frame(i) = idx(1);
    end_frame(i) = idx(end);
end
ti.stim_on = start_frame;
ti.stim_off = end_frame;

%% calculate baseline frames

% length of baseline (before stim on) in seconds
ti.baseline_length = 3;
ti.baseline_start_time = ti.stim_ict-ti.baseline_length;


start_frame = zeros(1,ti.num_trials);
end_frame = zeros(1,ti.num_trials);
for i = 1:ti.num_stim
    idx = find(ti.frame_times>=ti.baseline_start_time(i) & ti.frame_times<ti.stim_ict(i));    
    start_frame(i) = idx(1);
    end_frame(i) = idx(end);
end
ti.baseline_start = start_frame;
ti.baseline_end = end_frame;
% TODO TODO print out what this is if going to keep doing this
% + better format? flag to suppress?
disp([ti.baseline_start;ti.baseline_end]);

%% calculate peak response frames
ti.peakresp_length = 3;

start_frame = zeros(1,ti.num_trials);
end_frame = zeros(1,ti.num_trials);
for i = 1:ti.num_stim
    idx = find(ti.frame_times>=ti.stim_ict(i) & ti.frame_times<=(ti.stim_ict(i)+ti.peakresp_length));    
    start_frame(i) = idx(1);
    end_frame(i) = idx(end);
end
ti.peakresp_start = start_frame;
ti.peakresp_end = end_frame;
% TODO TODO print out what this is if going to keep doing this
% + better format?
disp([ti.peakresp_start;ti.peakresp_end]);
%% calculate response frames
ti.resp_length = 25;

start_frame = zeros(1,ti.num_trials);
end_frame = zeros(1,ti.num_trials);
for i = 1:ti.num_stim
    idx = find(ti.frame_times>=ti.stim_ict(i) & ti.frame_times<=(ti.stim_ict(i)+ti.resp_length));    
    start_frame(i) = idx(1);
    end_frame(i) = idx(end);
end
ti.resp_start = start_frame;
ti.resp_end = end_frame;
% TODO TODO print out what this is if going to keep doing this
% + better format?
disp([ti.resp_start;ti.resp_end]);

%% stimulus info: arduino pins, stimulus index, pin_odors

% TODO TODO remove all stuff dealing w/ hard coded odor / pin info
%{
ti.pin_odors = cell(13,1);
ti.pin_odors{2} = 'ethyl acetate';
ti.pin_odors{3} = 'butanoic acid';
ti.pin_odors{4} = 'methyl acetate';
ti.pin_odors{5} = 'ethanol';
ti.pin_odors{6} = 'propyl acetate';
ti.pin_odors{8} = 'paraffin';
ti.pin_odors{9} = '3-methylbutan-1-ol';
ti.pin_odors{10} = '3-methylbutyl acetate';
ti.pin_odors{11} = ' propan-1-ol';

ti.channelA =  [2,2,9,2,9,2,9,2,2];
ti.channelB =  [8,9,8,8,8,9,8,9,8];
%%

odor_pairs = [ti.channelA' ti.channelB'];
odor_pairs = sort(odor_pairs,2);
[U, ia, ib]=unique(odor_pairs, 'stable', 'rows');

ti.unique_odor_pairs = U;    % unique odor pairs
ti.si = ib;
ti.pair_list = cell(size(U,1), 1);
for i = 1:size(U,1)
    ti.pair_list{i} =  [ti.pin_odors{U(i,1)} ' + ' ti.pin_odors{U(i,2)}];
end
%}

%%
% 
% ti.stim_list = cell(numel(ti.unique_odor_pairs), 2);
% for i = 1:numel(ti.unique_odor_pairs)
%     ti.stim_list{i,1} = ti.pin_odors(ti.unique_odor_pairs(i,1));
%     ti.stim_list{i,2} = ti.pin_odors(ti.unique_odor_pairs(i,2));
% end

% si = pins;  %stimulus index
% si(si==31)=13;
% si=si-1;
% 
% ti.channelA = pins;
% ti.si = si;
% 
% ti.pin_odors = odor_list;
%%

% TODO at least include a flag to disable plotting when used
% non-interactively...
time = ai.time - ti.block_ict(1);
fsync = figure;
plot(time, ai.scopePin);
% TODO remove?
hold on;
plot(time, ai.olfDispPin);

switch nargin
    case 5
        str = sprintf('%s, fly %d: %s', date, fly_num, thorimage_id);
    case 4
        % Assumed the last argument is the date str in this case.
        str = sprintf('%s: %s', date, thorimage_id);
    case 3
        str = sprintf('%s', thorimage_id);
end
title(str);
ylim([-1 6]);


fig_output_dir = fullfile(output_dir, 'figures');
if ~exist(fig_output_dir, 'dir')
    mkdir(fig_output_dir);
end

fname = sprintf('%s_trial_structure.tif', thorimage_id);
print(fsync, fullfile(fig_output_dir, fname), '-dtiffn');


mat_output_dir = fullfile(output_dir, 'cnmf');
if ~exist(mat_output_dir, 'dir')
    mkdir(mat_output_dir);
end

mat_filename = sprintf('_%s_cnmf.mat', thorimage_id);
mat_filepath = fullfile(mat_output_dir, mat_filename);

% TODO maybe return from fn, and save outside? (to avoid needing files
% around...)

if exist(mat_filepath, 'file')
    save(mat_filepath, 'ti', '-append');
else
    save(mat_filepath, 'ti');
end

end
