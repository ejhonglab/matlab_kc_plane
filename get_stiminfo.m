
function [updated] = get_stiminfo(thorimage_dir, thorsync_dir, ...
    output_dir, date, fly_num, update)

% Can be called with all arguments, fly_num missing, or date and fly_num
% missing.

% Gets only the final directory in the path, whether or not there is a trailing
% filesep character.
[~, thorimage_id, ~] = fileparts(strip(thorimage_dir, filesep));
thorimage_id = strip(thorimage_id, '_');

mat_output_dir = fullfile(output_dir, 'cnmf');
if ~exist(mat_output_dir, 'dir')
    mkdir(mat_output_dir);
end

mat_filename = sprintf('_%s_cnmf.mat', thorimage_id);
mat_filepath = fullfile(mat_output_dir, mat_filename);

if exist(mat_filepath, 'file')
    % TODO automate compiling of matwho. maybe fail rather than fallback to
    % default who.
    % Default 'who' is too slow with large .mat files.
    if (~update) & any(contains(matwho(mat_filepath), 'ti'))
        % TODO put behind verbose flag
        fprintf('\nti was already defined in %s\n', mat_filepath);
        updated = false;
        return;
    end

    create_mat = false;
else
    create_mat = true;
end

thorimage_xml_filepath = fullfile(thorimage_dir, 'Experiment.xml');

temp = xml2struct(thorimage_xml_filepath);
ThorImageExperiment = temp.ThorImageExperiment;

% TODO strip /, split on /, and take last part, in case full path passed in
tsync_filename = sprintf('%s.mat', thorsync_dir);

% Probably won't work w/ Windows specific filesep, and filesep probably wouldn't
% work with forward slashed paths in Windows.
[parent_dir, ~, ~] = fileparts(['/' strip(thorimage_dir, '/')]);
tsync_filepath = fullfile(parent_dir, 'thorsync', tsync_filename);
ai = load(tsync_filepath);
% TODO just load from hdf5?

% TODO had I commented this? should it be uncommented?
%ai.fpid = smoothdata(ai.pid, 'gaussian', 100);
%%

%%
ti.ts = str2double(ThorImageExperiment.Streaming.Attributes.frames);
ti.averageMode = str2double(ThorImageExperiment.LSM.Attributes.averageMode);
ti.averageNum = str2double(ThorImageExperiment.LSM.Attributes.averageNum);

fr = str2double(ThorImageExperiment.LSM.Attributes.frameRate);
if ti.averageMode == 1
   fr = fr / ti.averageNum; 
end
ti.fr = fr;

% TODO does pulsewidth always have the initcross either up or down?
% maybe don't use it if not? may also want fixed ~2.5v thresh...
[~,ic,fc,~] = pulsewidth(ai.scopePin);

% Filter out very short pulses, which are likely an artifact of some sort.
min_acquisition_length_s = 5.0;
% TODO get thorsync sampling interval (from where?)
% in ThorSync xml, need to find the (presumably just 1) active AcquireBoard,
% then find the one active SampleRate
% so far, this seems to always have been 30kHz. maybe check against successive
% ai.time readings (mode? first two?)
thorsync_sampling_rate = 30000;
min_acquisition_samples = round(thorsync_sampling_rate * ...
    min_acquisition_length_s);

ti.block_ic_idx = round(ic);
ti.block_fc_idx = round(fc);

% Vector of length equal to number of acquisitions, where each entry is the
% number of frames in that acquisition.
acquisition_frame_nums = ti.block_fc_idx - ti.block_ic_idx;

good_acquisitions = acquisition_frame_nums >= min_acquisition_samples;

ti.block_ic_idx = ti.block_ic_idx(good_acquisitions);
ti.block_fc_idx = ti.block_fc_idx(good_acquisitions);

ti.block_ict = ai.time(ti.block_ic_idx);
ti.block_fct = ai.time(ti.block_fc_idx);
% TODO just rename to num_blocks then... make central set of key term defs
% and always use those consistently
% # of blocks
ti.num_trials = length(ti.block_fct);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TODO TODO fix what caused ti.block_ict to be empty in 4-10/2/_001 case.
% something seriously wrong with that data?
if numel(ti.block_ict) > 0
    time = ai.time - ti.block_ict(1);
else
    time = ai.time;
end

interactive_plots = false;
if interactive_plots
    fsync = figure;
else
    fsync = figure('visible', 'off');
end
plot(time, ai.scopePin);
hold on;
plot(time, ai.olfDispPin);

% TODO handle this better... arguments are getting pretty fragile now and
% it's not like i'm really taking advantage of the fact they are optional
switch nargin
    case 6
        str = sprintf('%s, fly %d: %s', date, fly_num, thorimage_id);
    case 5
        % Assumed the last argument is the date str in this case.
        str = sprintf('%s: %s', date, thorimage_id);
    case 4
        str = sprintf('%s', thorimage_id);
    % TODO empty str basecase
end
title(str);
ylim([-1 6]);

% TODO include a legend
% TODO also plot mirror signal?

fig_output_dir = fullfile(output_dir, 'figures');
if ~exist(fig_output_dir, 'dir')
    mkdir(fig_output_dir);
end

fname = sprintf('%s_trial_structure.tif', thorimage_id);
print(fsync, fullfile(fig_output_dir, fname), '-dtiffn');

if numel(ti.block_ict) == 0
    error('No pulses detected in scopePin!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

last_block_last_frame = ti.block_fc_idx(end);

% TODO TODO ask Remy if she thinks subsetting olfDispPin like this might cause
% other issues
[~,ic,fc,~] = pulsewidth(ai.olfDispPin(1:last_block_last_frame));
ti.stim_ic_idx = round(ic);
ti.stim_fc_idx = round(fc);
ti.stim_ict = ai.time(ti.stim_ic_idx);
ti.stim_fct = ai.time(ti.stim_fc_idx);
ti.num_stim = length(ti.stim_fct);                  % # of odor pulses

ti.spt = ti.num_stim / ti.num_trials;               % stimuli per trial
ti.fpt = ti.ts / ti.num_trials;                       % frames per trial

%%

% TODO what exactly is Frame_Out? just one pulse per frame? other timing
% information meaningful (e.g. does high period include all scan time for
% frame?)
% TODO TODO ask Remy if she thinks subsetting olfDispPin like this might cause
% other issues
% TODO might need to check Frame_Out phase w.r.t. scopePin. frames might come
% out after scope pins goes low, and would want to account for that
% TODO TODO in that case, use some timepoint just before the scope pin goes high
% in a partial block, in the case where there is one partial block at the end
% TODO or add some fixed amount of timepoints beyond last_block_last_frame, but
% then might need to be careful if acquisition ends shortly after real last
% block
frameOutLogical = logical(ai.Frame_Out(1:last_block_last_frame));
% TODO maybe use pulsewidth to be consistent w/ rest of code? why not?
% or maybe eliminate pulsewidth?
frameOutDiff = diff(frameOutLogical);
% TODO more meaningful names
risingEdge = find(frameOutDiff > 0);
fallingEdge = find(frameOutDiff < 0);

% TODO more meaningful names
rLT = ai.time(risingEdge);
fLT = ai.time(fallingEdge);

% TODO more meaningful name
gap = diff(rLT);
% 0.0169 is the time between consective times in frame out counter.
% TODO TODO TODO may change w/ dimensions or some other (?) imaging settings.
% compute each time?
% Checking if gap is >10% more than expected.
% TODO TODO give idx a more meaningful name. what is it?
% TODO why use both this and scopePin?
idx = find(gap > (1.1 * .0169));
idx = [0 idx' numel(rLT)];

% TODO case where noise peaks are not the last peaks, but are still in idx, will
% probably not be handled correctly right now.

idx_by_stim = cell(ti.num_trials, 1);
for i = 1:ti.num_trials
    idx_by_stim{i} = (idx(i)+1):idx(i+1);
end

% TODO meaning of this variable? "av_"?
av_frames_by_stim = cell(ti.num_trials, 1);
times_by_stim = cell(ti.num_trials, 1);
% TODO TODO TODO frame_times used to be same length as movie, but after
% subsetting Frame_Out, it is not. fix frame_times / provide other information
% to match it up to whole movie.
frame_times = [];
if ti.averageMode == 1
   for i = 1:ti.num_trials
       % TODO rename idx here and below
       idx = idx_by_stim{i};
       % TODO is ti.averageNum a correct base case?
       % and isn't idx indexed as in ThorSync (with the much higher 30KHz), or
       % is it not?
       av_frames_by_stim{i} = idx(ti.averageNum):ti.averageNum:idx(end);
       times_by_stim{i} = fLT(av_frames_by_stim{i});
       frame_times = [frame_times; times_by_stim{i}];
   end
end
ti.frames_by_stim = idx_by_stim;
ti.av_frames_by_stim = av_frames_by_stim;
ti.times_by_stim = times_by_stim;
ti.frame_times = frame_times;

% TODO also rename these to blocks?
%% calculate trial_start and trial_end frames
ti.trial_end = cumsum(cellfun(@numel, ti.av_frames_by_stim'));
ti.trial_start = [0 ti.trial_end(1:end-1)] + 1;

%% calculate frames when olfDispPin is high
start_frame = zeros(1, ti.num_trials);
end_frame = zeros(1, ti.num_trials);
for i = 1:ti.num_stim
    idx = find(ti.frame_times >= ti.stim_ict(i) & ...
        ti.frame_times <= ti.stim_fct(i));

    start_frame(i) = idx(1);
    end_frame(i) = idx(end);
end
ti.stim_on = start_frame;
ti.stim_off = end_frame;

%% calculate baseline frames

% length of baseline (before stim on) in seconds
ti.baseline_length = 3;
ti.baseline_start_time = ti.stim_ict - ti.baseline_length;


start_frame = zeros(1, ti.num_trials);
end_frame = zeros(1, ti.num_trials);
for i = 1:ti.num_stim
    idx = find(ti.frame_times >= ti.baseline_start_time(i) & ...
        ti.frame_times < ti.stim_ict(i));

    start_frame(i) = idx(1);
    end_frame(i) = idx(end);
end
ti.baseline_start = start_frame;
ti.baseline_end = end_frame;
% TODO TODO print out what this is if going to keep doing this
% + better format? flag to suppress?
%disp([ti.baseline_start;ti.baseline_end]);

% TODO huh?
%% calculate peak response frames
ti.peakresp_length = 3;

start_frame = zeros(1,ti.num_trials);
end_frame = zeros(1,ti.num_trials);
for i = 1:ti.num_stim
    idx = find(ti.frame_times >= ti.stim_ict(i) & ...
        ti.frame_times <= (ti.stim_ict(i) + ti.peakresp_length));

    start_frame(i) = idx(1);
    end_frame(i) = idx(end);
end
ti.peakresp_start = start_frame;
ti.peakresp_end = end_frame;
% TODO TODO print out what this is if going to keep doing this
% + better format?
%disp([ti.peakresp_start;ti.peakresp_end]);
%% calculate response frames
ti.resp_length = 25;

start_frame = zeros(1,ti.num_trials);
end_frame = zeros(1,ti.num_trials);
for i = 1:ti.num_stim
    idx = find(ti.frame_times >= ti.stim_ict(i) & ...
        ti.frame_times <= (ti.stim_ict(i) + ti.resp_length));    

    start_frame(i) = idx(1);
    end_frame(i) = idx(end);
end
ti.resp_start = start_frame;
ti.resp_end = end_frame;
% TODO TODO print out what this is if going to keep doing this
% + better format?
%disp([ti.resp_start;ti.resp_end]);

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

% TODO maybe return from fn, and save outside? (to avoid needing files
% around...)

if create_mat
    save(mat_filepath, 'ti');
else
    save(mat_filepath, 'ti', '-append');
end
updated = true;

end
