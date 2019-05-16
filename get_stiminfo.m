
function [updated] = get_stiminfo(thorimage_dir, thorsync_dir, ...
    output_dir, date, fly_num, update)

% TODO TODO this is definitely failing somehow in 2019-01-18/2/* cases
% figure out why. may impact future / other current experiments.

% Can be called with all arguments, fly_num missing, or date and fly_num
% missing.

% Gets only the final directory in the path, whether or not there is a trailing
% filesep character.
[~, thorimage_id, ~] = fileparts(strip(thorimage_dir, filesep));
%thorimage_id = strip(thorimage_id, '_');

mat_output_dir = fullfile(output_dir, 'cnmf');
if ~exist(mat_output_dir, 'dir')
    mkdir(mat_output_dir);
end

mat_filename = sprintf('%s_cnmf.mat', thorimage_id);
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

% <v4 XML format
if isstruct(ThorImageExperiment.Streaming)
    ti.ts = str2double(ThorImageExperiment.Streaming.Attributes.frames);

% v4+ XML format
elseif iscell(ThorImageExperiment.Streaming)
    ti.ts = str2double(ThorImageExperiment.Streaming{1}.Attributes.frames);

% TODO else error
end

% <v4 XML format
if isstruct(ThorImageExperiment.LSM)
    % TODO get LSM object -> share more code
    ti.averageMode = str2double(ThorImageExperiment.LSM.Attributes.averageMode);
    ti.averageNum = str2double(ThorImageExperiment.LSM.Attributes.averageNum);
    fr = str2double(ThorImageExperiment.LSM.Attributes.frameRate);

% v4+ XML format
elseif iscell(ThorImageExperiment.LSM)
    ti.averageMode = str2double(...
        ThorImageExperiment.LSM{1}.Attributes.averageMode);

    ti.averageNum = str2double(...
        ThorImageExperiment.LSM{1}.Attributes.averageNum);

    fr = str2double(ThorImageExperiment.LSM{1}.Attributes.frameRate);

% TODO else fail
end

if ti.averageMode == 1
   fr = fr / ti.averageNum; 
end
ti.fr = fr;

% TODO does pulsewidth always have the initcross either up or down?
% maybe don't use it if not? may also want fixed ~2.5v thresh...
[~, ic, fc, ~] = pulsewidth(ai.scopePin);

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

% TODO off by one? add one?
last_block_end_index = ti.block_fc_idx(end);
last_block_end_time = ti.block_fct(end);

max_trigger_to_last_frame_sec = 1;
last_fall_check_time = last_block_end_time + max_trigger_to_last_frame_sec;

if any(strcmp(fieldnames(ai), 'Frame_Out'))
    frameOutLogical = logical(ai.Frame_Out(ai.time <= last_fall_check_time));
elseif any(strcmp(fieldnames(ai), 'FrameOut'))
    frameOutLogical = logical(ai.FrameOut(ai.time <= last_fall_check_time));
% TODO else error
end
% TODO maybe use pulsewidth to be consistent w/ rest of code? why not?
% or maybe eliminate pulsewidth?
frameOutDiff = diff(frameOutLogical);
rising_edge_indices_frame_out = find(frameOutDiff > 0);
falling_edge_indices_frame_out = find(frameOutDiff < 0);

rising_edge_time_frame_out = ai.time(rising_edge_indices_frame_out);
falling_edge_time_frame_out = ai.time(falling_edge_indices_frame_out);

% TODO maybe just filter out beginning artifacts if this case is encountered?
% frame_out should start low (rise before falling)
assert(rising_edge_time_frame_out(1) < falling_edge_time_frame_out(1));
assert(numel(rising_edge_time_frame_out) == numel(falling_edge_time_frame_out));
%frame_out_length = falling_edge_time_frame_out - rising_edge_time_frame_out;

low2high_times = (rising_edge_time_frame_out(2:end) - ...
                  falling_edge_time_frame_out(1:(end - 1)));

% TODO maybe hardcode / approximate this for speed
max_frame_out_low2high_interval = mode(low2high_times) * 1.5;

falls_to_check = falling_edge_time_frame_out >= last_block_end_time & ...
                 falling_edge_time_frame_out <= last_fall_check_time;
fall_times_to_check = falling_edge_time_frame_out(falls_to_check);
fall_indices_to_check = falling_edge_indices_frame_out(falls_to_check);

% TODO delete
%{
disp('last_block_end_time')
disp(last_block_end_time)
disp('last_fall_check_time')
disp(last_fall_check_time)
disp('max_frame_out_low2high_interval');
disp(max_frame_out_low2high_interval);
disp('fall_times_to_check');
disp(fall_times_to_check);
%}

% TODO vectorized way to do this?
last_block_fall_time = nan;
for i = 1:length(fall_times_to_check)
    fall_time = fall_times_to_check(i);
    max_time = fall_time + max_frame_out_low2high_interval;

    rising_edge_logical = rising_edge_time_frame_out >= fall_time ...
        & rising_edge_time_frame_out <= max_time;

    %{
    disp('fall time:');
    disp(fall_time);
    disp('max_time')
    disp(max_time)
    disp('rise times in window:');
    disp(rising_edge_time_frame_out(rising_edge_logical));
    %}

    if sum(rising_edge_logical) == 0
        last_block_fall_time = fall_time;
        last_block_fall_index = fall_indices_to_check(i);
        break
    end
end

% TODO TODO why is this getting triggered in the v4 case?
% (just the stuff that had only 3000 frames because max frames issue?)
% TODO TODO TODO and why, similarly often, are no scopePin pulses detected?
if isnan(last_block_fall_time)
    error(['Frame Out continued pulsing right until end of recording. ' ...
          'The recording may have stopped incorrectly.']);
end

% TODO delete
%{
disp('last_block_fall_time')
disp(last_block_fall_time)
disp('last_block_fall_index')
disp(last_block_fall_index)

figure;
buff = 4000;
start = last_block_end_index - buff;
stop = last_block_end_index + buff;
plot_times = ai.time(start:stop);
plot(plot_times, frameOutLogical(start:stop))
hold on;
plot(plot_times, ai.scopePin(start:stop))
uiwait;
%}

% TODO right now, this does not include the double pulse i've seen at the
% last pulse of frame_out in at least one recording. change if necessary.
last_block_end_index = last_block_fall_index;

% TODO more meaningful name
gap = diff(rising_edge_time_frame_out);

% TODO delete. for debugging.
%{
disp('max(gap)')
disp(max(gap))
disp('min(gap)')
disp(min(gap))
disp('mode(gap)')
disp(mode(gap))
% TODO TODO check against something else that this is actually the # of frames
disp('numel(rising_edge_time_frame_out)')
disp(numel(rising_edge_time_frame_out))
disp('numel(falling_edge_time_frame_out)')
disp(numel(falling_edge_time_frame_out))
%}
%

% 0.0169 is the time between consective times in frame out counter.
% TODO TODO TODO may change w/ dimensions or some other (?) imaging settings.
% compute each time?
% Checking if gap is >10% more than expected.
% TODO TODO give idx a more meaningful name. what is it?
% TODO why use both this and scopePin?

% TODO TODO rename to indicate which side of block_boundary this is on
% (index of start/end of last frame/first frame of next block)
idx = find(gap > (1.1 * .0169));
idx = [0 idx' numel(rising_edge_time_frame_out)];

% TODO delete
%{
disp('idx(idx >= start & idx <= stop)')
disp(idx(idx >= start & idx <= stop))
%}

% TODO is there always a double peak at end of acquistion?
% (high period ~2x as long) is that intentional? should it be counted?
% counted twice?
%

% TODO case where noise peaks are not the last peaks, but are still in idx, will
% probably not be handled correctly right now.

idx_by_stim = cell(ti.num_trials, 1);
for i = 1:ti.num_trials
    idx_by_stim{i} = (idx(i) + 1):idx(i + 1);
end
% TODO delete. for debugging.
%{
disp('idx_by_stim')
disp(idx_by_stim)
%}
%

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
       % TODO TODO TODO are frameout pulses always in number that is a multiple
       % of averageNum??? need to handle case where there might be spurious
       % extra frames?
       % TODO is ti.averageNum a correct base case?
       % and isn't idx indexed as in ThorSync (with the much higher 30KHz), or
       % is it not?
       av_frames_by_stim{i} = idx(ti.averageNum):ti.averageNum:idx(end);

       % TODO delete. for debugging.
       %{
       disp('size(av_frames_by_stim{i})')
       disp(size(av_frames_by_stim{i}))
       disp('size(av_frames_by_stim)')
       disp(size(av_frames_by_stim))
       disp('size(falling_edge_time_frame_out)')
       disp(size(falling_edge_time_frame_out))
       disp('size(times_by_stim)')
       disp(size(times_by_stim))
       disp('i')
       disp(i)
       disp('max(av_frames_by_stim{i})')
       disp(max(av_frames_by_stim{i}))
       %}
       %
       times_by_stim{i} = falling_edge_time_frame_out(av_frames_by_stim{i});
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

%{
disp('size(ai.olfDispPin)')
disp(size(ai.olfDispPin))
disp('size(ai.olfDispPin(1:last_block_end_index))')
disp(size(ai.olfDispPin(1:last_block_end_index)))
%}

%% calculate frames when olfDispPin is high
[~, ic, fc, ~] = pulsewidth(ai.olfDispPin(1:last_block_end_index));
ti.stim_ic_idx = round(ic);
ti.stim_fc_idx = round(fc);
ti.stim_ict = ai.time(ti.stim_ic_idx);
ti.stim_fct = ai.time(ti.stim_fc_idx);
% # of odor pulses
ti.num_stim = length(ti.stim_fct);

start_frame = zeros(1, ti.num_trials);
end_frame = zeros(1, ti.num_trials);

%{
disp('ti.num_stim')
disp(ti.num_stim)
%}
for i = 1:ti.num_stim
    idx = find(ti.frame_times >= ti.stim_ict(i) & ...
        ti.frame_times <= ti.stim_fct(i));

    start_frame(i) = idx(1);
    end_frame(i) = idx(end);
end
ti.stim_on = start_frame;
ti.stim_off = end_frame;

if create_mat
    save(mat_filepath, 'ti');
else
    save(mat_filepath, 'ti', '-append');
end
updated = true;

end
