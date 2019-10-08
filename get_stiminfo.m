
function [updated] = get_stiminfo(thorimage_dir, thorsync_dir, ...
    output_dir, date, fly_num, update)
% Can be called with all arguments, fly_num missing, or date and fly_num
% missing.

% Output values used by populate_db.py/gui.py (as of 2019-08-02):
% - frame_times
% - trial_start
% - trial_end
% - block_start_sample (as block_ic_idx)
% - stim_on
% - stim_off

write_anything = true;
debug = false;

%{
thorimage_dir = '/mnt/nas/mb_team/raw_data/2019-08-27/9/fn_0002';
thorsync_dir = 'SyncData002';
output_dir = '/mnt/nas/mb_team/analysis_output/2019-08-27/9';
update = false;
write_anything = true;
%}

% TODO TODO this is definitely failing somehow in 2019-01-18/2/* cases
% figure out why. may impact future / other current experiments.

if debug
    disp('');
    disp(sprintf("thorimage_dir='%s';", thorimage_dir));
    disp(sprintf("thorsync_dir='%s';", thorsync_dir));
    disp(sprintf("output_dir='%s';", output_dir));
    disp(sprintf("update=%s;", string(update)));
end

% Gets only the final directory in the path, whether or not there is a trailing
% filesep character.
[~, thorimage_id, ~] = fileparts(strip(thorimage_dir, filesep));

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

%ai.fpid = smoothdata(ai.pid, 'gaussian', 100);

% <v4 XML format
if isstruct(ThorImageExperiment.Streaming)
    streaming = ThorImageExperiment.Streaming;
% v4+ XML format
elseif iscell(ThorImageExperiment.Streaming)
    streaming = ThorImageExperiment.Streaming{1};
% TODO else error
end
ti.ts = str2double(streaming.Attributes.frames);
num_frames = str2double(streaming.Attributes.frames);
% TODO need to subtract flyback frames to get effective #?
% or done automatically?

% <v4 XML format
if isstruct(ThorImageExperiment.LSM)
    lsm = ThorImageExperiment.LSM;
% v4+ XML format
elseif iscell(ThorImageExperiment.LSM)
    lsm = ThorImageExperiment.LSM{1};
% TODO else fail
end
ti.averageMode = str2double(lsm.Attributes.averageMode);
ti.averageNum = str2double(lsm.Attributes.averageNum);
fr = str2double(lsm.Attributes.frameRate);

if ti.averageMode == 0
    ti.averageNum = 1;
elseif ti.averageMode == 1
   fr = fr / ti.averageNum;
   % 2019-09-09: Was something like this division happening before?
   % It seems num_frames is already the count post-averaging.
   % Different source for this variable before or something?
   %assert(~mod(num_frames, ti.averageNum));
   %num_frames = num_frames / ti.averageNum;
end
ti.fr = fr;
ti.num_frames = num_frames;

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

% TODO TODO TODO rename all _idx stuff to either _frame or _sample, depending on
% whether it's thorimage or thorsync, respectively
% TODO (would need to accept both names in ti loading python stuff)
ti.block_start_sample = round(ic);
ti.block_end_sample = round(fc);

% Vector of length equal to number of acquisitions, where each entry is the
% number of samples in that acquisition.
acquisition_samples = ti.block_end_sample - ti.block_start_sample;

good_acquisitions = acquisition_samples >= min_acquisition_samples;

ti.block_start_sample = ti.block_start_sample(good_acquisitions);
ti.block_end_sample = ti.block_end_sample(good_acquisitions);

% TODO are these guaranteed to be same length? what does pulsewidth return
% if it stops high (want to err / toss that block in that case)?
ti.block_start_time = ai.time(ti.block_start_sample);
ti.block_end_time = ai.time(ti.block_end_sample);
ti.num_blocks = length(ti.block_end_time);

% TODO TODO fix what caused ti.block_start_time to be empty in 4-10/2/_001 case.
% something seriously wrong with that data?
if numel(ti.block_start_time) > 0
    time = ai.time - ti.block_start_time(1);
else
    error('No pulses detected in scopePin!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
interactive_plots = false;
if interactive_plots
    fsync = figure;
else
    fsync = figure('visible', 'off');
end
plot(time, ai.scopePin);
hold on;
plot(time, ai.olfDispPin);

% So that this script can be run not as a function too, where calls to
% nargin will err.
try
    % TODO handle this better... arguments are getting pretty fragile now
    % and it's not like i'm really taking advantage of the fact they are
    % optional
    switch nargin %#ok
        case 6
            str = sprintf('%s, fly %d: %s', date, fly_num, thorimage_id);
        case 5
            % Assumed the last argument is the date str in this case.
            str = sprintf('%s: %s', date, thorimage_id);
        case 4
            str = sprintf('%s', thorimage_id);
        % TODO empty str basecase
    end
catch err
    assert(strcmp(err.identifier, 'MATLAB:narginout:NoStack'));
    str = '';
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
if write_anything
    print(fsync, fullfile(fig_output_dir, fname), '-dtiffn');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if any(strcmp(fieldnames(ai), 'Frame_Out'))
    frame_out = ai.Frame_Out;
elseif any(strcmp(fieldnames(ai), 'FrameOut'))
    frame_out = ai.FrameOut;
% TODO else error
end

% TODO seems this fails on data from resonant system w/ frame averaging
% (why?) there probably shouldn't be two paths here anyway...
% TODO TODO TODO fix!
% (was hardcoded to true when analyzing data w/o frame averaging from the 
% galvo/galvo system)
only_count_frames_in_scopepin_high = false;
if only_count_frames_in_scopepin_high
    % TODO TODO test this w/ data from resonant scanner system

    assert(length(ti.block_start_sample) == length(ti.block_end_sample));

    frame_out_logical = logical(frame_out);
    % TODO maybe try to share more code w/ below. two ifs? scope shared across
    % them?
    frame_out_diff = diff(frame_out_logical);
    frame_out_rising_edge_samples = find(frame_out_diff > 0);
    frame_out_falling_edge_samples = find(frame_out_diff < 0);

    blocks_rising = [];
    blocks_falling = [];
    for i = 1:length(ti.block_start_sample)
        block_frame_out_rising_edge_samples = frame_out_rising_edge_samples(...
            (frame_out_rising_edge_samples >= ti.block_start_sample(i)) & ...
            (frame_out_rising_edge_samples <= ti.block_end_sample(i)));

        block_frame_out_falling_edge_samples = ...
            frame_out_falling_edge_samples(...
            (frame_out_falling_edge_samples >= ti.block_start_sample(i)) & ...
            (frame_out_falling_edge_samples <= ti.block_end_sample(i)));

        n_rising = length(block_frame_out_rising_edge_samples);
        n_falling = length(block_frame_out_falling_edge_samples);
        if n_falling < n_rising
            assert(n_falling + 1 == n_rising);
            next_fall = find(frame_out_falling_edge_samples > ...
                ti.block_end_sample(i), 1);
            next_fall_sample = frame_out_falling_edge_samples(next_fall);
            block_frame_out_falling_edge_samples = [...
                block_frame_out_falling_edge_samples; next_fall_sample];
        end
        assert(length(block_frame_out_rising_edge_samples) == ...
            length(block_frame_out_falling_edge_samples));

        % TODO or concat along other dimension to get grouped by block
        blocks_rising = [blocks_rising; block_frame_out_rising_edge_samples];
        blocks_falling = [blocks_falling; block_frame_out_falling_edge_samples];
    end

    frame_out_rising_edge_samples = blocks_rising;
    frame_out_falling_edge_samples = blocks_falling;

    last_block_end_index = ti.block_end_sample(end);

    % TODO try to share below code across branches of this if
    frame_out_rising_edge_time = ai.time(frame_out_rising_edge_samples);
    frame_out_falling_edge_time = ai.time(frame_out_falling_edge_samples);
    assert(frame_out_rising_edge_time(1) < frame_out_falling_edge_time(1));
    assert(numel(frame_out_rising_edge_time) == ...
           numel(frame_out_falling_edge_time));
else
    last_block_end_time = ti.block_end_time(end);
    max_trigger_to_last_frame_sec = 1;
    last_fall_check_time = last_block_end_time + max_trigger_to_last_frame_sec;
    frame_out_logical = logical(frame_out(ai.time <= last_fall_check_time));
    % TODO maybe use pulsewidth to be consistent w/ rest of code? why not?
    % or maybe eliminate pulsewidth?
    frame_out_diff = diff(frame_out_logical);
    % TODO TODO also rename 'indices' to 'samples' or 'frames' as appropriate
    frame_out_rising_edge_samples = find(frame_out_diff > 0);
    frame_out_falling_edge_samples = find(frame_out_diff < 0);

    frame_out_rising_edge_time = ai.time(frame_out_rising_edge_samples);
    frame_out_falling_edge_time = ai.time(frame_out_falling_edge_samples);

    % TODO maybe check all? or check samples and check time is sorted?
    % TODO maybe just filter out beginning artifacts if this case is
    % encountered?  frame_out should start low (rise before falling)
    assert(frame_out_rising_edge_time(1) < frame_out_falling_edge_time(1));
    assert(numel(frame_out_rising_edge_time) == ...
           numel(frame_out_falling_edge_time));
    %frame_out_length = frame_out_falling_edge_time - ...
    %                   frame_out_rising_edge_time;

    low2high_times = (frame_out_rising_edge_time(2:end) - ...
                      frame_out_falling_edge_time(1:(end - 1)));

    max_frame_out_low2high_interval = mode(low2high_times) * 1.5;

    falls_to_check = frame_out_falling_edge_time >= last_block_end_time & ...
                     frame_out_falling_edge_time <= last_fall_check_time;
    fall_times_to_check = frame_out_falling_edge_time(falls_to_check);
    fall_indices_to_check = frame_out_falling_edge_samples(falls_to_check);

    % TODO vectorized way to do this?
    last_block_fall_time = nan;
    for i = 1:length(fall_times_to_check)
        fall_time = fall_times_to_check(i);
        max_time = fall_time + max_frame_out_low2high_interval;

        rising_edge_logical = frame_out_rising_edge_time >= fall_time ...
            & frame_out_rising_edge_time <= max_time;

        if sum(rising_edge_logical) == 0
            last_block_fall_time = fall_time;
            last_block_fall_index = fall_indices_to_check(i);
            break
        end
    end

    % TODO why is this getting triggered in the v4 case?
    % (just the stuff that had only 3000 frames because max frames issue?)
    % TODO TODO and why, similarly often, are no scopePin pulses detected?
    if isnan(last_block_fall_time)
        error(['Frame Out continued pulsing right until end of recording. ' ...
              'The recording may have stopped incorrectly.']);
    end

    % TODO right now, this does not include the double pulse i've seen at the
    % last pulse of frame_out in at least one recording. change if necessary.
    last_block_end_index = last_block_fall_index;

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
    plot(plot_times, frame_out_logical(start:stop))
    hold on;
    plot(plot_times, ai.scopePin(start:stop))
    uiwait;
    %}
end

between_frame_dt = diff(frame_out_rising_edge_time);

% TODO why use both this and scopePin?

% Was previously hardcoded to 0.0169 for analysis of data from the resonant
% scanner.
% TODO check this mode is also ~hardcoded value above in resonant data case
consecutive_frame_dt = mode(between_frame_dt);

block_end_frames = find(between_frame_dt > (1.1 * consecutive_frame_dt));

% TODO TODO TODO or check that all blocks are the same length w/in some margin
% of error
block_end_frames = [0 block_end_frames' numel(frame_out_rising_edge_time)];
assert((length(block_end_frames) - 1) == ti.num_blocks);

% TODO is there always a double peak at end of acquistion?
% (high period ~2x as long) is that intentional? should it be counted?
% counted twice?
%

% TODO case where noise peaks are not the last peaks, but are still in 
% block_end_frames, will probably not be handled correctly right now.

% TODO TODO TODO rename to _by_block? was i using it as "by stim"
% recalc if so? need something "by block"?
idx_by_stim = cell(ti.num_blocks, 1);
for i = 1:ti.num_blocks
    idx_by_stim{i} = (block_end_frames(i) + 1):block_end_frames(i + 1);
end
ti.frames_by_stim = idx_by_stim;

av_frames_by_stim = cell(ti.num_blocks, 1);
times_by_stim = cell(ti.num_blocks, 1);
% TODO TODO TODO frame_times used to be same length as movie, but after
% subsetting Frame_Out, it is not. fix frame_times / provide other information
% to match it up to whole movie.
% TODO does this code work in averageNum == 1 case? check
frame_times = [];
for i = 1:ti.num_blocks
   % TODO TODO TODO isn't all the *_by_stim stuff wrong here, since loop was
   % over 1:ti.num_blocks (previously ti.num_trials but same meaning),
   % and n_stim > n_blocks in some cases. change?

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
   disp('size(frame_out_falling_edge_time)')
   disp(size(frame_out_falling_edge_time))
   disp('size(times_by_stim)')
   disp(size(times_by_stim))
   disp('i')
   disp(i)
   disp('max(av_frames_by_stim{i})')
   disp(max(av_frames_by_stim{i}))
   %}
   %
   times_by_stim{i} = frame_out_falling_edge_time(av_frames_by_stim{i});
   % TODO TODO TODO was remy using this frame_times as something other than
   % literally just the times for ***each*** frame? why is this in a loop over
   % either # blocks or # odor presentations then...?
   frame_times = [frame_times; times_by_stim{i}];
end
% TODO TODO TODO does this variable, and how it's computed, even make sense?
ti.av_frames_by_stim = av_frames_by_stim;
ti.times_by_stim = times_by_stim;

if debug
    disp(sprintf('length(frame_times): %d', length(frame_times)));
    disp(sprintf('ti.num_frames: %d', ti.num_frames));
end
assert(length(frame_times) == ti.num_frames);
ti.frame_times = frame_times;

% TODO TODO also rename these to blocks
%% calculate trial_start and trial_end frames
% TODO TODO TODO TODO fix... seems wrong in 2019-07-25/2/_007 case
% where both of these are [1,2,3]
ti.trial_end = cumsum(cellfun(@numel, ti.av_frames_by_stim'));
ti.trial_start = [0 ti.trial_end(1:end-1)] + 1;

%% calculate frames when olfDispPin is high
% TODO rename ic / fc to not shadow earlier def for debugging
[~, ic, fc, ~] = pulsewidth(ai.olfDispPin(1:last_block_end_index));
ti.stim_ic_idx = round(ic);
ti.stim_fc_idx = round(fc);
ti.stim_ict = ai.time(ti.stim_ic_idx);
ti.stim_fct = ai.time(ti.stim_fc_idx);
% # of odor pulses
ti.num_stim = length(ti.stim_fct);

start_frame = zeros(1, ti.num_blocks);
end_frame = zeros(1, ti.num_blocks);

% TODO TODO should this just be like all frames? ai.time?
% cause right now, frame_times seems to be of length equal to # blocks...
for i = 1:ti.num_stim
    idx = find(ti.frame_times >= ti.stim_ict(i) & ...
        ti.frame_times <= ti.stim_fct(i));

    start_frame(i) = idx(1);
    end_frame(i) = idx(end);
end
ti.stim_on = start_frame;
ti.stim_off = end_frame;

% TODO delete this after updating downstream code to use new names
ti.block_ic_idx = ti.block_start_sample;
ti.block_fc_idx = ti.block_end_sample;
%

if write_anything 
    if create_mat
        save(mat_filepath, 'ti');
    else
        save(mat_filepath, 'ti', '-append');
    end
    updated = true;
end

end
