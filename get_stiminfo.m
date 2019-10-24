
function [updated] = get_stiminfo(thorimage_dir, thorsync_dir, ...
    output_dir, update)

% Output values used by populate_db.py/gui.py (as of 2019-08-02):
% - frame_times
% - block_start_frame
% - block_end_frame
% - block_start_sample
% - stim_on
% - stim_off

if nargin == 0
    close all;
    clear all; %#ok<CLALL>
    matlab_testing = true;
else
    % TODO some generaly MATLAB way to check we have all args, without having to
    % hardcode # of args?
    full_n_args = 4;
    narginchk(full_n_args, full_n_args);
    matlab_testing = false;
end

% Set to true to get arguments for matlab_testing case from an external
% invocation of this function.
print_inputs = false;
if matlab_testing
    write_anything = false;
    interactive_plots = true;
    update = true;
    print_inputs = false;
    
    %{
    % Known good v4:
    thorimage_dir='/mnt/nas/mb_team/raw_data/2019-08-27/9/fn_0002';
    thorsync_dir='SyncData002';
    output_dir='/mnt/nas/mb_team/analysis_output/2019-08-27/9';
    %}
    %{
    % Max num frames hit in each block (v4, frame averaging):
    thorimage_dir='/mnt/nas/mb_team/raw_data/2019-10-04/1/fn_0002';
    thorsync_dir='SyncData002';
    output_dir='/mnt/nas/mb_team/analysis_output/2019-10-04/1';
    %}
    
    % TODO this is failing somehow in 2019-01-18/2/* cases figure out why. may
    % impact future / other current experiments. (still failing?)
else
    write_anything = true;
    interactive_plots = false;
    % TODO maybe somehow automatically print inputs on error
    % (something like atexit / some global catch-all error handler?)
    % builtin stuff (MException.last) seems to not support being called
    % from within a function...
    % there is: www.mathworks.com/matlabcentral/fileexchange/\
    % 15059-handleerror-generic-error-handling-function , but I'm not sure
    % how it works / don't want the extra dependency
    % could set some global flags or something? possible w/o global?
end

if print_inputs
    disp('');
    fprintf("thorimage_dir='%s';\n", thorimage_dir);
    fprintf("thorsync_dir='%s';\n", thorsync_dir);
    fprintf("output_dir='%s';\n", output_dir);
    fprintf("update=%s;\n", string(update));
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
    if (write_anything && (~update) && ...
        any(contains(matwho(mat_filepath), 'ti')))
    
        % TODO put behind verbose flag?
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

xml_err_msg = 'Unexpected XML format.';
% <v4 XML format
if isstruct(ThorImageExperiment.Streaming)
    streaming = ThorImageExperiment.Streaming;
% v4+ XML format
elseif iscell(ThorImageExperiment.Streaming)
    streaming = ThorImageExperiment.Streaming{1};
else
    error(xml_err_msg);
end
% This is the number of frames after any averaging (looking at some v4
% data w/ averaging and the TIFF output).
num_frames = str2double(streaming.Attributes.frames);
% TODO need to subtract flyback frames to get effective #?
% or done automatically?

% It seemed at one point like Thor's software was wrong about amount of
% time that could be aquired at a given value of this variable,
% since it didn't factor in averaging (still true? report?).
stim_max_frames = str2double(streaming.Attributes.stimulusMaxFrames);

% <v4 XML format
if isstruct(ThorImageExperiment.LSM)
    lsm = ThorImageExperiment.LSM;
% v4+ XML format
elseif iscell(ThorImageExperiment.LSM)
    lsm = ThorImageExperiment.LSM{1};
else
    error(xml_err_msg);
end
averageMode = str2double(lsm.Attributes.averageMode);
averageNum = str2double(lsm.Attributes.averageNum);
fr = str2double(lsm.Attributes.frameRate);

if averageMode == 0
    averageNum = 1;
elseif averageMode == 1
   fr = fr / averageNum;
   stim_max_frames = stim_max_frames / averageNum;
end

% From these assertions we establish that ai.time:
% - starts at time 0
% - is in units of seconds
% - all frames could be captured in its range
assert(min(ai.time) == 0);
% (If ThorImage were acquiring frames at frame rate (fr) for the entire
% length of time contained in the ThorSync data.)
max_possible_n_frames = max(ai.time) * fr;
assert(num_frames <= max_possible_n_frames);
% Partially to establish that ai.time is in units of seconds, but also 
% checks an assumption that at least some fraction of time is spent taking
% frames.
% Most experiments have this value below ~0.4, but some bad ones have ~0.7
% (like 2019-10-04/4/fn_0004). Could still warn if > ~0.4?
min_time_frac_acquiring_frames = 0.8;
assert((max(ai.time) * fr - num_frames) / num_frames <= ...
    min_time_frac_acquiring_frames);

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

ti.block_start_sample = round(ic);
block_end_sample = round(fc);

% Vector of length equal to number of acquisitions, where each entry is the
% number of samples in that acquisition.
acquisition_samples = block_end_sample - ti.block_start_sample;

good_acquisitions = acquisition_samples >= min_acquisition_samples;

ti.block_start_sample = ti.block_start_sample(good_acquisitions);
block_end_sample = block_end_sample(good_acquisitions);

% TODO are these guaranteed to be same length? what does pulsewidth return
% if it stops high (want to err / toss that block in that case)?
block_start_time = ai.time(ti.block_start_sample);
block_end_time = ai.time(block_end_sample);
num_blocks = length(block_end_time);

% TODO TODO fix what caused block_start_time to be empty in 4-10/2/_001 case.
% something seriously wrong with that data?
if numel(block_start_time) > 0
    time = ai.time - block_start_time(1);
else
    error('No pulses detected in scopePin!');
end

allow_reaching_stim_max_frames = true;
max_stimulus_frames_case = false;
% In general, would want to compare each block's # of frames to this,
% in case some blocks are longer than others.
if stim_max_frames * num_blocks == num_frames
    err_msg = ['Max number of stimulus frames reached in each block! ' ...
        'Some odor presentations may have been missed!'];
    
    if allow_reaching_stim_max_frames
        warning(err_msg);
        max_stimulus_frames_case = true;
    else
        error(err_msg); %#ok<*UNRCH>
    end
end

if any(strcmp(fieldnames(ai), 'Frame_Out'))
    frame_out = ai.Frame_Out;
elseif any(strcmp(fieldnames(ai), 'FrameOut'))
    frame_out = ai.FrameOut;
else
    error(xml_err_msg);
end
frame_out_logical = logical(frame_out);

if interactive_plots
    fsync = figure;
else
    fsync = figure('visible', 'off');
end
plot(time, ai.scopePin);
hold on;
plot(time, ai.olfDispPin);
% To check visually that that frames are captured throughout the full
% extent of each block.
plot(time, frame_out_logical);

parts = split(thorimage_dir, filesep);
date = char(parts(end-2));
fly_num = char(parts(end-1));

% The last two args ('Interpreter', 'none') prevent underscores from just 
% subscripting the next character, rather than displaying.
title(sprintf('%s, fly %s: %s', date, fly_num, thorimage_id), ...
    'Interpreter', 'none');
legend('olfDispPin', 'scopePin', 'frame_out_logical', ...
    'Interpreter', 'none');
ylim([-1 6]);

% TODO also plot mirror signal (or maybe just check it should have time to
% switch at start of scopePin pulses)?

fig_output_dir = fullfile(output_dir, 'figures');
if ~exist(fig_output_dir, 'dir')
    mkdir(fig_output_dir);
end

fname = sprintf('%s_trial_structure.tif', thorimage_id);
if write_anything
    print(fsync, fullfile(fig_output_dir, fname), '-dtiffn');
end

assert(length(ti.block_start_sample) == length(block_end_sample));

% TODO seems this fails on data from resonant system w/ frame averaging
% (why?) there probably shouldn't be two paths here anyway...
% TODO TODO TODO fix!
% (was hardcoded to true when analyzing data w/o frame averaging from the 
% galvo/galvo system)

% 2019-10-24: and when flag was false, failed on all my 10-4 data from 
% resonant system, w/ "Frame Out continued pulsing..." error

% while with this flag true on the same data:
%  seemed to work (no errors here):
%  - [2019-10-04/]1/SyncData001
%  - 2/SyncData002
%  - 4/SyncData004 (fn_0003 here, apparent num mismatch)
%  failed:
%  - 1/SyncData002
%  - 2/SyncData001
%  - 4/SyncData005 (fn_0004 here, apparent num mismatch)
% This above that failed did so with this error:
% Index exceeds the number of array elements (0).
% start_frame(i) = idx(1); (changed this part so not going to be exactly
% the same error anymore)

only_count_frames_in_scopepin_high = false;
% TODO TODO test this first branch w/ data from resonant scanner system
if only_count_frames_in_scopepin_high
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
            (frame_out_rising_edge_samples <= block_end_sample(i)));

        block_frame_out_falling_edge_samples = ...
            frame_out_falling_edge_samples(...
            (frame_out_falling_edge_samples >= ti.block_start_sample(i)) & ...
            (frame_out_falling_edge_samples <= block_end_sample(i)));

        n_rising = length(block_frame_out_rising_edge_samples);
        n_falling = length(block_frame_out_falling_edge_samples);
        if n_falling < n_rising
            assert(n_falling + 1 == n_rising);
            next_fall = find(frame_out_falling_edge_samples > ...
                block_end_sample(i), 1);
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

    last_block_end_index = block_end_sample(end);

    % TODO try to share below code across branches of this if
    frame_out_rising_edge_time = ai.time(frame_out_rising_edge_samples);
    frame_out_falling_edge_time = ai.time(frame_out_falling_edge_samples);
    assert(frame_out_rising_edge_time(1) < frame_out_falling_edge_time(1));
    assert(numel(frame_out_rising_edge_samples) == ...
           numel(frame_out_falling_edge_samples));
else
    last_block_end_time = block_end_time(end);
    max_trigger_to_last_frame_sec = 1;
    last_fall_check_time = last_block_end_time + max_trigger_to_last_frame_sec;
    % TODO TODO TODO is this the line causing "pulsing until end of
    % recording issue"?
    frame_out_logical = frame_out_logical(ai.time <= last_fall_check_time);
    % TODO maybe use pulsewidth to be consistent w/ rest of code? why not?
    % or maybe eliminate pulsewidth?
    frame_out_diff = diff(frame_out_logical);
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

    % Time differences in seconds (since ai.time is in seconds).
    low2high_times = (frame_out_rising_edge_time(2:end) - ...
                      frame_out_falling_edge_time(1:(end - 1)));

    max_frame_out_low2high_interval = mode(low2high_times) * 1.5;

    falls_to_check = frame_out_falling_edge_time >= last_block_end_time & ...
                     frame_out_falling_edge_time <= last_fall_check_time;
    
    % TODO TODO TODO fix!! probably need base case when no falls_to_check
    % sum(falls_to_check) == 0, but is that because no blocks are detected
    % or because no frames seem to be taken after scopePin goes low?
    % and if the latter, is it *true* that no frame are taken after scope
    % pin goes low? am i counting the wrong edge or something?
    
    fall_times_to_check = frame_out_falling_edge_time(falls_to_check);
    
    % TODO maybe move this check earlier. most of stuff inside this
    % top-level else may not need to be run here?
    if max_stimulus_frames_case
        last_block_end_index = block_end_sample(end);
    else
        % This ever happen outside of when max stimulus frames is reached?
        if isempty(fall_times_to_check)
            error('No frame pulses when scopePin was low, despite expectation.');
        end

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

        % TODO get a test case for when this should raise
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
end

% Need to look at frame out information in addition to scope acquisition 
% trigger pin (what we call "scopePin"), since in some cases, the two
% don't match up exactly (like if the max number of stimulus frames is
% reached, while the scopePin is still high).

between_frame_dt = diff(frame_out_rising_edge_time);
% Was previously hardcoded to 0.0169 for analysis of data from the resonant
% scanner.
% TODO check this mode is also ~hardcoded value above in resonant data
% case?
consecutive_frame_dt = mode(between_frame_dt);
block_end_nonavgd_frames = find(...
    between_frame_dt > (1.1 * consecutive_frame_dt));

% TODO TODO check that all blocks are the same length w/in some margin
% of error?
block_end_nonavgd_frames = [0 block_end_nonavgd_frames' ...
    numel(frame_out_rising_edge_time)];
assert((length(block_end_nonavgd_frames) - 1) == num_blocks);

%floor(block_end_nonavgd_frames / ti.averageNum)

% TODO is there always a double peak at end of acquistion?
% (high period ~2x as long) is that intentional? should it be counted?
% counted twice?

% TODO case where noise peaks are not the last peaks, but are still in 
% block_end_nonavgd_frames, will probably not be handled correctly right now.

% TODO TODO rename to _by_block? was i using it as "by stim"
% recalc if so? need something "by block"?
frames_by_stim = cell(num_blocks, 1);
for i = 1:num_blocks
    frames_by_stim{i} = (block_end_nonavgd_frames(i) + 1 ...
        ):block_end_nonavgd_frames(i + 1);
end

av_frames_by_stim = cell(num_blocks, 1);
times_by_stim = cell(num_blocks, 1);
% TODO TODO TODO frame_times used to be same length as movie, but after
% subsetting Frame_Out, it is not. fix frame_times / provide other information
% to match it up to whole movie.
% TODO does this code work in averageNum == 1 case? check
frame_times = [];
for i = 1:num_blocks
   % TODO rename idx
   idx = frames_by_stim{i};
   % TODO TODO TODO are frameout pulses always in number that is a multiple
   % of averageNum??? need to handle case where there might be spurious
   % extra frames?
   % (pretty sure not always multiple. just floor() it, probably)
   % TODO is averageNum a correct base case?
   % and isn't idx indexed as in ThorSync (with the much higher 30KHz), or
   % is it not?
   % TODO TODO TODO does this variable, and how it's computed, make sense?
   av_frames_by_stim{i} = idx(averageNum):averageNum:idx(end);

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

assert(length(frame_times) == num_frames);
ti.frame_times = frame_times;

% TODO TODO TODO fix... seems wrong in 2019-07-25/2/_007 case
% where both of these are [1,2,3] (still?)
ti.block_end_frame = cumsum(cellfun(@numel, av_frames_by_stim'));
ti.block_start_frame = [0 ti.block_end_frame(1:end-1)] + 1;

%% calculate frames when olfDispPin is high
[~, ic, fc, ~] = pulsewidth(ai.olfDispPin(1:last_block_end_index));
stim_start_time = ai.time(round(ic));
stim_stop_time = ai.time(round(fc));
% # of odor pulses
num_stim = length(stim_stop_time);

start_frame = [];
end_frame = [];

% TODO TODO should this just be like all frames? ai.time?
% cause right now, frame_times seems to be of length equal to # blocks...
for i = 1:num_stim
    stim_frames = find(ti.frame_times >= stim_start_time(i) & ...
        ti.frame_times <= stim_stop_time(i));
    
    if isempty(stim_frames)
        % Otherwise there should be frames for the time range above.
        assert(max_stimulus_frames_case);

        warning(sprintf('Odor presentation %d/%d has no frames!!!', ...
            i, num_stim));

        continue;
    end
    
    start_frame = [start_frame; stim_frames(1)];
    end_frame = [end_frame; stim_frames(end)];
end
ti.stim_on = start_frame;
ti.stim_off = end_frame;

assert(mod(length(start_frame), num_blocks) == 0);
assert(length(start_frame) == length(end_frame));
% TODO if returning stim_on/stim_off less than number of odor
% pulses, maybe also check:
% 1) equal number in each block (not just that their lengths are divisible
%    by num_blocks
% 3) (?) they are only included if the full trial is there?

% TODO checks that frame_times are correct even in case where max
% num stimulus frames is reached well before end of recording

if write_anything
    if create_mat
        save(mat_filepath, 'ti');
    else
        save(mat_filepath, 'ti', '-append');
    end
    updated = true;
end

if matlab_testing
    % Because things will still go out of scope at the end of the function
    % with no arguments, otherwise.
    keyboard;
end

end
