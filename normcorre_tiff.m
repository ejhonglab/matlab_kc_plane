
function [rig_updated, nr_updated] = normcorre_tiff(input_tif_path, output_dir)

% TODO delete
%{
raw_dir = '/mnt/nas/mb_team/raw_data';
analysis_dir = '/mnt/nas/mb_team/analysis_output';
fly_dir = '2019-07-25/2';
input_tif_path = fullfile(raw_dir, fly_dir, 'tif_stacks', '_008.tif');
% TODO just get this to auto-increment if exists
% TODO write params (named similarly) alongside tiff?
num = 5;
output_dir = fullfile(analysis_dir, fly_dir);
%}

rig_updated = false;
nr_updated = false;

% TODO gcp ("get current parallel pool") necessary?
%gcp;

% TODO TODO need to make this dir first? if that why some of the tiff saving
% earlier seemed to fail?
output_subdir = 'tif_stacks';

% Gets only the final directory in the path, whether or not there is a trailing
% filesep character.
[~, thorimage_id, ~] = fileparts(input_tif_path);

% TODO join common prefix once first
rig_tif = fullfile(output_dir, output_subdir, [thorimage_id '_rig.tif']);
avg_rig_tif = fullfile(output_dir, output_subdir, 'AVG', 'rigid', ...
    ['AVG' thorimage_id '_rig.tif']);

% TODO delete suffixed int for playing w/ params
%nr_suf = sprintf('_nr%d.tif', num);
nr_suf = '_nr.tif';
nr_tif = fullfile(output_dir, output_subdir, [thorimage_id nr_suf]);
avg_nr_tif = fullfile(output_dir, output_subdir, 'AVG', 'nonrigid',...
    ['AVG' thorimage_id nr_suf]);

need_rig_tif = ~ exist(rig_tif, 'file');
need_avg_rig_tif = ~ exist(avg_rig_tif, 'file');
need_nr_tif = ~ exist(nr_tif, 'file');
need_avg_nr_tif = ~ exist(avg_nr_tif, 'file');

if ~ (need_rig_tif | need_avg_rig_tif | need_nr_tif | need_avg_nr_tif)
    disp('All registration already done.');
    return;
end

% TODO make a script that looks over a certain number of names, calling the
% new one-folder function on each
% start of old for loop over ThorImage folders
% TODO reformat with more information
disp(['normcorre: ' thorimage_id]);

interactive_plots = false;

try
    % Remy: this seems like it might just be reading in the first frame?
    %%%Y = input_tif_path;
    % TODO TODO also only read this if at least one motion correction would be
    % run
    Y = imread_big(input_tif_path);

    % rigid moco (normcorre)
    % TODO just pass filename instead of Y, and compute dimensions or whatever
    % separately, so that normcorre can (hopefully?) take up less memory
    if need_rig_tif
        MC_rigid = MotionCorrection(Y);

        options_rigid = NoRMCorreSetParms('d1',MC_rigid.dims(1),...
            'd2',MC_rigid.dims(2),...
            'bin_width',50,...
            'max_shift',15,...
            'phase_flag', 1,...
            'us_fac', 50,...
            'init_batch', 100,...
            'plot_flag', false,...
            'iter', 2); 

        % TODO so is nothing actually happening in parallel?
        %% rigid moco
        MC_rigid.motionCorrectSerial(options_rigid);  % can also try parallel
        MC_rigid.computeMean();
        MC_rigid.correlationMean();
        MC_rigid.crispness();
        disp('normcorre done');

        %% plot shifts
        if interactive_plots
            figure;
        else
            figure('visible', 'off');
        end
        
        subplot(2,1,1)
        plot(MC_rigid.shifts_x);
        title('MC_rigid.shifts_x', 'Interpreter', 'None');
        subplot(2,1,2)
        plot(MC_rigid.shifts_y);
        title('MC_rigid.shifts_y', 'Interpreter', 'None');

        % TODO add date + fly_num for this too? just include last bit of path?
        str = sprintf('%s: rigid normcorre', thorimage_id);
        tt = suptitle(str);
        tt.Interpreter = 'None';
        %%
        if interactive_plots
            figure;
        else
            figure('visible', 'off');
        end
        subplot(2,1,1)
        % TODO TODO fix harcoded indices
        %%%%%plot(MC_rigid.corrY(1:2992));
        plot(MC_rigid.corrY);
        title('MC_rigid.corrY', 'Interpreter', 'None');
        subplot(2,1,2)
        %%%%%plot(MC_rigid.corrM(1:2992));
        plot(MC_rigid.corrM);
        title('MC_rigid.corrM', 'Interpreter', 'None');

        % TODO add date + fly_num for this too?
        str = sprintf('%s: corr. results of rigid normcorre', thorimage_id);
        tt = suptitle(str);
        tt.Interpreter = 'None';
        %%
        % save .tif
        M = MC_rigid.M;
        M = uint16(M);  
        tiffoptions.overwrite = true;

        disp(['saving tiff to ' rig_tif]);
        saveastiff(M, rig_tif, tiffoptions);
    end
    rig_updated = true;

    if need_avg_rig_tif
        %%
        % save average image
        %AVG = single(mean(MC_rigid.M,3));
        AVG = single(MC_rigid.template);
        tiffoptions.overwrite = true;

        disp(['saving tiff to ' avg_rig_tif]);
        saveastiff(AVG, avg_rig_tif, tiffoptions);
    end

    %%    
    % save MC_rigid to .mat file
    % MC_rigid_copy = MC_rigid;
    % 
    % MC_rigid.file = [];
    % MC_rigid.Y = [];
    % MC_rigid.M = [];
    % 
    % mat_filepath = fullfile(output_dir, 'mat_files', [thorimage_id '.mat']);
    % save(mat_filepath, 'MC_rigid');
    % disp('normcorre done')
    % disp('.mat saved');

    %% pw-rigid motion correction (in parallel)
    % end of old for loop over ThorImage folders

    % TODO maybe del all variables above that are no longer needed, to avoid
    % memory errors

    if need_nr_tif
        MC_nonrigid = MotionCorrection(Y);
        % TODO maybe try upping max_dev for bad warping?
        % or decreasing size of regions over which rigi stuff calced?
        options_nonrigid = NoRMCorreSetParms('d1',MC_nonrigid.dims(1),...
            'd2',MC_nonrigid.dims(2),...
            'grid_size',[64,64],...
            'mot_uf',4,...
            'bin_width',50,...
            'max_shift',[15 15],...
            'max_dev',3,...
            'us_fac',50,...
            'init_batch',200,...
            'iter', 2 ...
        );

        MC_nonrigid.motionCorrectParallel(options_nonrigid);
        MC_nonrigid.computeMean();
        MC_nonrigid.correlationMean();
        MC_nonrigid.crispness();
        disp('non-rigid normcorre done')

        %%
        if interactive_plots
            figure;
        else
            figure('visible', 'off');
        end
            subplot(2,1,1)
            plot(MC_nonrigid.shifts_x);
            title('MC_nonrigid.shifts_x', 'Interpreter', 'None');
            subplot(2,1,2)
            plot(MC_nonrigid.shifts_y);
            title('MC_nonrigid.shifts_y', 'Interpreter', 'None');
            
            str = sprintf('%s: non-rigid normcorre', thorimage_id);
            tt = suptitle(str);
            tt.Interpreter = 'None';
        %%
        if interactive_plots
            figure;
        else
            figure('visible', 'off');
        end
            subplot(2,1,1)
            plot(MC_nonrigid.corrY);
            title('MC_nonrigid.corrY', 'Interpreter', 'None');
            subplot(2,1,2)
            plot(MC_nonrigid.corrM);
            title('MC_nonrigid.corrM', 'Interpreter', 'None');
            
            str = sprintf('%s: corr. results of non-rigid normcorre', ...
                thorimage_id);

            tt = suptitle(str);
            tt.Interpreter = 'None';
        %%
        % save .tif
        M = uint16(MC_nonrigid.M);
        tiffoptions.overwrite  = true;
        disp(['saving tiff to ' nr_tif]);
        saveastiff(M, nr_tif, tiffoptions);
    end
    nr_updated = true;

    % TODO need to escape underscore in figure titles so that it doesnt make
    % next number subscript (and not displaying the underscore).
    % actually it doesn't always seem to be a problem... why?

    % TODO should really just recompute these at least whenever other thing
    % changes. (what happens in ~need_tif & need_avg case now?)
    if need_avg_nr_tif
        % TODO flag to disable saving this average
        %AVG = single(mean(MC_nonrigid.M,3));
        AVG = single(MC_nonrigid.template);
        tiffoptions.overwrite = true;
        disp(['saving tiff to ' avg_nr_tif]);
        saveastiff(AVG, avg_nr_tif, tiffoptions);
    end

catch err
    if strcmp(err.identifier, 'MATLAB:nomem')
        % To understand memory usage of variables *before* offending call.
        % TODO TODO fix display here... if this is showing variables, it's not
        % showing bytes
        % output once seems to have been
        % "name,size,bytes,class,global,sparse,complex,nesting,persistent", each
        % on one line
        disp(whos);
    end
    rethrow(err);
end

end
