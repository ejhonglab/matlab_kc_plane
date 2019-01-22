%%
clear, clc;
gcp;
%%
date = '2019-01-18';

folder = sprintf(['/media/remy/remy-storage/Remy''s Dropbox Folder/',...
    'HongLab @ Caltech Dropbox/Remy/2019 data analysis/%s_analysis'], date);

nams = {'_001', '_003', '_006', '_007', '_008'};
is_memmaped = 0;

%%

for ni = 1:numel(nams)
    nam = nams{ni};
    disp(['normcorre: ' nam])
    
    if is_memmaped
        filename = [nam '_rig.mat'];    
        filepath = fullfile(outfolder, 'tif stacks', filename);
        data = matfile(filepath,'Writable',true);

    else
        filename = [nam '.tif'];    
        filepath = fullfile(folder, 'tif stacks', filename);
    end
    
    %%
    Y = imread_big(filepath);
    ts = size(Y,3);
    %%
    % rigid moco (normcorre)
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
    %% rigid moco
    MC_rigid.motionCorrectSerial(options_rigid);  % can also try parallel
    MC_rigid.computeMean();
    MC_rigid.correlationMean();
    MC_rigid.crispness();
    disp('normcorre done')
    %% plot shifts
    figure
    subplot(2,1,1)
    plot(MC_rigid.shifts_x);
    title('MC_rigid.shifts_x', 'Interpreter', 'None');
    subplot(2,1,2)
    plot(MC_rigid.shifts_y);
    title('MC_rigid.shifts_y', 'Interpreter', 'None');
    
    str = sprintf('%s: rigid normcorre', nam);
    tt = suptitle(str);
    tt.Interpreter = 'None';
    %%
    figure
    subplot(2,1,1)
    plot(MC_rigid.corrY(1:2992));
    title('MC_rigid.corrY', 'Interpreter', 'None');
    subplot(2,1,2)
    plot(MC_rigid.corrM(1:2992));
    title('MC_rigid.corrM', 'Interpreter', 'None');
    
    str = sprintf('%s: corr. results of rigid normcorre', nam);
    tt = suptitle(str);
    tt.Interpreter = 'None';
%%
    % save .tif
    M = MC_rigid.M;
    M = uint16(M);  
    tiffoptions.overwrite = true;
    saveastiff(M, fullfile(folder, 'tif stacks', [nam '_rig.tif']), tiffoptions);
%%
    % save average image
    %AVG = single(mean(MC_rigid.M,3));
    AVG = single(MC_rigid.template);
    tiffoptions.overwrite = true;
    saveastiff(AVG, fullfile(folder, 'tif stacks', 'AVG', 'rigid', ['AVG' nam '_rig.tif']), tiffoptions);
    %%    
    % save MC_rigid to .mat file
%     MC_rigid_copy = MC_rigid;
% 
%     MC_rigid.file = [];
%     MC_rigid.Y = [];
%     MC_rigid.M = [];
% 
%     mat_filepath = fullfile(outfolder, 'mat files', [nam '.mat']);
%     save(mat_filepath, 'MC_rigid');
%     disp('normcorre done')
%     disp('.mat saved');
end

%% pw-rigid motion correction (in parallel)

MC_nonrigid = MotionCorrection(Y);
options_nonrigid = NoRMCorreSetParms('d1',MC_nonrigid.dims(1),...
    'd2',MC_nonrigid.dims(2),...
    'grid_size',[64,64],...
    'mot_uf',4,...
    'bin_width',50,...
    'max_shift',[15 15],...
    'max_dev',3,...
    'us_fac',50,...
    'init_batch',200,...
    'iter', 2);

MC_nonrigid.motionCorrectParallel(options_nonrigid);
MC_nonrigid.computeMean();
MC_nonrigid.correlationMean();
MC_nonrigid.crispness();
disp('non-rigid normcorre done')
%%
figure
    subplot(2,1,1)
    plot(MC_nonrigid.shifts_x);
    title('MC_nonrigid.shifts_x', 'Interpreter', 'None');
    subplot(2,1,2)
    plot(MC_nonrigid.shifts_y);
    title('MC_nonrigid.shifts_y', 'Interpreter', 'None');
    
    str = sprintf('%s: non-rigid normcorre', nam);
    tt = suptitle(str);
    tt.Interpreter = 'None';
%%
figure
    subplot(2,1,1)
    plot(MC_nonrigid.corrY);
    title('MC_nonrigid.corrY', 'Interpreter', 'None');
    subplot(2,1,2)
    plot(MC_nonrigid.corrM);
    title('MC_nonrigid.corrM', 'Interpreter', 'None');
    
    str = sprintf('%s: corr. results of non-rigid normcorre', nam);
    tt = suptitle(str);
    tt.Interpreter = 'None';
%%
% save .tif
M = uint16(MC_nonrigid.M);
tiffoptions.overwrite = true;
saveastiff(M, fullfile(folder, 'tif stacks', [nam '_nr.tif']), tiffoptions);
%%
    %AVG = single(mean(MC_nonrigid.M,3));
    AVG = single(MC_nonrigid.template);
    tiffoptions.overwrite = true;
    saveastiff(AVG, fullfile(folder, 'tif stacks', 'AVG', 'nonrigid', ['AVG' nam '_nr.tif']), tiffoptions);

