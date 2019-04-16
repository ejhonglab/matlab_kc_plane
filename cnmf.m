
% TODO any conflicts w/ installed CaImAn? rename?
function cnmf(tiff_path, output_dir)
% TODO maybe take some of cnmf params too?

% Delete these harcoded paths once I start using this as a function, and
% not interactively.
tiff_path = '/mnt/nas/mb_team/analysis_output/2019-01-18/2/tif_stacks/_003_nr.tif';
output_dir = '/mnt/nas/mb_team/analysis_output/2019-01-18/2';
%

output_dir = fullfile(output_dir, 'cnmf');

% Gets only the final directory in the path, whether or not there is a trailing
% filesep character.
[~, tiff_prefix, ~] = fileparts(tiff_path);

assert(tiff_prefix(1) == '_');
% Will need to change this if our ThorImage output naming convention changes.
thorimage_id = tiff_prefix(1:4);

%filename = [thorimage_id '_nr.tif'];

% TODO TODO make new tif_stacks in output dir for output of normcorre?
% can't be relative to output_dir if only tif_stacks is going to live alongside
% the raw data...
%filepath = fullfile(output_dir, 'tif_stacks', filename);

% TODO TODO why is tiff loaded twice? first time just to compute d? Y doesn't
% seem to be used... is d?
Y = loadtiff(tiff_path);

if ndims(Y) == 4
    % Last dimension is # of timepoints (T)
    % Second to last dimension is d3 (Z?)
    [d1,d2,~,~] = size(Y);                            % dimensions of dataset
else
    [d1,d2,~] = size(Y);
end


CNM = CNMF;

% K gives number of cells of each expected size.
K = 550;
% And tau gives the sizes of each (along each of X and Y dims)
% (actually, non-scalar tau breaks the ability to add components in
% edit_footprints)
tau = 2; %3;

% Order of autoregressive system
% (p=0 no dynamics, p=1 just decay, p=2, both rise and decay)

% TODO these parameters depend on frame rate / image size / noise / etc
% one setting doesnt tend to work globally
% most sensitive: tau, merge_thr, energy_thresh, thr_method

p = 0;
% Other values: 0.95
% TODO is this only used by merge function?
merge_thr = 0.85;                                  % merging threshold

options = CNMFSetParms(...   
    'd1',d1,'d2',d2,...                         % dimensions of datasets
    'fr', 4.55,...
    'p',p,...                                   % order of AR dynamics    
    'gSig',tau,...                                % half size of neuron
    'merge_thr',merge_thr,...                        % merging threshold  
    'nb', 2,...                                  % number of background components    
    'min_SNR',3,...                             % minimum SNR threshold
    % TODO seems to correspond to quality/rval_thr in current python
    'space_thresh',0.3,...                      % space correlation threshold
    % TODO TODO does this correspond to anything in current python? merge/thr?
    'time_thresh',0.3,...                      % space correlation threshold
    % TODO = min_cnn_thr? or cnn_lowest?
    % TODO how do those two interact?
    'cnn_thr',0.2,...                            % threshold for CNN classifier    
    % TODO this is init/method_init, right?
    % TODO greedy_roi or greedy_pnr?
    'init_method', 'greedy',...
    % TODO seems to be spatial/method_exp
    'search_method', 'ellipse',...
    % spatial
    'thr_method', 'nrg',...
    % spatial
    'nrgthr', .99,...
    % TODO TODO this not in python version? it get removed / diff name / not
    % impl? if latter, implement?
    'min_pixel', 10,...
    % TODO python equiv?
    'noise_norm', true,...
    % TODO why is this false?
    % temporal
    'bas_nonneg', false,...
    % TODO python equiv?
    'cont_threshold', .95);
%'cluster_pixels', true,...   
%'dist', 3,...
%'max_size', 3,...
%'min_size', 1,...

%
% Below is the standard processing pipeline. This processing can be
% executed in one shot using the CNM.fit function:
%CNM.fit(Y,options,K)

% load the dataset and create the object
CNM.readFile(tiff_path);                         % insert path to file here  
CNM.optionsSet(options);                        % setup the options structure
CNM.gSig = options.gSig;

% Process the dataset
% TODO is this not supposed to be a function call?
CNM.preprocess;             % preprocessing (compute some quantities)
CNM.initComponents(K);      % initialization 
CNM.plotCenters()           % plot center of ROIs detected during initialization

% TODO why plotCenters before any updates?
% POSTPROCESSING BEGINS
% Update components
CNM.updateSpatial();        % update spatial components
CNM.updateTemporal(0);      % update temporal components (do not deconvolve at this point)

%% Plot contours
% TODO any point to this cell? not information we'd get in manual
% refinement process?
% Apparently, at least CNM.contours needs to be set to empty for contours
% to update from some other state variable (according to Remy)?
CNM.cm = [];
CNM.contours = [];
CNM.COM();
CNM.K = size(CNM.A,2);
CNM.CI = [];
%CNM.correlationImage();

% TODO is this figure call necessary or not?
%figure
% TODO what is 'level' when hovering over contours?
CNM.plotContours();

%% If you want to edit components, run this cell

% TODO flag to skip, among w/ other interactive stuff
%CNM.manuallyRefineComponents();
CNM.updateSpatial();        % update spatial components
CNM.updateTemporal(0);      % update temporal components (do not deconvolve at this point)


%% component classification
% Remy generally just uses this for 3D. Or if lots of stuff is
% oversegmented.
%{
CNM.evaluateComponents();   % evaluate spatial components based on their correlation with the data
CNM.eventExceptionality();  % evaluate traces

ind_throw = find(~CNM.keep_eval);

ind_keep = find(CNM.keep_eval);

figure
CNM.plotContours(1,[],[],ind_keep); % plot components to get rid of
%

% TODO this keep_eval only comes from evaluateComponents, right?
%get rid of components if sure
CNM.keepComponents(CNM.keep_eval);      % keep the components that are above certain thresholds
%}

%% merge components
% TODO should this also be done before first refinement?
% TODO i forget, was remy saying this needed to be done before each
% plotContours or something? (just put it before that?)
% R: "never actually does anything" "maybe if really bad?"
CNM.merge();                            %% merge found components
CNM.displayMerging();

%% repeat processing and extract dF/F values
CNM.updateSpatial();
CNM.updateTemporal();
CNM.extractDFF();           % extract DF/F values.

%% plot components in GUI
% TODO disable this and all other interactive parts w/ flag so first pass can be
% done offline (then script to just load -> view/edit in gui?)
CNM.plotComponentsGUI();     % display all components

%% 

% From here to the end is "refining" the traces, including calculating dFF
% within blocks.
% TODO TODO TODO factor that refinement stuff into another fn (that needs timing
% info, etc, as additional arguments)

%{ refactor + uncomment
%sC_df = smoothdata(CNM.C_df, 2, 'loess', 25);

% TODO what is the meaning of these hardcoded cids?
% "for plotting"
% TODO random sample instead
%{
cids = 51:80;
figure, iosr.figures.multiwaveplot(sC_df(cids,:),...
    'gain', 5,...
    'reverseY', true);
%}

mat_filepath = fullfile(output_dir, [thorimage_id '_cnmf.mat']);
load(mat_filepath, 'ti');

disp(['Loading timing information from ' mat_filepath]);
% TODO sanity check these somehow
n_blocks = ti.num_trials;
fprintf('Reshaping CNMF output into %d blocks.\n', n_blocks);
cnm.T = CNM.T / n_blocks;
% TODO doesn't reshaping like this assumpe same number of frames per block?
% and that's not true, right?
% TODO what's empty list here?
cnm.C = reshape(full(CNM.C), size(CNM.C,1), [], n_blocks);
cnm.f = reshape(full(CNM.f), size(CNM.f,1), [], n_blocks);
cnm.R = reshape(full(CNM.R), size(CNM.R,1), [], n_blocks);

cnm.C_df = [];
cnm.F0 = [];
cnm.sC_df = [];
for i = 1:n_blocks
    % TODO TODO why is the block number already meaningful on the third
    % index?
    C = cnm.C(:,:,i);
    f = cnm.f(:,:,i);
    R = cnm.R(:,:,i);
    % TODO TODO will it turn out the CNM copy CNMc, that Remy was using
    % here earlier, was necessary? what was phenotype of error?
    [C_df, F0] = detrend_df_f(CNM.A, CNM.b,...
        C, f, R, CNM.options);
    sC_df = smoothdata(C_df, 2, 'loess', 20);
    cnm.C_df = [cnm.C_df C_df];
    cnm.sC_df = [cnm.sC_df sC_df];
    cnm.F0 = [cnm.F0 F0];
end

S.F0 = cnm.F0;
S.DFF = cnm.C_df;
S.sDFF = cnm.sC_df;
K = size(S.DFF,1);
S.K = K;
%%
%{
responseOptions.std_thr = 3;
responseOptions.thr = 1;

% TODO TODO fix... (just factor out?)
response = compute_response_vectors(S.sDFF, ti, responseOptions);

str = sprintf('%s: %s', date, thorimage_id);
ft = plot_S_traces(S, ti, flag(1:50), response, str);
%}

%% Save output
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% TODO just save 'ti' to different file and delete this hack. it is not clear
% that _cnmf could *only* hold ti (which is completely unrelated to CNMF) if ti
% is computed first

% TODO if going to keep hack, might as well skip append check, because it
% needs to exist already to load timing info anyway...

if exist(mat_filepath, 'file')
    save(mat_filepath, 'CNM', '-append');
else
    save(mat_filepath, 'CNM');
end

sCNM = struct(CNM);

% TODO what was A_pre? useful?
%sCNM.A_pre = full(sCNM.A_pre);

save(mat_filepath, 'sCNM', '-append');

save(mat_filepath, 'S', '-append');

end
