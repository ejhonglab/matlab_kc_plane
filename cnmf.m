
% TODO any conflicts w/ installed CaImAn? rename?
function cnmf(tiff_path, output_dir)
% TODO TODO just take path to tif and strip to get thorimage id
% TODO maybe take some of cnmf params too?

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

%%
if ndims(Y) == 4
    [d1,d2,d3,T] = size(Y);                            % dimensions of dataset
else
    [d1,d2,T] = size(Y);
    d3 = 1;
end
d = d1*d2*d3;                                          % total number of pixels

%%
CNM = CNMF;
K = 600;        % # of expected components
tau = 3;        % 

% Order of autoregressive system
% (p=0 no dynamics, p=1 just decay, p=2, both rise and decay)

% TODO these parameters depend on frame rate / image size / noise / etc
% one setting doesnt tend to work globally
% most sensitive: tau, merge_thr, energy_thresh, thr_method

p = 0;                                            
merge_thr = 0.95;                                  % merging threshold

options = CNMFSetParms(...   
    'd1',d1,'d2',d2,...                         % dimensions of datasets
    'fr', 4.55,...
    'p',p,...                                   % order of AR dynamics    
    'gSig',tau,...                                % half size of neuron
    'merge_thr',merge_thr,...                        % merging threshold  
    'nb', 2,...                                  % number of background components    
    'min_SNR',3,...                             % minimum SNR threshold
    'space_thresh',0.3,...                      % space correlation threshold
    'time_thresh',0.3,...                      % space correlation threshold
    'cnn_thr',0.2,...                            % threshold for CNN classifier    
    'init_method', 'greedy',...
    'search_method', 'ellipse',...
    'thr_method', 'nrg',...
    'nrgthr', .95,...
    'noise_norm', true,...
    'bas_nonneg', false,...
    'cont_threshold', .95);
%'cluster_pixels', true,...   
%'dist', 3,...
%'max_size', 3,...
%'min_size', 1,...

%%
% Below is the standard processing pipeline. This processing can be
% executed in one shot using the CNM.fit function:
 %CNM.fit(Y,options,K)
%% load the dataset and create the object
CNM.readFile(tiff_path);                         % insert path to file here  
CNM.optionsSet(options);                        % setup the options structure
CNM.gSig = options.gSig;
%% Process the dataset

CNM.preprocess;             % preprocessing (compute some quantities)
CNM.initComponents(K);      % initialization 
CNM.plotCenters()           % plot center of ROIs detected during initialization
%% POSTPROCESSING BEGINS

%% Update components
CNM.updateSpatial();        % update spatial components
CNM.updateTemporal(0);      % update temporal components (do not deconvolve at this point)

%% Plot contours
CNM.cm = [];
CNM.contours=[];
CNM.COM();
CNM.K = size(CNM.A,2);
CNM.CI=[];
%CNM.correlationImage();
figure,CNM.plotContours();
%% If you want to edit components, run this cell

% TODO flag to skip, among w/ other interactive stuff
CNM.manuallyRefineComponents();
CNM.updateSpatial();        % update spatial components
CNM.updateTemporal(0);      % update temporal components (do not deconvolve at this point)


%% component classification
CNM.evaluateComponents();   % evaluate spatial components based on their correlation with the data
CNM.eventExceptionality();  % evaluate traces

ind_throw = find(~CNM.keep_eval);

ind_keep = find(CNM.keep_eval);

figure
CNM.plotContours(1,[],[],ind_keep); % plot components to get rid of
%%
%get rid of components if sure
CNM.keepComponents(CNM.keep_eval);      % keep the components that are above certain thresholds
%% merge components
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
sC_df = smoothdata(CNM.C_df, 2, 'loess', 25);
% TODO what is the meaning of these hardcoded cids?
% "for plotting"
cids = 51:80;

%cids([1 24])=[];

figure, iosr.figures.multiwaveplot(sC_df(cids,:),...
    'gain', 5,...
    'reverseY', true);

n_blocks = ti.num_trials;
cnm.T = CNM.T / n_blocks;
cnm.C = reshape(full(CNM.C), size(CNM.C,1), [], n_blocks);
cnm.f = reshape(full(CNM.f), size(CNM.f,1), [], n_blocks);
cnm.R = reshape(full(CNM.R), size(CNM.R,1), [], n_blocks);

cnm.C_df = [];
cnm.F0 = [];
cnm.sC_df = [];
for i = 1:n_blocks
    C = cnm.C(:,:,i);
    f = cnm.f(:,:,i);
    R = cnm.R(:,:,i);
    [C_df,F0] = detrend_df_f(CNMc.A,CNMc.b,...
        C,f,R,CNMc.options);
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
responseOptions.std_thr = 3;
responseOptions.thr = 1;

% TODO TODO fix... (just factor out?)
response = compute_response_vectors(S.sDFF, ti, responseOptions);

str = sprintf('%s: %s', date, thorimage_id);
ft = plot_S_traces(S, ti, flag(1:50), response, str);
%}

output_dir = fullfile(output_dir, 'cnmf');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% save CNM
mat_filepath = fullfile(output_dir, [thorimage_id '_cnmf.mat']);

% TODO just save 'ti' to different file and delete this hack. it is not clear
% that _cnmf could *only* hold ti (which is completely unrelated to CNMF) if ti
% is computed first
if exist(mat_filepath, 'file')
    save(mat_filepath, 'CNM', '-append');
else
    save(mat_filepath, 'CNM');
end

% TODO if i'm just going to use the MATLAB engine to read the MAT file
% anyway, maybe get rid of the sCNM output?
% TODO TODO test that output generated this way can actually be read by my
% python stuff...
sCNM = struct(CNM);
% i forget, was converting from the sparse arrays actually necessary in the end
% to load successfully?
sCNM.A = full(sCNM.A);
sCNM.A_pre = full(sCNM.A_pre);

save(mat_filepath, 'sCNM', '-append');

end
