clear;
gcp;

% same demo as demo_script.m but using the class @CNMF
%% load file

date = '2019-01-18';

folder = sprintf(['/media/remy/remy-storage/Remy''s Dropbox Folder/',...
   'HongLab @ Caltech Dropbox/Remy/2019 data analysis/%s_analysis'], date);
% folder = sprintf('D:\Remy\2018 imaging data\%', date);

nam = '_006';
filename = [nam '_nr.tif'];
filepath = fullfile(folder, 'tif stacks', filename);
Y = loadtiff(filepath);

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
K = 700;        % # of expected components
tau = 2;        % 

p = 2;                                            % order of autoregressive system (p = 0 no dynamics, p=1 just decay, p = 2, both rise and decay)
merge_thr = 0.95;                                  % merging threshold

options = CNMFSetParms(...   
    'd1',d1,'d2',d2,...                         % dimensions of datasets
    'fr', 4.55,...
    'p',p,...                                   % order of AR dynamics    
    'gSig',tau,...                                % half size of neuron
    'merge_thr',merge_thr,...                        % merging threshold  
    'nb', 2,...                                  % number of background components    
    'min_SNR',3,...                             % minimum SNR threshold
    'space_thresh',0.5,...                      % space correlation threshold
    'cnn_thr',0.2,...                            % threshold for CNN classifier    
    'init_method', 'greedy',...
    'search_method', 'ellipse',...
    'thr_method', 'nrg',...
    'nrgthr', .95,...
    'noise_norm', true);
%'cluster_pixels', true,...   
%'dist', 3,...
%'max_size', 3,...
%'min_size', 1,...

%%
% Below is the standard processing pipeline. This processing can be
% executed in one shot using the CNM.fit function:
 CNM.fit(Y,options,K)
%% load the dataset and create the object
CNM.readFile(filepath);                         % insert path to file here  
CNM.optionsSet(options);                        % setup the options structure
CNM.gSig = options.gSig;
%% Process the dataset
close all
CNM.preprocess;             % preprocessing (compute some quantities)
CNM.initComponents(K);      % initialization
CNM.plotCenters()           % plot center of ROIs detected during initialization
%% Update components
CNM.updateSpatial();        % update spatial components
CNM.updateTemporal(0);      % update temporal components (do not deconvolve at this point)
%% Plot contours
CNM.cm = [];
CNM.contours=[];
CNM.COM();
CNM.K = size(CNM.A,2);
CNM.CI=[];
figure,CNM.plotContours();
%% component classification
CNM.evaluateComponents();   % evaluate spatial components based on their correlation with the data
CNM.eventExceptionality();  % evaluate traces

ind_throw = find(~CNM.keep_eval);
ind_keep = find(CNM.keep_eval);

figure
CNM.plotContours(1,[],[],ind_keep); % plot components to get rid of

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
CNM.plotComponentsGUI();     % display all components