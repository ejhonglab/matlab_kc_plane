close all, clear, clc

date = '2019-01-18';

folder = sprintf(['/media/remy/remy-storage/Remy''s Dropbox Folder/',...
   'HongLab @ Caltech Dropbox/Remy/2019 data analysis/%s_analysis'], date);
% folder = sprintf('D:\Remy\2018 imaging data\%', date);

nam = '_007';
mat_filepath = fullfile(folder, 'cnmf', [nam '_cnmf.mat']);
disp(matwho(mat_filepath))

%% 
load(mat_filepath);

%% edit traces
CNMc = copy(CNM);
CNMc.options.df_window=ti.fpt;
CNMc.extractDFF();

cnm.T = ti.fpt;
cnm.C = reshape(full(CNMc.C), size(CNMc.C,1), [], ti.num_trials);
cnm.f = reshape(full(CNMc.f), size(CNMc.f,1), [], ti.num_trials);
cnm.R = reshape(full(CNMc.R), size(CNMc.R,1), [], ti.num_trials);

cnm.C_df = [];
cnm.F0 = [];
cnm.sC_df = [];
for i = 1:ti.num_trials
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
%% compute response vectors
responseOptions.std_thr = 3;
responseOptions.thr = 1;
response = compute_response_vectors(S.sDFF, ti, responseOptions);

flag = sum(response.bin_response,2);
flag = find(flag>1);
%%
str = sprintf('%s: %s', date, nam);
ft = plot_S_traces(S, ti, flag, response,str);