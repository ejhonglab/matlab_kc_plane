clear, clc

close all
date = '2018-10-21';

folder = sprintf(['/media/remy/remy-storage/Remy''s Dropbox Folder/',...
   'HongLab @ Caltech Dropbox/Remy/2018 data analysis/%s_analysis'], date);
%%
nams = {'003', '004', '005'};
Tc = table;
%%

for i = 1:numel(nams)
    nam = nams{i};    
    filename = [nam '_updated_cnmf3d.mat'];
    mat_filepath = fullfile(folder, 'cnmf3d', filename);
    data = load(mat_filepath);
    
    siz = size(data.corresp.binC);
    uniquely_assigned = data.corresp.uniquely_assigned;
    labels = uniquely_assigned(:,1);
    comp = uniquely_assigned(:,2);

    components = nan(siz(1),1);
    components(labels) = comp;

    var_nam = ['x' nam];
    Tc.(var_nam) = components;
    disp(nam)
end



%%

listing = dir(fullfile(folder, 'correspondence', '*.mat'));
disp({listing.name}')
%%
fn = load(fullfile(listing(2).folder, listing(2).name));
z2z = load(fullfile(listing(9).folder, listing(9).name));
%%
% M1: [components x labels]
% corresp.binC is [labels (X) x components(Y)]

M1 = zeros(size(fn.corresp.binC))';
for i = 1:numel(fn.corresp.one_to_one)
   M1(fn.corresp.uniquely_assigned(i,2),fn.corresp.uniquely_assigned(i,1))=1;
end
%%
M2 = zeros(size(z2z.corresp.binC));
for i = 1:numel(z2z.corresp.one_to_one)
    M2(z2z.corresp.uniquely_assigned(i,1), z2z.corresp.uniquely_assigned(i,2))=1;
end
%%
M = M1 * M2;
corresp = get_correspondence(M');
%%
nam = '002';
siz = size(corresp.binC);
uniquely_assigned = corresp.uniquely_assigned;
labels = uniquely_assigned(:,1);
comp = uniquely_assigned(:,2);

components = nan(siz(1),1);
components(labels) = comp;

var_nam = ['x' nam];
Tc.(var_nam) = components;
disp(nam)
%%
mat_filepath = fullfile(folder, 'correspondence', 'correspondence_zstackC.mat');
save(mat_filepath, 'Tc');
