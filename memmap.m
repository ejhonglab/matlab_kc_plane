clear, clc
date = '2018-10-21';

folder = sprintf(['/media/remy/remy-storage/Remy''s Dropbox Folder/',...
   'HongLab @ Caltech Dropbox/Remy/2018 data analysis/%s_analysis'], date);

nams = {'002', '003', '004', '005'};
%%
for i = numel(nams)
    nam = nams{i};
    filename = [nam '_rig.tif'];
    filepath = fullfile(folder, 'tif stacks', filename);

    Y = single(loadtiff(filepath));
    Y = reshape(Y, 256, 256, 11, []);
    Y = Y(:,:,2:10,:);

    sizY = size(Y);
    Yr = reshape(Y,prod(sizY(1:end-1)),[]);
    nY = min(Yr(:));
    %Yr = Yr - nY;
    %save([filename(1:end-3),'mat'],'Yr','Y','nY','sizY','-v7.3');
    savefast([filepath(1:end-3),'mat'],'Yr','Y','nY','sizY');
end

%%
file=[filepath(1:end-3),'mat'];
data = matfile(file,'Writable',true);
%[d1,d2,d3,T] = size(data,'Y');
[d1,d2,T] = size(data,'Y');
is_memmaped = 1;

