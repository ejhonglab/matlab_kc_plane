
% TODO TODO rewrite this to just call hdf5_to_mat
clear, clc
%folder = 'D:\MB_team\Natural_odors\2019-01-24\3';

folder = '/home/tom/HongLab @ Caltech Dropbox/Remy/2018 data/2018_02_02';
%folder = "/media/remy/remy-storage/Remy's Dropbox Folder/HongLab @ Caltech Dropbox/Remy/2018 data/2018_02_02";
%folder = '/run/user/1000/gvfs/smb-share:server=10.0.0.3,share=main/mb_team/2019-01-23/2';

%listing = dir(fullfile(folder, '**\Episode001.h5'));
listing = dir(fullfile(folder, '*/Episode001.h5'));

disp(listing)

% TODO what is this doing?
disp({listing.folder}')

% TODO fix loop
for i = 1
    % TODO check if exists first?
    filepath = fullfile(listing(i).folder, listing(i).name);
    S = FnLoadSyncEpisode(filepath);
    
    folders = regexp(filepath, filesep, 'split');
    disp(folders)
    thorSyncName = folders(end-1);
    disp(thorSyncName)
    fname = [thorSyncName{1} '.mat'];
    save(fullfile(outfolder, fname), '-struct', 'S');
    fprintf('%s saved\n', fullfile(outfolder, fname));
end


disp('All ThorSync HDF5 files saved as .mat files.')

