clear, clc
%folder = 'D:\Remy\2018 imaging data\2018-07-05';
%folder = 'D:\Remy\2018 pid data\2018-07-12';
%folder = 'D:\Remy\2018 imaging data\2018-10-26';
%folder = 'D:\Remy\2018 pid data\2018-09-26_pid';
%folder = 'D:\Kellan\2018-10-29';
%folder = 'D:\Remy\2019 data\2019-01-16';

folder = 'D:\MB_team\Natural_odors\2019-01-24\3';
%folder = 'D:\MB_team\Odor_space\2019-01-24\1';
listing = dir(fullfile(folder, '**\Episode001.h5'));
disp(listing)

outfolder = fullfile(folder, 'thorSync');
if ~exist(outfolder, 'dir')
    mkdir(outfolder)
end

disp({listing.folder}')

%%
for i = 1
   filepath = fullfile(listing(i).folder, listing(i).name);
   S = FnLoadSyncEpisode(filepath);
   
   folders = regexp(filepath, filesep, 'split');
   thorSyncName = folders(end-1);
   fname = [thorSyncName{1} '.mat'];
   fprintf('%s saved\n', fullfile(outfolder, fname));
   save(fullfile(outfolder, fname), '-struct', 'S');
end


disp('All thorSync files saved as .mat files.')
%%
% figure
% yyaxis left
% plot(time, Piezo_Monitor)
% %yyaxis right
% %plot(time, Frame_Out)
% %%

% end