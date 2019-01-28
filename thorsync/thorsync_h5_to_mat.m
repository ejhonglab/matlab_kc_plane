
function thorsync_h5_to_mat(thorsync_dir, output_dir)

hdf5_path = fullfile(thorsync_dir, 'Episode001.h5');
if exist(hdf5_path, 'file')
    return
end

fprintf('converting ThorSync HDF5 in %s to .mat in %s...', thorsync_dir, ...
    output_dir)

outfolder = fullfile(output_dir, 'thorsync');
if ~exist(outfolder, 'dir')
    mkdir(outfolder);
end

S = FnLoadSyncEpisode(hdf5_path);

% TODO redundant? simpler way?
folders = regexp(hdf5_path, filesep, 'split');
thorSyncName = folders(end-1);
fname = [thorSyncName{1} '.mat'];
save(fullfile(outfolder, fname), '-struct', 'S');
fprintf('%s saved\n', fullfile(outfolder, fname));

end
