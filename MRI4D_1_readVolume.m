clearvars

dataPath = 'Q:\Taeho\ZZZZ_FreeMax\MRI4D_04302023\Tumor Image';
path_MatData = fullfile(fileparts(dataPath), 'MatData');

d= dir(dataPath);
dfolders = d([d(:).isdir]);
dfolders = dfolders(~ismember({dfolders(:).name},{'.','..'}));

FolderList = {dfolders.name}';

NN = numel(FolderList);
for n = 1:NN
    ffn_mat = fullfile(path_MatData, [FolderList{n}, '.mat']);
    if exist(ffn_mat, 'file')
        disp([num2str(n), '/', num2str(NN), ' already processed...', FolderList{n}]);
    else
        VolumePath = fullfile(dataPath, FolderList{n});
        disp([num2str(n), '/', num2str(NN), '...', FolderList{n}]);
        [V, SI] = dicomreadVolume(VolumePath);
        save(ffn_mat, 'V', 'SI');
    end
end