clearvars

dataPath = 'Q:\Taeho\ZZZZ_FreeMax\MRI4D_04302023\Tumor Image';
path_MatData = fullfile(fileparts(dataPath), 'MatData');
ffn = fullfile(path_MatData, 'MTable.csv');
TT = readtable(ffn);

[val, ind] = unique(TT.iGroup);

NN = numel(val);
nR = ceil(NN/2);
nC = ceil(NN/nR);
hF1 = figure(1); clf
for n = 1:NN
    T = TT(TT.iGroup == val(n), :);
    for m = 1:height(T)
        idx = strfind(T.Folder{m}, '%');
        perc(m) = str2num(T.Folder{m}(idx-2:idx-1));
%         TTS{m} = T.Folder{m}(7:end-17);
    end

    hA1 = subplot(nR, nC, n, 'parent', hF1); hold on
    yyaxis(hA1, "left")
    plot(perc, T.avgY, 'o', 'MarkerSize', 6, 'LineWidth', 2, 'parent', hA1);
%     plot(perc, T.corY, 'ro', 'MarkerSize', 4, 'LineWidth', 2, 'parent', hA1);
%     plot(perc, T.segY, 'go', 'MarkerSize', 4, 'LineWidth', 2, 'parent', hA1);

    xlabel('%')
    ylabel('avgY')

    yyaxis(hA1, "right")
    plot(perc, T.avgR, 'o',  'MarkerSize', 6, 'LineWidth', 2,  'parent', hA1);
    ylabel('avgR')

    hA1.Title.String = T.Folder{m}(7:end-17);
    hA1.Title.Interpreter = 'none';
end