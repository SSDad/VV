function Callback_Pushbutton_LoadImagePanel_LoadImage(src, evnt)

hFig = ancestor(src, 'Figure');
data = guidata(hFig);

%% load image info
td = tempdir;
path_info = fullfile(td, 'VV');
ffn_info = fullfile(path_info, 'info.mat');
if ~exist(path_info, 'dir')
    mkdir(path_info);
end
if ~exist(ffn_info, 'file')
    VolumePath = uigetdir();
    [dataPath] = fileparts(VolumePath);
    save(ffn_info, 'dataPath');
else
    load(ffn_info);
    VolumePath = uigetdir(dataPath);
end

data.FileInfo.VolumePath = VolumePath;

%% load image data
hWB = waitbar(0, 'Loading Volume...');

path_MatData = fullfile(dataPath, 'MatData');
[~, fd_VolumePath] = fileparts(VolumePath);
ffn_mat = fullfile(path_MatData, [fd_VolumePath, '.mat']);

if exist(ffn_mat, 'file')
    load(ffn_mat);
else
    [V, SI] = dicomreadVolume(VolumePath);
    save(ffn_mat, 'V', 'SI');
end

data.V = V;
data.SI = SI;

waitbar(1/3, hWB, 'Initializing View...');

% %%  first slice in image viewer
% % first slice
% iSlice = 1;
% I = data.Image.Images(:,:,iSlice);
%     
% % fill Image Info 
% data.Panel.LoadImage.Comp.hEdit.ImageInfo(2).String = [num2str(mImgSize), 'x', num2str(nImgSize)];
% data.Panel.LoadImage.Comp.hEdit.ImageInfo(2).ForegroundColor = 'c';
% 
% data.Panel.LoadImage.Comp.hEdit.ImageInfo(3).String = num2str(Image.RA.PixelExtentInWorldX);
% data.Panel.LoadImage.Comp.hEdit.ImageInfo(3).ForegroundColor = 'c';
% 
% %  image viewer
% hA = data.Panel.View.Comp.hAxis.Image;
% hPlotObj.Image = imshow(I, Image.RA, [], 'parent', hA);
% axis(data.Panel.View.Comp.hAxis.Image, 'tight', 'equal')
% 
% % template box
% hPlotObj.LBox = rectangle(hA, 'Position', [0 0 0 0], 'EdgeColor', 'g', 'LineWidth', 1);
% hPlotObj.TBox = rectangle(hA, 'Position', [0 0 0 0], 'EdgeColor', 'c', 'LineWidth', 1);
% hPlotObj.CBox = rectangle(hA, 'Position', [0 0 0 0], 'EdgeColor', 'b', 'LineWidth', 1);
% 
% data.Panel.View.Comp.hPlotObj = hPlotObj;
% 
% %% Slice Slider
% hSS = data.Panel.SliceSlider.Comp.hSlider.Slice;
% hSS.Min = 1;
% hSS.Max = nSlices;
% hSS.Value = iSlice;
% hSS.SliderStep = [1 10]/(nSlices-1);
% 
% data.Panel.SliceSlider.Comp.hText.nImages.String = [num2str(iSlice), ' / ', num2str(nSlices)];

waitbar(1, hWB, 'Volume loaded!');
pause(1);
close(hWB);

% %% contrast
% yc = histcounts(I, max(I(:))+1);
% yc = log10(yc);
% yc = yc/max(yc);
% xc = 1:length(yc);
% xc = xc/max(xc);
% 
% data.Panel.ContrastBar.Comp.hPlotObj.Hist.XData = xc;
% data.Panel.ContrastBar.Comp.hPlotObj.Hist.YData = yc;
% 
% data.Snake.SlitherDone = false;

guidata(hFig, data);