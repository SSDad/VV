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

waitbar(1/2, hWB, 'Initializing View...');

%%  middle slice in image viewer
iSliceA = round(SI.ImageSize(3)/2);
IA = V(:,:,iSliceA);
xmin = SI.PatientPositions(iSliceA, 1);
ymin = SI.PatientPositions(iSliceA, 2);
xL = SI.PixelSpacings(1)*(SI.ImageSize(2)-1);
yL = SI.PixelSpacings(2)*(SI.ImageSize(1)-1);
RA = imref2d(SI.ImageSize(1:2), [xmin xmin+xL], [ymin ymin+yL]);
hA = data.Panel.View.Comp.hAxis.Image(1);
hPlotObj.Image(1) = imshow(IA, RA, [], 'parent', hA);
axis(hA, 'tight', 'equal')

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