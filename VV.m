function VV

%% main window
hFig = figure('MenuBar',            'none', ...
                    'Toolbar',              'none', ...
                    'HandleVisibility',  'callback', ...
                    'Name',                'Volume Viewer - Department of Radiation Oncology, Washington University in St. Louis', ...
                    'NumberTitle',      'off', ...
                    'Units',                 'normalized',...
                    'Position',             [0.1 0.1 0.6 0.8],...
                    'Color',                 'black', ...
                    'CloseRequestFcn', @figCloseReq, ...
                    'Visible',               'on');

addToolbar(hFig);
                
data.Panel = addPanel(hFig);
data.Panel.LoadImage.Comp = addComponents2Panel_LoadImage(data.Panel.LoadImage.hPanel);
% data.Panel.Selection.Comp = addComponents2Panel_Selection(data.Panel.Selection.hPanel);
% data.Panel.View.Comp = addComponents2Panel_View(data.Panel.View.hPanel);
% data.Panel.Snake.Comp = addComponents2Panel_Snake(data.Panel.Snake.hPanel);
% data.Panel.Body.Comp = addComponents2Panel_Body(data.Panel.Body.hPanel);
% data.Panel.ContrastBar.Comp = addComponents2Panel_ContrastBar(data.Panel.ContrastBar.hPanel);
% data.Panel.SliceSlider.Comp = addComponents2Panel_SliceSlider(data.Panel.SliceSlider.hPanel);
% 
% data.Panel.Point.Comp = addComponents2Panel_Point(data.Panel.Point.hPanel);
% 
% data.Panel.About.Comp = addComponents2Panel_About(data.Panel.About.hPanel);

guidata(hFig, data);