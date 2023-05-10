function hToolbar = addToolbar(hFig)

%% toolbar          
hToolbar = uitoolbar ('parent', hFig);

% standard
uitoolfactory(hToolbar,'Exploration.ZoomIn');
uitoolfactory(hToolbar,'Exploration.ZoomOut');
uitoolfactory(hToolbar,'Exploration.Pan');                
uitoolfactory(hToolbar,'Exploration.DataCursor');                
% uitoolfactory(hToolbar, 'Exploration.Rotate');