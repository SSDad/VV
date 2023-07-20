function Callback_Pushbutton_LoadImagePanel_Seg(src, evnt)

hFig = ancestor(src, 'Figure');
data = guidata(hFig);

%% segmentation
bPlot = 1;
V = data.Image.V;
SI = data.Image.SI;

[AXL, COR, SEG] = fun_Sphere(V, SI, bPlot);

%% Axial
% update view
hPlotObj = data.Panel.View.Comp.hPlotObj;
hPlotObj(1).Image.CData = AXL.I; 
set(hPlotObj(1).bd1, 'XData', AXL.bdxw, 'YData', AXL.bdyw);
set(hPlotObj(1).bd2, 'XData', AXL.bdx2w, 'YData', AXL.bdy2w);
set(hPlotObj(1).cent2, 'XData', AXL.xcent2w, 'YData', AXL.ycent2w);

% update crosshair
[~, y] = intrinsicToWorld(data.Image.RC,  1, AXL.iSliceA); 
hPlotObj(3).X.Position(2) = y;
hPlotObj(4).X.Position(2) = y;

%% coronal and sagittal planes
% update view
hPlotObj(3).Image.CData = COR.IC; 
hPlotObj(4).Image.CData = SEG.IS;

% update Xhair
hPlotObj(1).X.Position = [AXL.xcent2w AXL.ycent2w]; 
hPlotObj(3).X.Position(1) = AXL.xcent2w; 
hPlotObj(4).X.Position(1) = AXL.ycent2w; 

%% cor gourd and ball
set(hPlotObj(3).bd1, 'XData', COR.bdxwC, 'YData', COR.bdywC);
set(hPlotObj(3).cent, 'XData', COR.xc2, 'YData', COR.yc2);
set(hPlotObj(3).ball, 'XData', COR.xxc, 'YData', COR.yyc);

%% seg gourd and ball
set(hPlotObj(4).bd1, 'XData', SEG.bdxwS, 'YData', SEG.bdywS);
set(hPlotObj(4).cent, 'XData', SEG.xc2, 'YData', SEG.yc2);
set(hPlotObj(4).ball, 'XData', SEG.xxc, 'YData', SEG.yyc);

guidata(hFig, data);


function   [xc,yc,R,a] = circfit(x,y)
%
%   [xc yx R] = circfit(x,y)
%
%   fits a circle  in x,y plane in a more accurate
%   (less prone to ill condition )
%  procedure than circfit2 but using more memory
%  x,y are column vector where (x(i),y(i)) is a measured point
%
%  result is center point (yc,xc) and radius R
%  an optional output is the vector of coeficient a
% describing the circle's equation
%
%   x^2+y^2+a(1)*x+a(2)*y+a(3)=0
%
%  By:  Izhak bucher 25/oct /1991, 
    x=x(:); y=y(:);
   a=[x y ones(size(x))]\[-(x.^2+y.^2)];
   xc = -.5*a(1);
   yc = -.5*a(2);
   R  =  sqrt((a(1)^2+a(2)^2)/4-a(3));

end

end