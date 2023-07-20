function Callback_Pushbutton_LoadImagePanel_Seg(src, evnt)

hFig = ancestor(src, 'Figure');
data = guidata(hFig);

%% segmentation
bPlot = 0;

hPlotObj = data.Panel.View.Comp.hPlotObj;

ImageSize = data.Image.SI.ImageSize;
iSliceA = round(ImageSize(3)/2);  % middle slice
V = squeeze(data.Image.V);
I = V(:, :, iSliceA);
hPlotObj(1).Image.CData = I; % update view

[~, y] = intrinsicToWorld(data.Image.RC,  1, iSliceA); % update crosshair
hPlotObj(3).X.Position(2) = y;
hPlotObj(4).X.Position(2) = y;

J = rescale(I);
%T = graythresh(J);
BW = imbinarize(J);
BW2 = imfill(BW, 'holes');
BW3 = bwareaopen(BW2, round(ImageSize(1)*ImageSize(2)/8));

if bPlot
    hF = figure(11);
    subplot(221),     imshow(J, []);
    subplot(222),    imshow(BW);
    subplot(223),    imshow(BW2);
    subplot(224),    imshow(BW3);
end

B = bwboundaries(BW3);

    boundary = B{1};
    polyin = polyshape(boundary(:, 2), boundary(:, 1));
    [xlim, ylim] = boundingbox(polyin);
    bdx = boundary(:,2);
    bdy = boundary(:,1);

    cc(1) = mean(xlim);
    cc(2) = mean(ylim);
    bX = diff(xlim); % box X
    D1 = bX*(60/196); % diameter 1
    A = 3;
    x0 = cc(1) - D1/2 - A;
    y0 = cc(2) - D1/2 - A;
    
    rect = [x0 y0 D1+A*2 D1+A*2];
    rect = round(rect);

    if bPlot
        subplot(221), hold on
        plot(bdx, bdy, 'g', 'LineWidth', 2);
        plot(cc(2), cc(1), 'x', 'LineWidth', 2);
        rectangle('Position', rect, 'EdgeColor', 'r')
    end

    [bdxw, bdyw] = intrinsicToWorld(data.Image.RA,  bdx, bdy); 
    set(hPlotObj(1).bd1, 'XData', bdxw, 'YData', bdyw);

%% cylinder crop
JC = imcrop(J, rect);
JC = rescale(JC);
JCBW = imbinarize(JC);
JCBW2 = bwareaopen(~JCBW, round(prod(size(JC))/6));
JCB = bwboundaries(JCBW2);
JCboundary1 = JCB{1};

[xc,yc,R,a] = circfit(JCboundary1(:,2), JCboundary1(:,1));  % find center and radii

if bPlot
    hF = figure(12); clf
    subplot(221), imshow(JC, []);
    hold on, plot(JCboundary1(:,2), JCboundary1(:,1), 'g', 'LineWidth', 2);
    viscircles([xc yc], R)

    subplot(222), imshow(JCBW, []);
    subplot(223), imshow(~JCBW, []);
    subplot(224), imshow(JCBW2, []);
end

% circle boundary
bdx2 = JCboundary1(:,2)+x0-1;
bdy2 = JCboundary1(:,1)+y0-1;
[bdx2w, bdy2w] = intrinsicToWorld(data.Image.RA,  bdx2, bdy2);
set(hPlotObj(1).bd2, 'XData', bdx2w, 'YData', bdy2w);

% circle center
xcent2 = xc+x0-1;
ycent2 = yc+y0-1;
[xcent2w, ycent2w] = intrinsicToWorld(data.Image.RA,  xcent2, ycent2);
set(hPlotObj(1).cent2, 'XData', xcent2w, 'YData', ycent2w);


% find coronal and sagittal planes
iSliceC = round(ycent2);
IC = squeeze(V(iSliceC, :, :));
IC = IC';
hPlotObj(3).Image.CData = IC; % update view

iSliceS = round(xcent2);
IS = squeeze(V(:, iSliceS, :));
IS = IS';
hPlotObj(4).Image.CData = IS; % update view

% update Xhair
hPlotObj(1).X.Position = [xcent2w ycent2w]; 
hPlotObj(3).X.Position(1) = xcent2w; 
hPlotObj(4).X.Position(1) = ycent2w; 

%% gourd crop
R2 = round(R*(20/29));

% coronal
KC = IC(:, iSliceC-R2:iSliceC+R2);
KC = rescale(KC);
KCBW = imbinarize(1-KC);
KCBW2 = bwareaopen(KCBW, round(prod(size(KC))/16));
KCBW3 = imfill(KCBW2, 'holes');

BB = bwboundaries(KCBW3);
bd = BB{1};
% cc = mean(bd);

[bdxwC, bdywC] = intrinsicToWorld(data.Image.RC,  bd(:, 2)+iSliceC-R2-1, bd(:, 1)); 
set(hPlotObj(3).bd1, 'XData', bdxwC, 'YData', bdywC);

% ball
polyin = polyshape(bdxwC, bdywC);
[xlim, ylim] = boundingbox(polyin);
R2 = diff(xlim)/2;
ind = find(bdywC < ylim(1)+R2*3/2);
bdx2 = bdxwC(ind);
bdy2 = bdywC(ind);
[xc2, yc2, R2, a] = circfit(bdx2, bdy2);

set(hPlotObj(3).cent, 'XData', mean(xc2), 'YData', mean(yc2));

theta = linspace(0, pi*2, 100);
xx = xc2+R2*cos(theta);
yy = yc2+R2*sin(theta);

set(hPlotObj(3).ball, 'XData', xx, 'YData', yy);

if bPlot
    figure(13)
    subplot(221), imshow(KC, []); axis xy on;
    subplot(222), imshow(KCBW, []); axis xy on;
    subplot(223), imshow(KCBW2, []); axis xy on;
    subplot(224), imshow(KCBW3, []); axis xy on;
    
    subplot(221), hold on
    plot(bd(:,2), bd(:,1), 'g', 'LineWidth', 2);
    plot(cc(2), cc(1), 'rx', 'LineWidth', 2);
end


% sagittal
KS = IS(:, iSliceS-R2:iSliceS+R2);
KS = rescale(KS);
KSBW = imbinarize(1-KS);
KSBW2 = bwareaopen(KSBW, round(prod(size(KS))/8));
KSBW3 = imfill(KSBW2, 'holes');

BB = bwboundaries(KSBW3);
bd = BB{1};
cc = mean(bd);

if bPlot
    figure(14)
    subplot(221), imshow(KS, []); axis xy on;
    subplot(222), imshow(KSBW, []); axis xy on;
    subplot(223), imshow(KSBW2, []); axis xy on;
    subplot(224), imshow(KSBW3, []); axis xy on;
    
    subplot(221), hold on
    plot(bd(:,2), bd(:,1), 'g', 'LineWidth', 2);
    plot(cc(2), cc(1), 'rx', 'LineWidth', 2);
end

[bdxwS, bdywS] = intrinsicToWorld(data.Image.RS,  bd(:, 2)+iSliceS-R2-1, bd(:, 1)); 
set(hPlotObj(4).bd1, 'XData', bdxwS, 'YData', bdywS);

% ball
polyin = polyshape(bdxwS, bdywS);
[xlim, ylim] = boundingbox(polyin);
R2 = diff(xlim)/2;
ind = find(bdywS < ylim(1)+R2*3/2);
bdx2 = bdxwS(ind);
bdy2 = bdywS(ind);
[xc2, yc2, R2, a] = circfit(bdx2, bdy2);

set(hPlotObj(4).cent, 'XData', mean(xc2), 'YData', mean(yc2));

theta = linspace(0, pi*2, 100);
xx = xc2+R2*cos(theta);
yy = yc2+R2*sin(theta);

set(hPlotObj(4).ball, 'XData', xx, 'YData', yy);

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