function [AXL, COR, SEG] = fun_Sphere(V, SI, bPlot)

    ImageSize = SI.ImageSize;

    %% RA, RS, RC
    xmin = SI.PatientPositions(1, 1);
    ymin = SI.PatientPositions(1, 2);
    xL = SI.PixelSpacings(1)*(ImageSize(2)-1);
    yL = SI.PixelSpacings(2)*(ImageSize(1)-1);
    RA = imref2d(SI.ImageSize(1:2), [xmin xmin+xL], [ymin ymin+yL]);
    
    ymin = SI.PatientPositions(1, 2);
    zmin = SI.PatientPositions(1, 3);
    zmax = SI.PatientPositions(end, 3);
    RS = imref2d(SI.ImageSize([3, 1]), [ymin ymin+yL], [zmin zmax]);
    
    RC = imref2d(SI.ImageSize([3, 2]), [xmin xmin+xL], [zmin zmax]);

    AXL.RA = RA;
    COR.RC = RC;
    SEG.RS = RS;

    %% axial 
    iSliceA = round(SI.ImageSize(3)/2);  % middle slice
    V = squeeze(V);
    I = V(:, :, iSliceA);
       
    [BW] = fun_findTank(I);

    B = bwboundaries(BW);
    boundary = B{1};
    bdx = boundary(:,2);
    bdy = boundary(:,1);
    
    polyin = polyshape(bdx, bdy);
    [xlim, ylim] = boundingbox(polyin);
    cc(1) = mean(xlim);
    cc(2) = mean(ylim);

    [bdxw, bdyw] = intrinsicToWorld(RA,  bdx, bdy); 

    AXL.I = I;
    AXL.iSliceA = iSliceA;
    AXL.bdxw = bdxw;
    AXL.bdyw = bdyw;

    if bPlot
        figure(11); hA1 = subplot(221);
        imshow(I, RA, [], 'parent', hA1); hold on

        plot(bdxw, bdyw, 'g', 'LineWidth', 2);
    end


    %% cylinder crop
    % cylinder box
    bX = diff(xlim); % box X
    D1 = bX*(60/196); % diameter 1
    A = 3;
    x0 = cc(1) - D1/2 - A;
    y0 = cc(2) - D1/2 - A;
    
    rect = [x0 y0 D1+A*2 D1+A*2];
    rect = round(rect);

    % cylinder crop
    J = rescale(I);
    JC = imcrop(J, rect);
    JC = rescale(JC);
    JCBW = imbinarize(JC);
    JCBW2 = bwareaopen(~JCBW, round(prod(size(JC))/6));
    
    JCB = bwboundaries(JCBW2);
    JCboundary1 = JCB{1};
    
    [xc,yc,R1,a] = circfit(JCboundary1(:,2), JCboundary1(:,1));  % cylinder center and radii

    % circle boundary
    bdx2 = JCboundary1(:,2)+x0-1;
    bdy2 = JCboundary1(:,1)+y0-1;
    [bdx2w, bdy2w] = intrinsicToWorld(RA,  bdx2, bdy2);
    
    % circle center
    xcent2 = xc+x0-1;
    ycent2 = yc+y0-1;
    [xcent2w, ycent2w] = intrinsicToWorld(RA,  xcent2, ycent2);

    AXL.bdx2w = bdx2w;
    AXL.bdy2w = bdy2w;
    AXL.xcent2w = xcent2w;
    AXL.ycent2w = ycent2w;

    if bPlot
        hF = figure(12); clf
        subplot(221), imshow(JC, []); hold on
        plot(JCboundary1(:,2), JCboundary1(:,1), 'g', 'LineWidth', 2);
        viscircles([xc yc], R1)
    
        subplot(222), imshow(JCBW, []);
        subplot(223), imshow(~JCBW, []);
        subplot(224), imshow(JCBW2, []);

        figure(11), subplot(221)
        plot(bdx2w, bdy2w, 'g', 'LineWidth', 2);
        plot(xcent2w, ycent2w, 'rx', 'LineWidth', 2);

    end

    %% find coronal and sagittal planes
    xcent2 = xc+x0-1;
    ycent2 = yc+y0-1;

    iSliceC = round(ycent2);
    IC = squeeze(V(iSliceC, :, :));
    IC = IC';

    iSliceS = round(xcent2);
    IS = squeeze(V(:, iSliceS, :));
    IS = IS';

    COR.IC = IC;
    SEG.IS = IS;

    if bPlot
        figure(11), hA3 = subplot(223);
        imshow(IC, RC, [], 'parent', hA3);
        hA3.YDir = 'normal';

        figure(11), hA4 = subplot(224);
        imshow(IS, RS, [], 'parent', hA4);
        hA4.YDir = 'normal';

        linkaxes([hA1 hA3], 'x');
    end

    %% coronal gourd
    R2 = round(R1*(20/29));

    KC = IC(:, iSliceC-R2:iSliceC+R2);
    KC = rescale(KC);
    KCBW = imbinarize(1-KC);
    KCBW2 = bwareaopen(KCBW, round(prod(size(KC))/16));
    KCBW3 = imfill(KCBW2, 'holes');
    
    BB = bwboundaries(KCBW3);
    bd = BB{1};
    
    [bdxwC, bdywC] = intrinsicToWorld(RC,  bd(:, 2)+iSliceC-R2-1, bd(:, 1)); 
    COR.bdxwC = bdxwC;
    COR.bdywC = bdywC;

    if bPlot
        figure(11), subplot(223), hold on
        plot(bdxwC, bdywC, 'g', 'LineWidth', 2);
    end

    % ball
    polyin = polyshape(bdxwC, bdywC);
    [xlim, ylim] = boundingbox(polyin);
    R2 = diff(xlim)/2;
    ind = find(bdywC < ylim(1)+R2*3/2);
    bdx2 = bdxwC(ind);
    bdy2 = bdywC(ind);
    [xc2, yc2, R2, a] = circfit(bdx2, bdy2);
    COR.xc2 = xc2;
    COR.yc2 = yc2;

    theta = linspace(0, pi*2, 100);
    xxc = xc2+R2*cos(theta);
    yyc = yc2+R2*sin(theta);
    COR.xxc = xxc;
    COR.yyc = yyc;
    COR.R = R2;
    
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

    %% segittal gourd

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
    
    [bdxwS, bdywS] = intrinsicToWorld(RS,  bd(:, 2)+iSliceS-R2-1, bd(:, 1)); 
    SEG.bdxwS = bdxwS;
    SEG.bdywS = bdywS;
    
    % ball
    polyin = polyshape(bdxwS, bdywS);
    [xlim, ylim] = boundingbox(polyin);
    R2 = diff(xlim)/2;
    ind = find(bdywS < ylim(1)+R2*3/2);
    bdx2 = bdxwS(ind);
    bdy2 = bdywS(ind);
    [xc2, yc2, R2, a] = circfit(bdx2, bdy2);
    
    theta = linspace(0, pi*2, 100);
    xxc = xc2+R2*cos(theta);
    yyC = yc2+R2*sin(theta);

    SEG.R = R2;
    SEG.xc2 = xc2;
    SEG.yc2 = yc2;
    SEG.xxc = xxc;
    SEG.yyc = yyc;

end

%%%%%%%%%%%%%%%%%%%%%%%%%
function [BW3] = fun_findTank(I)
    J = rescale(I);
    ImageSize = size(J);
    BW = imbinarize(J);
    BW2 = imfill(BW, 'holes');
    BW3 = bwareaopen(BW2, round(ImageSize(1)*ImageSize(2)/8));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = fun_findO(I)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%
function   [xc,yc,R,a] = circfit(x,y)
x=x(:); y=y(:);
a=[x y ones(size(x))]\[-(x.^2+y.^2)];
xc = -.5*a(1);
yc = -.5*a(2);
R  =  sqrt((a(1)^2+a(2)^2)/4-a(3));

end
