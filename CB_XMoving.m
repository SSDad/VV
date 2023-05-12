function CB_XMoving(src, evnt)

hFig = ancestor(src, 'Figure');
data = guidata(hFig);

pos = src.Position;
hPlotObj = data.Panel.View.Comp.hPlotObj;

if strcmp(src.Tag, 'A')
    x = pos(1);
    y = pos(2);
    [n, m] = worldToIntrinsic(data.Image.RA, x, y);

    % update Coronal
    p3 = hPlotObj(3).X.Position;
    hPlotObj(3).X.Position = [x p3(2)];

    iSliceC = round(m);
    IC = squeeze(data.Image.V(iSliceC, :, :));
    IC = IC';
    hPlotObj(3).Image.CData = IC;

    % update Sagittal
    p4 = hPlotObj(4).X.Position;
    hPlotObj(4).X.Position = [y p4(2)];

    iSliceS = round(n);
    IS = squeeze(data.Image.V(:, iSliceS, :));
    IS = IS';
    hPlotObj(4).Image.CData = IS;

elseif strcmp(src.Tag, 'C')
    x = pos(1);
    z = pos(2);
    [n, p] = worldToIntrinsic(data.Image.RC, x, z);

    % update Axial
    p1 = hPlotObj(1).X.Position;
    hPlotObj(1).X.Position = [x p1(2)];

    iSliceA = round(p);
    IA = squeeze(data.Image.V(:, :, iSliceA));
    hPlotObj(1).Image.CData = IA;

    % update Sagittal
    p4 = hPlotObj(4).X.Position;
    hPlotObj(4).X.Position = [p4(1) z];

    iSliceS = round(n);
    IS = squeeze(data.Image.V(:, iSliceS, :));
    IS = IS';
    hPlotObj(4).Image.CData = IS;

elseif strcmp(src.Tag, 'S')
    y = pos(1);
    z = pos(2);
    [m, p] = worldToIntrinsic(data.Image.RS, y, z);

    % update Axial
    p1 = hPlotObj(1).X.Position;
    hPlotObj(1).X.Position = [p1(1) y];

    iSliceA = round(p);
    IA = squeeze(data.Image.V(:, :, iSliceA));
    hPlotObj(1).Image.CData = IA;

    % update Coronal
    p3 = hPlotObj(3).X.Position;
    hPlotObj(3).X.Position = [p3(1) z];

    iSliceC = round(m);
    IC = squeeze(data.Image.V(iSliceC, :, :));
    IC = IC';
    hPlotObj(3).Image.CData = IC;

end