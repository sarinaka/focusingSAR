function [] = prettyPlot(X)
    colormap gray
    X = flipud(X);
    h = pcolor(X);
    h.FaceColor = 'interp';
    h.EdgeColor = 'none';
end