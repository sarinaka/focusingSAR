function [] = prettyPlot(varargin)
% function [] = prettyPlot(varargin) plots pretty versions of sounder data
% [] = prettyPlot(X) plots matrix X
% [] = prettyPlot(x,y,X) plots matrix X on axis x,y
colormap(flipud(gray))
    if(nargin == 1)
        X = flipud(varargin{1});
        h = pcolor(X);
      h.FaceColor = 'interp';
        h.EdgeColor = 'none';
    else if(nargin == 3)
        X = flipud(varargin{3});
        y = flip(varargin{2});
%         y = (varargin{2});
        h = pcolor(varargin{1},y,X);
        h.FaceColor = 'interp';
        h.EdgeColor = 'none';
    end
set(gca, 'YDir','reverse')

end