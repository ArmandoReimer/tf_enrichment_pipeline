function snakesTile(fun, paramRange, nTiles, varargin)

cmap = 'hsv';

for i = 1:2:(numel(varargin)-1)
    if i ~= numel(varargin)
        eval([varargin{i} '=varargin{i+1};']);
    end
end

if iscell(paramRange)
    nParams = numel(paramRange);
else
    nParams = 1;
end

figure('Units', 'normalized', 'Position', [.25, .25, .5, .5]);
tiledlayout('flow', 'TileSpacing', 'none', 'Padding', 'none')
colormap(cmap);

params = cell(1, nParams);
for k = 1:nParams
    params{k} = linspace(paramRange{k}(1),paramRange{k}(2), nTiles);
end

for b = params{1}
    for s = params{2}

        ax = nexttile;
        im = fun(b, s);
        imagescUpdate(ax,im, []);
        colorbar;

        ax.Visible = 'off';
        title(ax, ['bias=',num2str(b),'; $\sigma$=', num2str(s)], 'Interpreter', 'latex');
        set(findall(ax, 'type', 'text'), 'Visible', 'on')
        axis image

        drawnow;

    end
end

end