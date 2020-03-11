function stack = load_stacks(rawPath, src, frame, channel, varargin)


if ~isempty(varargin)
    xDim = varargin{1};
    yDim = varargin{2};
    zDim = varargin{3};
    stack = zeros(yDim, xDim, zDim, 'uint16');%NaN(size(mcp_stack));
else
    stack = [];%NaN(size(mcp_stack));
end

files = dir([rawPath src '/*_' sprintf('%03d',frame) '*_ch0' num2str(channel) '.tif']);
if ~isempty(dir)
    
    for im = 2:numel(files)-1
        stack(:,:,im-1) = imread([rawPath src '/' files(im).name]);
    end
    stack = double(stack);
    
elseif exist([rawPath, src, filesep, src, '_movieMat.mat'], 'file')
   
    warning('tif stacks not found. using movie stacks.');
    stack = loadMovieMat([rawPath, src, filesep, src, '_movieMat.mat'],'frameRange', frame, 'chRange', channel);

else
    
    error('no movie stacks found');
    
end