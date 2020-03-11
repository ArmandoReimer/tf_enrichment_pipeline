function segmentNuclei_main02(yDim, xDim, segment_indices, set_frame_array, set_ref, frame_ref, nc_x_ref, ...
    nc_y_ref, master_ind_ref, spot_x_ref, spot_y_ref, set_key, rawPath, proteinChannel, ...
    sm_kernel, nb_size, refPath, displayFigures, min_area, max_area, DropboxFolder, varargin)

maskingMethod = 'kSnakeCircles';
% maskingMethod = 'gradientOtsuHulls';
%

nSets = size(set_key, 1);
nuclear_mov = cell(1, nSets);

if ~isempty(varargin)
    
    nuclear_mov = varargin{1};
    
else
    
    for set = 1:nSets
        prefix = cell2mat(set_key{set, 2});
        hisPath = [rawPath, prefix, filesep, prefix, '_hisMat.mat'];
        if exist(hisPath, 'file')
            %         load(hisPath, 'hisMat'); % t y x set
            hisMat = loadHisMat(hisPath);
        else
            load([DropboxFolder, filesep, prefix, filesep, 'FrameInfo.mat'], 'FrameInfo')
            nWorkers = 1; [~, hisMat, ~, ~, ~] = makeMovieMats(prefix, rawPath, nWorkers, FrameInfo);
        end
        nuclear_mov{set} = double(hisMat); % {set(t y x)}
    end
    clear hisMat;
end

if displayFigures
    
    figure('Units', 'normalized', 'Position',[0.6441 0.0744 0.3184 0.3844]);
    tiledlayout('flow', 'TileSpacing', 'none', 'Padding', 'none')
    
end

disp('segmenting nuclei...')
tic
%%% Segment nuclei
% initialize arrays to store segmentation info
nucleus_frame_array = NaN(yDim,xDim,numel(segment_indices));
spot_frame_array = NaN(yDim,xDim,numel(segment_indices));
for w = 1:numel(segment_indices)
    i = segment_indices(w);
    setID_temp = set_frame_array(i,1);
    frame_temp = set_frame_array(i,2);
    % get nucleus
    frame_set_filter = set_ref==setID_temp&frame_ref==frame_temp;
    nc_x_vec_temp = round(nc_x_ref(frame_set_filter));
    nc_y_vec_temp = round(nc_y_ref(frame_set_filter));
    % indexing vectors
    nc_master_vec = master_ind_ref(frame_set_filter);
    % get list of unique indices
    [nc_master_vec_u,ia,~] = unique(nc_master_vec,'stable');
    % unique nucleus vectors
    nc_x_vec_u = round(nc_x_vec_temp(ia));
    nc_y_vec_u = round(nc_y_vec_temp(ia));
    % particle positions
    spot_x_vec = round(spot_x_ref(frame_set_filter));
    spot_y_vec = round(spot_y_ref(frame_set_filter));
    
    
    % generate protein gradient frame for segmentation
    
    src = set_key(set_key.setID==setID_temp,:).prefix{1};
    %         nuclear_stack = load_stacks(rawPath, src, frame_temp, proteinChannel);
    %         nuclear_image = mean(nuclear_stack, 3);
    
    nuclear_image = nuclear_mov{setID_temp}(:, :, frame_temp);
    
    yDim = size(nuclear_image, 1);
    xDim = size(nuclear_image, 2);
    
    nc_lin_indices = sub2ind([yDim, xDim],round(nc_y_vec_u),round(nc_x_vec_u));
    
    %% AR Addition- let's just operate on a denoised image from here out
    nuclear_image = wiener2(nuclear_image);
    
    nc_ref_frame = generateNuclearMask(nuclear_image,...
        'maskingMethod', maskingMethod, 'areaFilter', [min_area, max_area],...
        'nc_lin_indices', nc_lin_indices, 'nb_size', nb_size,...
        'sm_kernel', sm_kernel, 'nc_master_vec_u', nc_master_vec_u);
    
    if displayFigures, imagescUpdate(nexttile, nuclear_image, []); drawnow; end
    if displayFigures, imagescUpdate(nexttile, nc_ref_frame, []); drawnow; end
    
    % generate array indicating distance of each pixel from an active locus
    nc_indices = sub2ind(size(nc_ref_frame),spot_y_vec,spot_x_vec);
    spot_dist_frame_temp = zeros(size(nc_ref_frame));
    spot_dist_frame_temp(nc_indices(~isnan(nc_indices))) = 1;
    spot_dist_frame_temp = bwdist(spot_dist_frame_temp);
    
    if displayFigures, imagescUpdate(nexttile, spot_dist_frame_temp, []); drawnow; end
    
    % label regions within designated integration radius of a spot
    % store arrays
    nucleus_frame_array(:,:,w) = nc_ref_frame;
    spot_frame_array(:,:,w) = spot_dist_frame_temp;
    % save spot and nucleus reference frames
    
end

disp('saving segmentation results...')

% save arrays
for w = 1:numel(segment_indices)
    i = segment_indices(w);
    setID_temp = set_frame_array(i,1);
    frame_temp = set_frame_array(i,2);
    nc_ref_name = [refPath 'nc_ref_frame_set' sprintf('%02d',setID_temp) '_frame' sprintf('%03d',frame_temp) '.mat'];
    nc_ref_frame = nucleus_frame_array(:,:,w);
    save(nc_ref_name,'nc_ref_frame');
    spot_ref_name = [refPath 'spot_roi_frame_set' sprintf('%02d',setID_temp) '_frame' sprintf('%03d',frame_temp) '.mat'];
    spot_dist_frame = spot_frame_array(:,:,w);
    save(spot_ref_name,'spot_dist_frame');
end
disp('done.')
toc

end