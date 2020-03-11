function segmentNuclei_main02(yDim, xDim, segment_indices, set_frame_array, set_ref, frame_ref, nc_x_ref, ...
    nc_y_ref, master_ind_ref, spot_x_ref, spot_y_ref, set_key, rawPath, proteinChannel, ...
    sm_kernel, nb_size, refPath, display_figures, min_area, max_area, DropboxFolder, varargin)


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

if display_figures
    
    figure('Units', 'normalized', 'Position',[0.6441 0.0744 0.3184 0.3844]);
    tiledlayout('flow', 'TileSpacing', 'none', 'Padding', 'none')
    nPlots = 3;
    ax = {};
    for i = 1:nPlots
        ax{i} = nexttile;
    end
    drawnow;
    
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
    
    %% AR Addition- let's just operate on a denoised image from here out
    
    nuclear_image = wiener2(nuclear_image);
    
    %%
    if display_figures
        imagescUpdate(ax{1}, nuclear_image, []);  
        drawnow;
    end
    
    
%     %% Original Masking Method
%     protein_smooth = imgaussfilt(nuclear_image,round(sm_kernel/2));
%     protein_grad = imgradient(protein_smooth);
%     % flatten background
%     protein_bkg = imgaussfilt(protein_smooth, round(nb_size/2));
%     protein_grad_norm = protein_grad ./ protein_bkg;
%     
%     
%     
%     % try local thresholding
%     n_x_local = floor(xDim / nb_size);
%     n_y_local = floor(yDim / nb_size);
%     x_dim_local = round(xDim / n_x_local);
%     y_dim_local = round(yDim / n_y_local);
%     protein_bin_clean = false(size(protein_grad_norm));
%     for x = 1:n_x_local
%         for y = 1:n_y_local
%             x_start = (x-1)*x_dim_local+1;
%             x_stop = min([x*x_dim_local,xDim]);
%             y_start = (y-1)*y_dim_local+1;
%             y_stop = min([y*y_dim_local,yDim]);
%             pt_section = protein_grad_norm(y_start:y_stop,x_start:x_stop);
%             thresh = multithresh(pt_section);
%             section_bin = pt_section > thresh;
%             section_bin_clean = bwareaopen(section_bin,sm_kernel^2);
%             % record
%             protein_bin_clean(y_start:y_stop,x_start:x_stop) = bwmorph(section_bin_clean,'hbreak');
%         end
%     end
%     
    %% AR MASKING METHOD
    sigma = 7; b = -.4; s = 0;
    %assuming the right class is 3 here. sometimes that might be wrong.
    %it's hard to decide which one is right automatically. 
    kmask = imsegkmeans(single(imgaussfilt(nuclear_image,sigma)),3)==3;
    snakesFun = @(b, s, sigma) activecontour(imgaussfilt(nuclear_image, sigma), kmask, 'Chan-Vese', 'ContractionBias', b, 'SmoothFactor', s);
    kMaskRefined = snakesFun(b, s, sigma/2);
    %same issue here. sometimes the watershed is inverted. hard to pick
    %whether it should be or not automatically
    protein_bin_clean = bwareafilt(wshed(kMaskRefined), [min_area, max_area]);
    
    %fit with circles instead of convex hulls
    protein_bin_clean = fitCirclesToNuclei(protein_bin_clean);
    
    
    %%
    
    %%%%%%%%%%%%%diagnostic display
    if display_figures
        imagescUpdate(ax{2}, protein_bin_clean, []);
        drawnow;
    end
    %%%%%%%%%%%%%%%%%%%%%
    
    
    % label regions
%     nc_frame_labels = logical(protein_bin_clean);
    
    nc_frame_labels = bwlabel(protein_bin_clean);
    
    % frame info
%     nc_lin_indices = sub2ind(size(protein_bin_clean),round(nc_y_vec_u),round(nc_x_vec_u));
    % take convex hull
%     stats = regionprops(nc_frame_labels,'ConvexHull');
% %         nc_frame_labels = imfill(nc_frame_labels, 'holes');
%     nc_ref_frame = zeros(size(nc_frame_labels));
%     maskTotal = zeros(yDim, xDim); %for diagnostics
    
%     maskTotal = protein_bin_clean;
    nc_ref_frame = protein_bin_clean;

       
%     for j = 1:numel(stats)
        
%         hull_points = stats(j).ConvexHull;
%         mask = poly2mask(hull_points(:,1),hull_points(:,2),yDim,xDim);
%             
%             mask = nc_frame_labels == j;

%         mask = bwareafilt(mask, [min_area, max_area]); %throw out small oversegmentation products
    
%         maskTotal = maskTotal + mask;
        
%         nc_bin_ids = mask(nc_lin_indices);
%         if sum(nc_bin_ids) == 1 % enforce unique
%             nc_ref_frame(mask) = nc_master_vec_u(nc_bin_ids);
%         end
%     end
    
    %%%%%%%%%%%%%diagnostic display
    if display_figures
        imagescUpdate(ax{3}, nc_ref_frame, []);
        drawnow;
    end
    %%%%%%%%%%%%%
    
    
    % generate array indicating distance of each pixel from an active locus
    nc_indices = sub2ind(size(nc_ref_frame),spot_y_vec,spot_x_vec);
    spot_dist_frame_temp = zeros(size(nc_ref_frame));
    spot_dist_frame_temp(nc_indices(~isnan(nc_indices))) = 1;
    spot_dist_frame_temp = bwdist(spot_dist_frame_temp);
    
    
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