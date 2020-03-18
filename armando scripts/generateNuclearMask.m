function nuclearMask = generateNuclearMask(nuclearImage, varargin)

maskingMethod = 'gradientOtsuHulls';
displayFigures = false;
yDim = size(nuclearImage, 1);
xDim = size(nuclearImage, 2);
areaFilter = [0, Inf];
nb_size = [];
sm_kernel  = [];
nc_lin_indices = []; %indices for nuclear centroids from Ellipses.mat
nc_master_vec_u = [];



%options must be specified as name, value pairs. unpredictable errors will
%occur, otherwise.
for i = 1:2:(numel(varargin)-1)
    if i ~= numel(varargin)
        eval([varargin{i} '=varargin{i+1};']);
    end
end

if displayFigures
    figure('Units', 'normalized', 'Position',[0.6441 0.0744 0.3184 0.3844]);
    tiledlayout('flow', 'TileSpacing', 'none', 'Padding', 'none')
end

switch maskingMethod
    case 'kSnakeCircles'
        
        sigma = 7; b = -.4; s = 0;
        %assuming the right class is 3 here. sometimes that might be wrong.
        %it's hard to decide which one is right automatically.
        kmask = imsegkmeans(single(imgaussfilt(nuclearImage,sigma)),3)==3;
        snakesFun = @(b, s, sigma) activecontour(imgaussfilt(nuclearImage, sigma), kmask, 'Chan-Vese', 'ContractionBias', b, 'SmoothFactor', s);
        kMaskRefined = snakesFun(b, s, sigma/2);
        %same issue here. sometimes the watershed is inverted. hard to pick
        %whether it should be or not automatically
        nuclearMask = bwareafilt(wshed(kMaskRefined), areaFilter);
        
        %fit with circles instead of convex hulls
        protein_bin_clean = fitCirclesToNuclei(nuclearMask);
        
        lab = bwlabel(protein_bin_clean);
        
        nuclearMask = zeros(yDim, xDim);
        for j = 1:max(lab(:))
                        
            mask = lab==j;
            
            nc_bin_ids = mask(nc_lin_indices);
            if sum(nc_bin_ids) == 1 % enforce unique
                nuclearMask(mask) = nc_master_vec_u(nc_bin_ids);
            end
            
        end
            
    case 'gradientOtsuHulls'
        
        %% Original Masking Method
        protein_smooth = imgaussfilt(nuclearImage,round(sm_kernel/2));
        protein_grad = imgradient(protein_smooth);
        % flatten background
        protein_bkg = imgaussfilt(protein_smooth, round(nb_size/2));
        protein_grad_norm = protein_grad ./ protein_bkg;
        
        % try local thresholding
        n_x_local = floor(xDim / nb_size);
        n_y_local = floor(yDim / nb_size);
        x_dim_local = round(xDim / n_x_local);
        y_dim_local = round(yDim / n_y_local);
        protein_bin_clean = false(size(protein_grad_norm));
        for x = 1:n_x_local
            for y = 1:n_y_local
                x_start = (x-1)*x_dim_local+1;
                x_stop = min([x*x_dim_local,xDim]);
                y_start = (y-1)*y_dim_local+1;
                y_stop = min([y*y_dim_local,yDim]);
                pt_section = protein_grad_norm(y_start:y_stop,x_start:x_stop);
                thresh = multithresh(pt_section);
                section_bin = pt_section > thresh;
                section_bin_clean = bwareaopen(section_bin,sm_kernel^2);
                % record
                protein_bin_clean(y_start:y_stop,x_start:x_stop) = bwmorph(section_bin_clean,'hbreak');
            end
        end
        
        stats = regionprops(~~protein_bin_clean,'ConvexHull');
        nRegions = numel(stats);
        
        nuclearMask = zeros(yDim, xDim);
        for j = 1:nRegions
            
            hull_points = stats(j).ConvexHull;
            mask = poly2mask(hull_points(:,1),hull_points(:,2),yDim,xDim);
            
            mask = bwareafilt(mask, areaFilter); %throw out small oversegmentation products

            
            nc_bin_ids = mask(nc_lin_indices);
            if sum(nc_bin_ids) == 1 % enforce unique
                nuclearMask(mask) = nc_master_vec_u(nc_bin_ids);
            end
        end
        
    otherwise
        error('what?');
end


end