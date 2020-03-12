function [qualityControlStructure, nucleus_struct] =...
    ...
    main02subfunction(...
    rawPath, set_key, nucleusXVector, nucleusYVector, spot_x_vec, spot_y_vec, spot_z_vec,...
    spot_x_vec3D, spot_y_vec3D, spot_z_vec3D, proteinQualityControlVector, nc_ref_frame, nucleusMasterVector,...
    display_figures, spotRoiFrame,...
    xReferenceMesh,yReferenceMesh,zReferenceMesh,xy_sigma_px,z_sigma_px,...
    proteinSnippetSize_px, nucleusDistanceFrame, minSampleSeparation_px,...
    roiRadiusSpot_px, nucleus_struct, nc_lin_index_vec,  nc_sub_index_vec, minEdgeSeparation_px,...
    driftTolerance,PixelSize_um, setID, frame, particle_id_vec,...
    neighborhoodSize_px, minArea_px, maxArea_px,...
    mf_samp_rad, spot_dist_frame,  newSnippetFields, newVectorFields, DataPath, write_snip_flag, jIndex)

xDim = nucleus_struct(1).xDim;
yDim = nucleus_struct(1).yDim;
zDim = nucleus_struct(1).zDim;

source = set_key(set_key.setID==setID,:).prefix{1};

mcp_stack = load_stacks(rawPath, source, frame, mcp_channels(setID), xDim, yDim, zDim);
protein_stack = load_stacks(rawPath, source, frame, protein_channels(setID), xDim, yDim, zDim);

nNuclei = numel(nucleusXVector);
qc_mat = struct;


% initialize arrays to store relevant info
for i = 1:numel(newVectorFields)
    eval([newVectorFields{i} ' = NaN(size(spot_x_vec));']);
end
for j = 1:numel(newSnippetFields)
    eval([newSnippetFields{j} ' = NaN(2*proteinSnippetSize_px+1,2*proteinSnippetSize_px+1,numel(spot_x_vec));']);
end


% generate lookup table of inter-nucleus distances
xInterNuclearDistanceMatrix = repmat(nucleusXVector,nNuclei,1)-repmat(nucleusXVector',1,nNuclei);
yInterNuclearDistanceMatrix = repmat(nucleusYVector,nNuclei,1)-repmat(nucleusYVector',1,nNuclei);
interNuclearDisplacementMatrix = sqrt(double(xInterNuclearDistanceMatrix).^2 + double(yInterNuclearDistanceMatrix).^2);


% iterate through spots
qc_mat = struct;
for spotIndex = 1:nNuclei
    % get location info
    x_nucleus = round(nucleusXVector(spotIndex));
    y_nucleus = round(nucleusYVector(spotIndex));
    x_spot = round(spot_x_vec(spotIndex));
    y_spot = round(spot_y_vec(spotIndex));
    z_spot = round(spot_z_vec(spotIndex))-1;
    % get position info from 3D fit
    x_spot3D = spot_x_vec3D(spotIndex);
    y_spot3D = spot_y_vec3D(spotIndex);
    z_spot3D = spot_z_vec3D(spotIndex)-1;
    
    if isnan(x_spot) || ~proteinQualityControlVector(spotIndex)
        continue
    end
    
    % extract mask
    spot_nc_mask =  nc_ref_frame == nucleusMasterVector(spotIndex);
    if display_figures, imagescUpdate(nexttile, spot_nc_mask, []); drawnow; end
    
    
    %%%%%%%%%%%%%%%%%%%%%% Sample protein levels %%%%%%%%%%%%%%%%%%%%%%
    % get frames
    protein_frame = protein_stack(:,:,z_spot);
    if display_figures, imagescUpdate(nexttile, protein_frame, []); drawnow; end
    mcp_frame = mcp_stack(:,:,z_spot);
    if display_figures, imagescUpdate(nexttile, mcp_frame, []); drawnow; end
    int_id = spotRoiFrame(y_spot,x_spot);
    
    %         % regardless of quality control issues filtered for later on, take protein samples
    %         % in vicinity of spot
    spot_protein_vec(spotIndex) = nanmean(protein_frame(int_id==spotRoiFrame));
    spot_mcp_vec(spotIndex) = nanmean(mcp_frame(int_id==spotRoiFrame));
    
    % volume protein sampling
    spot_protein_vec_3d(spotIndex) = sample_protein_3D(x_spot3D,y_spot3D,z_spot3D,...
        xReferenceMesh,yReferenceMesh,zReferenceMesh,xy_sigma_px,z_sigma_px,protein_stack);
    
    
    % make sure size is reasonable and that spot is inside nucleus
    %beware- continue ahead
    isTooSmall = sum(spot_nc_mask(:)) < minArea_px;
    isTooBig = sum(spot_nc_mask(:)) > maxArea_px;
    isOutsideNucleus = ~spot_nc_mask(y_spot,x_spot);
    
    if isTooSmall || isTooBig || isOutsideNucleus%|| int_it==0
        edge_qc_flag_vec(spotIndex) = -1;
        serial_qc_flag_vec(spotIndex) = -1;
        continue
    end
    
    
    % sample snippets
    spot_protein_snips(:,:,spotIndex) = sample_snip(x_spot,y_spot,proteinSnippetSize_px,protein_frame,spot_nc_mask);
    spot_mcp_snips(:,:,spotIndex) = sample_snip(x_spot,y_spot,proteinSnippetSize_px,mcp_frame,spot_nc_mask);
    
    % Take average across all pixels within 1.5um of nucleus center
    dist_mat = bwdist(~spot_nc_mask);
    if display_figures, imagescUpdate(nexttile, dist_mat, []); drawnow; end
    
    meanFreeSampleMask = dist_mat*PixelSize_um >= mf_samp_rad;
    if display_figures, imagescUpdate(nexttile, meanFreeSampleMask, []); drawnow; end
    
    mf_null_protein_vec(spotIndex) = nanmean(protein_frame(meanFreeSampleMask));% / voxel_size;
    % Edge sampling
    spot_edge_dist = nucleusDistanceFrame(y_spot,x_spot);
    nc_edge_dist_vec = nucleusDistanceFrame(spot_nc_mask);
    spot_edge_dist_vec(spotIndex) = spot_edge_dist;
    spot_sep_vec = spot_dist_frame(spot_nc_mask);
    nc_indices = find(spot_nc_mask);
    
    % Now find control "spot" that is same distance from nucleus edge
    % as true spot
    [edge_null_x_vec(spotIndex), edge_null_y_vec(spotIndex), edge_null_nc_vec(spotIndex), edge_qc_flag_vec(spotIndex),~]...
        = find_control_sample(nc_edge_dist_vec, xReferenceMesh, yReferenceMesh, spot_sep_vec, spot_edge_dist,...
        spotIndex, minSampleSeparation_px, spot_nc_mask,0);
    
    
    
    % if initial attempt failed, try nearest neighbor nucleus
    null_mask = spot_nc_mask;
    if edge_qc_flag_vec(spotIndex) == 0
        % Find nearest neighbor nucleus
        interNuclearDisplacementVector = interNuclearDisplacementMatrix(:,spotIndex);
        interNuclearDisplacementVector(spotIndex) = Inf;
        [~, mi] = min(interNuclearDisplacementVector);
        
        % get nn nucleus mask
        x_spot_nn = round(spot_x_vec(mi));
        y_spot_nn = round(spot_y_vec(mi));
        
        nn_nc_mask = nc_ref_frame == nucleusMasterVector(mi);
        if display_figures, imagescUpdate(nexttile, nn_nc_mask, []); drawnow; end
        
        null_mask = nn_nc_mask; % reassign null mask
        if ~isnan(x_spot_nn)
            nan_flag = isnan(nn_nc_mask(y_spot_nn,x_spot_nn));
        end
        % make sure size is reasonable
        if sum(nn_nc_mask(:)) >= minArea_px && sum(nn_nc_mask(:)) <= maxArea_px
            nn_edge_dist_vec = nucleusDistanceFrame(nn_nc_mask);
            nn_sep_vec = spot_dist_frame(nn_nc_mask);
            [edge_null_x_vec(spotIndex), edge_null_y_vec(spotIndex), edge_null_nc_vec(spotIndex), edge_qc_flag_vec(spotIndex),~]...
                = find_control_sample(nn_edge_dist_vec, xReferenceMesh, yReferenceMesh, nn_sep_vec, spot_edge_dist,...
                mi, minSampleSeparation_px, null_mask,1);
            if  edge_qc_flag_vec(spotIndex) == 1
                edge_qc_flag_vec(spotIndex) =  2; % to flag cases when nn was used
            end
        end
    end
    
    
    
    
    % Draw control samples (as appropriate)
    if edge_qc_flag_vec(spotIndex) > 0
        
        xControl = edge_null_x_vec(spotIndex);
        yControl = edge_null_y_vec(spotIndex);
        
        null_dist_frame = zeros(size(protein_frame));
        null_dist_frame(yControl,xControl) = 1;
        null_dist_frame = bwdist(null_dist_frame);
        if display_figures, imagescUpdate(nexttile, null_dist_frame, []); drawnow; end
        
        
        edge_null_protein_vec(spotIndex) = nanmean(protein_frame(nc_ref_frame>0&null_dist_frame<roiRadiusSpot_px));% / voxel_size;
        edge_null_mcp_vec(spotIndex) = nanmean(mcp_frame(nc_ref_frame>0&null_dist_frame<roiRadiusSpot_px));% / voxel_size;
        
        % take 3D protein sample
        edge_null_protein_vec_3d(spotIndex) = sample_protein_3D(xControl,yControl,z_spot3D,xReferenceMesh,yReferenceMesh,zReferenceMesh,xy_sigma_px,z_sigma_px,protein_stack);
        
        % draw snips
        edge_null_protein_snips(:,:,spotIndex) = sample_snip(xControl,yControl,proteinSnippetSize_px,protein_frame,nc_ref_frame>0);
        edge_null_mcp_snips(:,:,spotIndex) = sample_snip(xControl,yControl,proteinSnippetSize_px,mcp_frame,nc_ref_frame>0);
        
    end
    
    
    
    % Draw serialized control
    nc_index = nc_lin_index_vec(spotIndex);
    nc_sub_index = nc_sub_index_vec(spotIndex);
    frame_vec_temp = nucleus_struct(nc_index).frames;
    serial_null_x = nucleus_struct(nc_index).serial_null_x_vec;
    serial_null_y = nucleus_struct(nc_index).serial_null_y_vec;
    % if this is the first sample for this spot, just find random
    % control snip. This will "seed" subsequent samples
    if all(isnan(serial_null_x))
        % Now take a random sample
        sample_index_vec = 1:numel(spot_sep_vec);
        % filter for regions far enough away from locus
        cr_filter = spot_sep_vec >= minSampleSeparation_px & nc_edge_dist_vec >= minEdgeSeparation_px;
        sample_index_vec = sample_index_vec(cr_filter);
        % if candidate found, then proceed. Else look to neighboring nuclei
        if ~isempty(sample_index_vec)
            new_index = randsample(sample_index_vec,1);
            x_pos_vec = xReferenceMesh(spot_nc_mask);
            y_pos_vec = yReferenceMesh(spot_nc_mask);
            xControl = x_pos_vec(new_index);
            yControl = y_pos_vec(new_index);
            ec = nc_edge_dist_vec(new_index);
        else
            error('Unable to draw random sample. Check "PixelSize" and "min_sample_sep" variables')
        end
        % otherwise, draw snip based on previous location
    else
        prev_frame = find(~isnan(serial_null_x),1,'last');
        n_frames = frame - frame_vec_temp(prev_frame); % used to adjust jump weights
        old_x = double(serial_null_x(prev_frame));
        old_y = double(serial_null_y(prev_frame));
        % possible locations
        x_pos_vec = xReferenceMesh(spot_nc_mask);
        y_pos_vec = yReferenceMesh(spot_nc_mask);
        drControl = double(sqrt((old_x-x_pos_vec).^2+(old_y-y_pos_vec).^2));
        %             edge_dev_vec = spot_edge_dist-nc_edge_dist_vec;
        % calculate weights
        wt_vec = exp(-.5*((drControl/double((n_frames*driftTolerance/PixelSize_um))).^2));%+((edge_dev_vec)/roi_spot).^2));
        % anything too close to locus or with an edge distance too different from locus is excluded
        wt_vec(spot_sep_vec<minSampleSeparation_px|nc_edge_dist_vec < minEdgeSeparation_px) = 0;
        % draw sample
        xControl = NaN;
        yControl = NaN;
        ec = NaN;
        if any(wt_vec>0)
            new_index = randsample(1:numel(x_pos_vec),1,true,wt_vec);
            xControl = x_pos_vec(new_index);
            yControl = y_pos_vec(new_index);
            ec = nc_edge_dist_vec(new_index);
        else
            new_index = randsample(1:numel(x_pos_vec),1,true);
            xControl = x_pos_vec(new_index);
            yControl = y_pos_vec(new_index);
            ec = nc_edge_dist_vec(new_index);
            %             else
            %                 error('This should not happen')
        end
    end
    % draw samples
    null_dist_frame = zeros(size(protein_frame));
    if ~isnan(xControl)
        null_dist_frame(yControl,xControl) = 1;
    end
    null_dist_frame = bwdist(null_dist_frame);
    % samples below default to NaN if no sample taken
    % sample protein
    serial_qc_flag_vec(spotIndex) = ~isnan(xControl);
    serial_null_protein_vec(spotIndex) = nanmean(protein_frame(nc_ref_frame>0&null_dist_frame<roiRadiusSpot_px));% / voxel_size;
    % record
    serial_null_x_vec(spotIndex) = xControl;
    serial_null_y_vec(spotIndex) = yControl;
    serial_null_edge_dist_vec(spotIndex) = ec;
    % 3D version
    serial_null_protein_vec_3d(spotIndex) = sample_protein_3D(xControl,yControl,z_spot3D,...
        xReferenceMesh,yReferenceMesh,zReferenceMesh,xy_sigma_px,z_sigma_px,protein_stack);% / sum(vol_denominator(:)) / voxel_size;
    % check for presence of sister spot
    x_spot_sister = NaN;
    y_spot_sister = NaN;
    if sum(nucleusMasterVector(spotIndex)==nucleusMasterVector) == 2
        indices = find(nucleusMasterVector(spotIndex)==nucleusMasterVector);
        si = indices(indices~=spotIndex);
        x_spot_sister = spot_x_vec(si);
        y_spot_sister = spot_y_vec(si);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % save quality control data
    q = nNuclei-spotIndex+1;
    
    qc_mat(q).setID = setID;
    qc_mat(q).frame = frame;
    qc_mat(q).nc_index = nc_lin_index_vec(spotIndex);
    qc_mat(q).nc_sub_index = nc_sub_index_vec(spotIndex);
    qc_mat(q).qc_flag = edge_qc_flag_vec(spotIndex);
    qc_mat(q).xp = x_spot;
    qc_mat(q).yp = y_spot;
    qc_mat(q).xp_sister = x_spot_sister;
    qc_mat(q).yp_sister = y_spot_sister;
    qc_mat(q).xc_edge = edge_null_x_vec(spotIndex);
    qc_mat(q).yc_edge = edge_null_y_vec(spotIndex);
    qc_mat(q).xc_serial = serial_null_x_vec(spotIndex);
    qc_mat(q).yc_serial = serial_null_y_vec(spotIndex);
    qc_mat(q).ParticleID = particle_id_vec(spotIndex);
    qc_mat(q).serial_qc_flag = serial_qc_flag_vec(spotIndex);
    qc_mat(q).edge_qc_flag = edge_qc_flag_vec(spotIndex);
    sz = neighborhoodSize_px;
    edge_dist_mat = nucleusDistanceFrame;
    edge_dist_mat(~spot_nc_mask&~null_mask) = 0;
    if edge_qc_flag_vec(spotIndex) == 2
        sz = max([neighborhoodSize_px,abs(x_nucleus - edge_null_x_vec(spotIndex)),abs(y_nucleus - edge_null_y_vec(spotIndex))...
            abs(x_nucleus - serial_null_x_vec(spotIndex)),abs(y_nucleus - serial_null_y_vec(spotIndex))]);
    end
    y_range = max(1,y_nucleus-sz):min(yDim,y_nucleus+sz);
    x_range = max(1,x_nucleus-sz):min(xDim,x_nucleus+sz);
    qc_mat(q).x_origin = x_range(1);
    qc_mat(q).y_origin = y_range(1);
    qc_mat(q).mcp_snip = mcp_frame(y_range,x_range);
    qc_mat(q).protein_snip = protein_frame(y_range,x_range);
    qc_mat(q).edge_dist_snip = edge_dist_mat(y_range,x_range);
    
end

qualityControlStructure(jIndex).qc_mat = fliplr(qc_mat);

%%

% save snip data
nNuclei = numel(nucleusMasterVector);
% initialize struct to store snip data
snip_data = struct;
% map data back to nucleus_struct
for indo = 1:nNuclei
    nc_index = nc_lin_index_vec(indo);
    nc_sub_index = nc_sub_index_vec(indo);
    frame = nucleus_struct(nc_index).frames(nc_sub_index);
    for k = 1:numel(newVectorFields)
        vec = eval(newVectorFields{k});
        nucleus_struct(nc_index).(newVectorFields{k})(nc_sub_index) = vec(indo);
    end
    % store snips
    for k = 1:numel(newSnippetFields)
        snip = eval([newSnippetFields{k} '(:,:,indo)']);
        snip_data.(newSnippetFields{k})(:,:,indo) = snip;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% store key ID variables
snip_data.frame = frame;
snip_data.setID = setID;
snip_data.particle_id_vec = particle_id_vec;
snip_data.spot_edge_dist_vec = spot_edge_dist_vec;
% indexing vectors
snip_data.nc_sub_index_vec = nc_sub_index_vec;
snip_data.nc_lin_index_vec = nc_lin_index_vec;
snip_data.nc_master_vec = nucleusMasterVector;
% specify name
% read snip file
snip_name = ['snip_data_F' sprintf('%03d',frame) '_S' sprintf('%02d',setID)];
if write_snip_flag
    blank = struct;
    save([DataPath 'snip_data.mat'],'blank','-v7.3')
    snip_file = matfile([DataPath 'snip_data.mat'],'Writable',true);
    snip_file.(snip_name)= snip_data;
    clear snip_file;
end


end