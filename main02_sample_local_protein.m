% main02_sample_local_protein(project, RawPath, keyword)
%
% DESCRIPTION
% Script to generate samples of local protein contentration at active gene
% loci and selected control locations
%
% ARGUMENTS
% project: master ID variable 
%
% RawPath: Full or relative path to PreProcessed folder (pr
%               equivalent)
%
%
% OUTPUT: ref_frame_struct: compiled data set

function nucleus_struct_protein = main02_sample_local_protein(project,RawPath,protein_channel,varargin)
tic
zeiss_flag = 0;
DataPath = ['../dat/' project '/'];
for i = 1:numel(varargin)
    if strcmpi(varargin{i}, 'zeiss')
        disp('accounting for zeiss offset')
        zeiss_flag = 1;
    elseif strcmpi(varargin{i}, 'DropboxFolder')        
        DataPath = [varargin{i+1} '/ProcessedEnrichmentData/' project '/'];
    end
end
% Load trace data

load([DataPath '/nucleus_struct.mat'],'nucleus_struct')
load([DataPath '/set_key.mat'],'set_key')
SnipPath = [DataPath '/qc_images/'];
mkdir(SnipPath)
addpath('./utilities')
% get MCP channel
options = [1 2];
mcp_channel = options(options~=protein_channel);
% generate indexing vectors
frame_ref = [nucleus_struct.frames];
nc_x_ref = [nucleus_struct.xPos];
nc_y_ref = [nucleus_struct.yPos];
spot_x_ref = [nucleus_struct.xPosParticle]-zeiss_flag;
spot_y_ref = [nucleus_struct.yPosParticle]-zeiss_flag;
spot_z_ref = [nucleus_struct.brightestZs];

set_ref = [];
ind_ref = [];
sub_ind_ref = [];
pt_ref = [];
for i = 1:numel(nucleus_struct)
    ParticleID = nucleus_struct(i).ParticleID;    
    pt_ref = [pt_ref repelem(ParticleID, numel(nucleus_struct(i).frames))];
    set_ref = [set_ref repelem(nucleus_struct(i).setID,numel(nucleus_struct(i).frames))];
    ind_ref = [ind_ref repelem(i,numel(nucleus_struct(i).frames))];
    sub_ind_ref = [sub_ind_ref 1:numel(nucleus_struct(i).frames)];
end

%%% set snip size to use
set_vec = [nucleus_struct.setID];
set_index = unique(set_vec);
px_sizes = [];
for s = 1:numel(set_index)
    px = [nucleus_struct(set_vec==set_index(s)).PixelSize];
    px_sizes = [px_sizes px(1)];
end% determine size of neighborhood to use
nb_sizes = round(10 ./ px_sizes);
% size of gaussian smoothing kernel 
sm_kernels = round(.2 ./ px_sizes);
% set min and max acceptable area
min_areas = round(pi*(2 ./ px_sizes).^2);
max_areas = round(pi*(4 ./ px_sizes).^2);
% set snippet to be 3um in size
pt_snippet_size_vec = round(1.5 ./ px_sizes);
% set min separation between control and locus to 2um
min_sample_sep_vec = round(2 ./ px_sizes);

% Designate fields ot be added to nucleus structure
new_vec_fields = {'spot_protein_vec', 'edge_null_protein_vec','centroid_null_protein_vec',...
    'spot_mcp_vec','edge_null_mcp_vec','centroid_null_mcp_vec',...
    'edge_qc_flag_vec', 'centroid_qc_flag_vec', 'spot_fov_edge_flag_vec','edge_fov_edge_flag_vec', ...
    'centroid_fov_edge_flag_vec','edge_null_x_vec', 'edge_null_y_vec', 'edge_null_nc_vec',...
    'centroid_null_x_vec', 'centroid_null_y_vec', 'centroid_null_nc_vec',...
    'spot_edge_dist_vec', 'spot_centroid_dist_vec'};

new_snip_fields = {'spot_protein_snips', 'edge_null_protein_snips','centroid_null_protein_snips',...
    'spot_mcp_snips','edge_null_mcp_snips','centroid_null_mcp_snips'};
% get most common snip size
default_snip_size = mode(pt_snippet_size_vec);
% Initialize fields
for i = 1:numel(nucleus_struct)
    ref = nucleus_struct(i).xPos;
    for j = 1:numel(new_vec_fields)
        nucleus_struct(i).(new_vec_fields{j}) = NaN(size(ref));
    end
    for j = 1:numel(new_snip_fields)
        nucleus_struct(i).(new_snip_fields{j}) = NaN(2*default_snip_size+1,2*default_snip_size+1,numel(ref));
    end
end
%%% make source key
src_cell = {1,numel(set_index)};
pt_id_vec = [nucleus_struct.ParticleID];
for s = 1:numel(set_index)
    ind = find([nucleus_struct.setID]==set_index(s)&~isnan(pt_id_vec),1);   
    src = nucleus_struct(ind).source_path;    
    src_cell{s} = src;     
end
%%% Generate reference array for set-frame combos
set_frame_array = unique([set_ref' frame_ref'],'row');
qc_structure = struct;
%%% iterate
for i = 10:15%1:size(set_frame_array,1)
    tic
    setID = set_frame_array(i,1);
    frame = set_frame_array(i,2);       
    % get nucleus and particle positions
    frame_set_filter = set_ref==setID&frame_ref==frame;
    nc_x_vec = nc_x_ref(frame_set_filter);
    nc_y_vec = nc_y_ref(frame_set_filter);        
    spot_x_vec = spot_x_ref(frame_set_filter);
    spot_y_vec = spot_y_ref(frame_set_filter);        
    spot_z_vec = spot_z_ref(frame_set_filter);  
    particle_id_vec = pt_ref(frame_set_filter);
    % find indices for spots
    if sum(~isnan(spot_x_vec)) == 0
        continue
    end
    % indexing variables
    nc_sub_index_vec = sub_ind_ref(frame_set_filter);     
    nc_index_vec = ind_ref(frame_set_filter);    
       
    %%%%% First load MCP frames and generate Nucleus Segmentation frame%%%%
    
    % get size params    
    pt_snippet_size = pt_snippet_size_vec(set_index==setID); % size of snippet
    min_sample_sep = min_sample_sep_vec(set_index==setID); % minimum distance between spot and control
    nb_sz = nb_sizes(set_index==setID); % size of nucleus neighborhood to use
    sm_kernel = sm_kernels(set_index==setID); %sigma for gaussian smoothing kerne;
    min_area = min_areas(set_index==setID); % lower bound on permitted nucleus size
    max_area = max_areas(set_index==setID); % upper bound
    max_r = round(sqrt(max_area/pi))'; % max nucleus neighborhood size
    % determine whether snips will need to be resampled
    snip_scale_factor =  default_snip_size / pt_snippet_size;
    % load and  MCP mCherry and protein stacks
    src = set_key(set_key.setID==setID,:).prefix{1};        
    [mcp_stack, protein_stack] = load_stacks(RawPath, src, frame, mcp_channel);
    
    % Invert image
    mcp_med = mat2gray(nanmedian(mcp_stack,3));
    mcp_med_inv = 1-mcp_med;
    % un-invert pixels around spot center
    spot_x_regions = repmat([spot_x_vec-1 spot_x_vec spot_x_vec+1],1,3);
    spot_y_regions = [repmat(spot_y_vec-1,1,3) repmat(spot_y_vec,1,3) repmat(spot_y_vec+1,1,3)];
    indices = sub2ind(size(mcp_med),spot_y_regions,spot_x_regions);
    indices = indices(~isnan(indices));
    mcp_med_inv(indices) = prctile(mcp_med_inv(:),99);
    
    % smooth and normalize
    mcp_sm = imgaussfilt(mcp_med_inv,sm_kernel);
    his_sm = mcp_sm / mean(mcp_sm(:));    
    
    [id_array, yDim, xDim, y_ref, x_ref] = assign_nc_neighborhoods(his_sm, nc_x_vec, nc_y_vec, max_r, nc_index_vec);
    % generate lookup table of inter-nucleus distances
    x_dist_mat = repmat(nc_x_vec,numel(nc_x_vec),1)-repmat(nc_x_vec',1,numel(nc_x_vec));
    y_dist_mat = repmat(nc_y_vec,numel(nc_y_vec),1)-repmat(nc_y_vec',1,numel(nc_y_vec));
    r_dist_mat = sqrt(x_dist_mat.^2 + y_dist_mat.^2);
    % for each spot, segment nearby nuclei and attempt to sample local
    % protein levels        
    % initialize arrays to store relevant info 
    for j = 1:numel(new_vec_fields)
        eval([new_vec_fields{j} ' = NaN(size(spot_x_vec));']);
    end
    for j = 1:numel(new_snip_fields)
        eval([new_snip_fields{j} ' = NaN(2*pt_snippet_size+1,2*pt_snippet_size+1,numel(spot_x_vec));']);
    end    
    
    % iterate through spots
    qc_mat = struct;
    for j = 1:numel(nc_x_vec)        
        % get location info
        xn = round(nc_x_vec(j));
        yn = round(nc_y_vec(j));   
        xp = round(spot_x_vec(j));
        yp = round(spot_y_vec(j));           
        zp = round(spot_z_vec(j))-1;   
        if isnan(xp)
            continue
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%% Sample protein levels %%%%%%%%%%%%%%%%%%%%%%
        % get frames
        protein_frame = protein_stack(:,:,zp);
        mcp_frame = mcp_stack(:,:,zp);
        % take locus samples
        spot_protein_vec(j) = protein_frame(yp,xp);
        spot_mcp_vec(j) = mcp_frame(yp,xp);
        
        nc_bw_final = segment_nc_neighborhood(his_sm, xn, yn, xp, ...
            yp, id_array, nb_sz, nc_index_vec(j));        
        % make sure size is reasonable and that spot is inside nucleus
        if sum(nc_bw_final(:)) < min_area || sum(nc_bw_final(:)) > max_area || ~nc_bw_final(yp,xp)   
            edge_qc_flag_vec(j) = 0;
            centroid_qc_flag_vec(j) = 0;
            continue
        end 
        
        % pull snippets     
        null_mask = nc_bw_final;
        temp_pt = protein_frame;
        temp_pt(~null_mask) = NaN;
        temp_mcp = mcp_frame;
        temp_mcp(~null_mask) = NaN;
        x_range = max(1,xp-pt_snippet_size):min(xDim,xp+pt_snippet_size);
        y_range = max(1,yp-pt_snippet_size):min(yDim,yp+pt_snippet_size);
        if numel(x_range) == 2*pt_snippet_size+1 && numel(y_range) == 2*pt_snippet_size+1
            spot_protein_snips(:,:,j) = temp_pt(y_range,x_range);
            spot_mcp_snips(:,:,j) = temp_mcp(y_range,x_range);
            spot_fov_edge_flag_vec(j) = 0;
        else
            spot_fov_edge_flag_vec(j) = 1;
        end
        % Now find locations from which to take control samples
        % Currently entertaining two alternative metrics: distance from
        % edge and distance from centroid. In principle, distance from edge
        % should likely be a better metric, however, it is subject to
        % systematic biases
              
        % Edge sampling first        
        edge_dist_mat = bwdist(~nc_bw_final);        
        spot_edge_dist = edge_dist_mat(yp,xp);        
        edge_dist_vec = edge_dist_mat(nc_bw_final);
        spot_edge_dist_vec(j) = spot_edge_dist;
        
        [edge_null_x_vec(j), edge_null_y_vec(j), edge_null_nc_vec(j), edge_qc_flag_vec(j), edge_null_mask, edge_dist_mat_nn]...
            = find_control_sample(edge_dist_vec, x_ref, y_ref, xp, yp, spot_edge_dist,...
                nc_x_vec, nc_y_vec, spot_x_vec, spot_y_vec, j, null_mask, r_dist_mat, min_sample_sep,his_sm,...
                    id_array, nb_sz, nc_index_vec, [min_area max_area]);    
               
        % Draw control samples (as appropriate)    
        if edge_qc_flag_vec(j) > 0
            xc = edge_null_x_vec(j);
            yc = edge_null_y_vec(j);            
            edge_null_protein_vec(j) = protein_frame(yc,xc);
            edge_null_mcp_vec(j) = mcp_frame(yc,xc);
            % draw snips
            temp_pt = protein_frame;
            temp_pt(~edge_null_mask) = NaN;
            temp_mcp = mcp_frame;
            temp_mcp(~edge_null_mask) = NaN;
            x_range = max(1,xc-pt_snippet_size):min(xDim,xc+pt_snippet_size);
            y_range = max(1,yc-pt_snippet_size):min(yDim,yc+pt_snippet_size);
            if numel(x_range) == 2*pt_snippet_size+1 && numel(y_range) == 2*pt_snippet_size+1
                edge_null_protein_snips(:,:,j) = temp_pt(y_range,x_range);
                edge_null_mcp_snips(:,:,j) = temp_mcp(y_range,x_range);
                edge_fov_edge_flag_vec(j) = 0;
            else
                edge_fov_edge_flag_vec(j) = 1;
            end
        end  
        
        % Now perform centroid sampling first  
        x_centroid = round(mean(x_ref(nc_bw_final)));
        y_centroid = round(mean(y_ref(nc_bw_final)));
        centroid_mat = zeros(size(x_ref));
        centroid_mat(y_centroid,x_centroid) = 1;
        centroid_dist_mat = bwdist(centroid_mat); 
%         centroid_dist_mat(~nc_bw_final) = NaN;
                       
        spot_centroid_dist = centroid_dist_mat(yp,xp);    
        centroid_dist_vec = centroid_dist_mat(nc_bw_final);
        spot_centroid_dist_vec(j) = spot_centroid_dist;
        
        [centroid_null_x_vec(j), centroid_null_y_vec(j), centroid_null_nc_vec(j), centroid_qc_flag_vec(j), centroid_null_mask, centroid_dist_mat_nn]...
            = find_control_sample(centroid_dist_vec, x_ref, y_ref, xp, yp, spot_centroid_dist,...
                nc_x_vec, nc_y_vec, spot_x_vec, spot_y_vec, j, null_mask, r_dist_mat, min_sample_sep,his_sm,...
                    id_array, nb_sz, nc_index_vec,[min_area max_area]);    
               
        % Draw control samples (as appropriate)
        xc = NaN;
        yc = NaN;
        if centroid_qc_flag_vec(j) > 0
            xc = centroid_null_x_vec(j);
            yc = centroid_null_y_vec(j);            
            centroid_null_protein_vec(j) = protein_frame(yc,xc);
            centroid_null_mcp_vec(j) = mcp_frame(yc,xc);
            % draw snips
            temp_pt = protein_frame;
            temp_pt(~centroid_null_mask) = NaN;
            temp_mcp = mcp_frame;
            temp_mcp(~centroid_null_mask) = NaN;
            x_range = max(1,xc-pt_snippet_size):min(xDim,xc+pt_snippet_size);
            y_range = max(1,yc-pt_snippet_size):min(yDim,yc+pt_snippet_size);
            if numel(x_range) == 2*pt_snippet_size+1 && numel(y_range) == 2*pt_snippet_size+1
                centroid_null_protein_snips(:,:,j) = temp_pt(y_range,x_range);
                centroid_null_mcp_snips(:,:,j) = temp_mcp(y_range,x_range);
                centroid_fov_edge_flag_vec(j) = 0;
            else
                centroid_fov_edge_flag_vec(j) = 1;
            end
        end        
        % save qc data                 
        qc_mat(j).setID = setID;
        qc_mat(j).frame = frame;
        qc_mat(j).nc_index = nc_index_vec(j);
        qc_mat(j).nc_sub_index = nc_sub_index_vec(j);
        qc_mat(j).qc_flag = edge_qc_flag_vec(j);            
        qc_mat(j).xp = xp;
        qc_mat(j).yp = yp;  
        qc_mat(j).xc_edge = edge_null_x_vec(j);
        qc_mat(j).yc_edge = edge_null_y_vec(j);
        qc_mat(j).xc_centroid = centroid_null_x_vec(j);
        qc_mat(j).yc_centroid = centroid_null_y_vec(j);
        qc_mat(j).ParticleID = particle_id_vec(j);
        sz = nb_sz;
        if edge_qc_flag_vec(j) == 2 || centroid_qc_flag_vec(j) == 2         
            sz = max([nb_sz,abs(xn - edge_null_x_vec(j)),abs(yn - edge_null_y_vec(j))...
                abs(xn - centroid_null_x_vec(j)),abs(yn - centroid_null_y_vec(j))]);
            edge_dist_mat = edge_dist_mat + edge_dist_mat_nn;
            centroid_dist_mat = centroid_dist_mat + centroid_dist_mat_nn;
        end
        y_range = max(1,yn-sz):min(yDim,yn+sz);
        x_range = max(1,xn-sz):min(xDim,xn+sz);
        qc_mat(j).x_center = median(x_range);
        qc_mat(j).y_center = median(y_range);
        qc_mat(j).mcp_snip = his_sm(y_range,x_range);
        qc_mat(j).edge_dist_snip = edge_dist_mat(y_range,x_range);        
        qc_mat(j).centroid_dist_snip = centroid_dist_mat(y_range,x_range);        
    end 
    qc_structure(i).qc_mat = qc_mat;
    
    % map data back to nucleus_struct    
    for j = 1:numel(nc_index_vec)
        nc_index = nc_index_vec(j);
        nc_sub_index = nc_sub_index_vec(j);
        for k = 1:numel(new_vec_fields)
            vec = eval(new_vec_fields{k});
            nucleus_struct(nc_index).(new_vec_fields{k})(nc_sub_index) = vec(j);
        end
        for k = 1:numel(new_snip_fields)
            snip = eval(new_snip_fields{k});
            if snip_scale_factor ~= 1
                snip = imresize(snip,snip_scale_factor);
            end
            nucleus_struct(nc_index).(new_snip_fields{k})(:,:,nc_sub_index) = snip(:,:,j);
        end
    end
    disp([num2str(i) ' of ' num2str(size(set_frame_array,1)) ' frames completed (' num2str(toc) 'seconds)'])
end

% save qc data
for i = 1:numel(qc_structure)
    qc_mat = qc_structure(i).qc_mat;
    for  j = 1:numel(qc_mat)
        qc_spot = qc_mat(j);
        if ~isfield(qc_spot,'ParticleID')
            continue
        end
        ParticleID = qc_spot.ParticleID;
        if isempty(ParticleID)
            continue
        end        
        frame = qc_spot.frame;      
        save_name = [SnipPath 'pt' num2str(1e4*ParticleID) '_frame' sprintf('%03d',frame) '.mat'];
        save(save_name,'qc_spot');
    end
end
toc
% save updated nucleus structure
nucleus_struct_protein = nucleus_struct;
save([DataPath 'nucleus_struct_protein.mat'],'nucleus_struct_protein','-v7.3') 