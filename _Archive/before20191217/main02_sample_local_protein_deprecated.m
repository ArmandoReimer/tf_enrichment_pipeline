% main02_sample_local_protein(project, RawPath)
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
% protein_channel: 
%
% OPTIONS
% DropboxFolder: Path to data folder where you wish to save
%                pipeline-generated data sets and figures. If this
%                var is not specified, output will be saved one level
%                above git repo in the folder structure 
%
% OUTPUT: nucleus_struct_protein: compiled data set with protein samples

function nucleus_struct_protein = main02_sample_local_protein(project,rawPath,proteinChannel,varargin)
tic
zeissFlag = 0;
DataPath = ['../dat/' project '/'];
for i = 1:numel(varargin)
    if strcmpi(varargin{i}, '')
        disp('accounting for zeiss offset')
        zeissFlag = 1;
    elseif strcmpi(varargin{i}, 'DropboxFolder')        
        DataPath = [varargin{i+1} '/ProcessedEnrichmentData/' project '/'];
    end
    if isstring(varargin{i})
        if ismember(varargin{i},{'dropboxFolder', 'zeissFlag'})       
            eval([varargin{i} '=varargin{i+1}']);
        end
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
mcp_channel = options(options~=proteinChannel);
% generate indexing vectors
frame_ref = [nucleus_struct.frames];
nc_x_ref = [nucleus_struct.xPos];
nc_y_ref = [nucleus_struct.yPos];
spot_x_ref = [nucleus_struct.xPosParticle]-zeissFlag;
spot_y_ref = [nucleus_struct.yPosParticle]-zeissFlag;
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
new_vec_fields = {'spot_protein_vec', 'edge_null_protein_vec','rand_null_protein_vec',...
    'spot_mcp_vec','edge_null_mcp_vec','rand_null_mcp_vec',...
    'edge_qc_flag_vec', 'rand_qc_flag_vec', 'spot_fov_edge_flag_vec','edge_fov_edge_flag_vec', ...
    'rand_fov_edge_flag_vec','edge_null_x_vec', 'edge_null_y_vec', 'edge_null_nc_vec',...
    'rand_null_x_vec', 'rand_null_y_vec', 'rand_null_nc_vec',...
    'spot_edge_dist_vec'};

new_snip_fields = {'spot_protein_snips', 'edge_null_protein_snips','rand_null_protein_snips',...
    'spot_mcp_snips','edge_null_mcp_snips','rand_null_mcp_snips'};
% get most common snip size
default_snip_size = mode(pt_snippet_size_vec);
% Initialize fields
for i = 1:numel(nucleus_struct)
    ref = nucleus_struct(i).xPos;
    for j = 1:numel(new_vec_fields)
        nucleus_struct(i).(new_vec_fields{j}) = NaN(size(ref));
    end
    for j = 1:numel(new_snip_fields)
        nucleus_struct(i).(new_snip_fields{j}) = [];%NaN(2*default_snip_size+1,2*default_snip_size+1,numel(ref));
    end
    nucleus_struct(i).snip_frame_vec = [];
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
for i = 20:size(set_frame_array,1)
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
    [mcp_stack, protein_stack] = load_stacks(rawPath, src, frame, mcp_channel);
    
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
        x_nucleus = round(nc_x_vec(j));
        y_nucleus = round(nc_y_vec(j));   
        x_spot = round(spot_x_vec(j));
        y_spot = round(spot_y_vec(j));           
        zp = round(spot_z_vec(j))-1;   
        if isnan(x_spot)
            continue
        end
                
        %%%%%%%%%%%%%%%%%%%%%% Sample protein levels %%%%%%%%%%%%%%%%%%%%%%
        % get frames
        protein_frame = protein_stack(:,:,zp);
        mcp_frame = mcp_stack(:,:,zp);
        % take locus samples
        spot_protein_vec(j) = protein_frame(y_spot,x_spot);
        spot_mcp_vec(j) = mcp_frame(y_spot,x_spot);
        
        spot_nc_mask = segment_nc_neighborhood(his_sm, x_nucleus, y_nucleus, x_spot, ...
            y_spot, id_array, nb_sz, nc_index_vec(j));        
        % make sure size is reasonable and that spot is inside nucleus
        if sum(spot_nc_mask(:)) < min_area || sum(spot_nc_mask(:)) > max_area || ~spot_nc_mask(y_spot,x_spot)   
            edge_qc_flag_vec(j) = 0;
            rand_qc_flag_vec(j) = 0;
            continue
        end 
        
        % pull snippets     
        null_mask = spot_nc_mask;
        temp_pt = protein_frame;
        temp_pt(~null_mask) = NaN;
        temp_mcp = mcp_frame;
        temp_mcp(~null_mask) = NaN;
        x_range = max(1,x_spot-pt_snippet_size):min(xDim,x_spot+pt_snippet_size);
        y_range = max(1,y_spot-pt_snippet_size):min(yDim,y_spot+pt_snippet_size);
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
        edge_dist_mat = bwdist(~spot_nc_mask);        
        spot_edge_dist = edge_dist_mat(y_spot,x_spot);        
        edge_dist_vec = edge_dist_mat(spot_nc_mask);
        spot_edge_dist_vec(j) = spot_edge_dist;
        
        [edge_null_x_vec(j), edge_null_y_vec(j), edge_null_nc_vec(j), edge_qc_flag_vec(j), edge_null_mask, edge_dist_mat_nn]...
            = find_control_sample(edge_dist_vec, x_ref, y_ref, x_spot, y_spot, spot_edge_dist,...
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
        
        % Now take a random sample           
        
         % get position vectors for nucleus mask
        x_pos_vec_spot = x_ref(spot_nc_mask);
        y_pos_vec_spot = y_ref(spot_nc_mask);
        % calculate distance from spot
        x_sep_vec = x_pos_vec_spot - x_spot;        
        y_sep_vec = y_pos_vec_spot - y_spot;
        r_sep_vec = sqrt(x_sep_vec.^2 + y_sep_vec.^2);
        sample_index_vec = 1:numel(x_sep_vec);
        % filter for regions far enough away from locus
        cr_filter = r_sep_vec >= min_sample_sep;
        sample_distances = r_sep_vec(cr_filter);
        sample_index_vec = sample_index_vec(cr_filter);
        % if candidate found, then proceed. Else look to neighboring nuclei
        if ~isempty(sample_distances)
            sample_index = randsample(sample_index_vec,1);
            rand_null_x_vec(j) = x_pos_vec_spot(sample_index);
            rand_null_y_vec(j) = y_pos_vec_spot(sample_index);
            rand_qc_flag_vec(j) = 1;               
        else
            error('unable to draw random sample. Check "PixelSize" and "min_sample_sep" variables')
        end
                
        % Draw control samples (as appropriate)
        if rand_qc_flag_vec(j) > 0
            xc = rand_null_x_vec(j);
            yc = rand_null_y_vec(j);            
            rand_null_protein_vec(j) = protein_frame(yc,xc);
            rand_null_mcp_vec(j) = mcp_frame(yc,xc);
            % draw snips
            temp_pt = protein_frame;
            temp_pt(~null_mask) = NaN;
            temp_mcp = mcp_frame;
            temp_mcp(~null_mask) = NaN;
            x_range = max(1,xc-pt_snippet_size):min(xDim,xc+pt_snippet_size);
            y_range = max(1,yc-pt_snippet_size):min(yDim,yc+pt_snippet_size);
            if numel(x_range) == 2*pt_snippet_size+1 && numel(y_range) == 2*pt_snippet_size+1
                rand_null_protein_snips(:,:,j) = temp_pt(y_range,x_range);
                rand_null_mcp_snips(:,:,j) = temp_mcp(y_range,x_range);
                rand_fov_edge_flag_vec(j) = 0;
            else
                rand_fov_edge_flag_vec(j) = 1;
            end
        end        
        % save qc data                 
        qc_mat(j).setID = setID;
        qc_mat(j).frame = frame;
        qc_mat(j).nc_index = nc_index_vec(j);
        qc_mat(j).nc_sub_index = nc_sub_index_vec(j);
        qc_mat(j).qc_flag = edge_qc_flag_vec(j);            
        qc_mat(j).xp = x_spot;
        qc_mat(j).yp = y_spot;  
        qc_mat(j).xc_edge = edge_null_x_vec(j);
        qc_mat(j).yc_edge = edge_null_y_vec(j);
        qc_mat(j).xc_rand = rand_null_x_vec(j);
        qc_mat(j).yc_rand = rand_null_y_vec(j);
        qc_mat(j).ParticleID = particle_id_vec(j);
        sz = nb_sz;
        if edge_qc_flag_vec(j) == 2         
            sz = max([nb_sz,abs(x_nucleus - edge_null_x_vec(j)),abs(y_nucleus - edge_null_y_vec(j))...
                abs(x_nucleus - rand_null_x_vec(j)),abs(y_nucleus - rand_null_y_vec(j))]);
            edge_dist_mat = edge_dist_mat + edge_dist_mat_nn;
        end
        y_range = max(1,y_nucleus-sz):min(yDim,y_nucleus+sz);
        x_range = max(1,x_nucleus-sz):min(xDim,x_nucleus+sz);
        qc_mat(j).x_center = median(x_range);
        qc_mat(j).y_center = median(y_range);
        qc_mat(j).mcp_snip = his_sm(y_range,x_range);
        qc_mat(j).edge_dist_snip = edge_dist_mat(y_range,x_range);        
        qc_mat(j).rand_dist_snip = edge_dist_mat(y_range,x_range);        
    end 
    qc_structure(i).qc_mat = qc_mat;
    
    % map data back to nucleus_struct    
    for j = 1:numel(nc_index_vec)
        nc_index = nc_index_vec(j);
        nc_sub_index = nc_sub_index_vec(j);
        frame = nucleus_struct(nc_index).frames(nc_sub_index);
        for k = 1:numel(new_vec_fields)
            vec = eval(new_vec_fields{k});
            nucleus_struct(nc_index).(new_vec_fields{k})(nc_sub_index) = vec(j);
        end
        if rand_qc_flag_vec(j) > 0 || edge_qc_flag_vec(j) > 0
            for k = 1:numel(new_snip_fields)
                snip = eval([new_snip_fields{k} '(:,:,j)']);
                ind = numel(nucleus_struct(nc_index).snip_frame_vec)+1;
                if snip_scale_factor ~= 1
                    snip = imresize(snip,snip_scale_factor);
                end
                nucleus_struct(nc_index).(new_snip_fields{k})(:,:,ind) = snip;
            end
            nucleus_struct(nc_index).snip_frame_vec(ind) = frame;
        end
    end
    disp([num2str(i) ' of ' num2str(size(set_frame_array,1)) ' frames completed'])
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