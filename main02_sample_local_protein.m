% main02_sample_local_protein(project, rawPath, proteinChannel, varargin)
%
% DESCRIPTION
% Script to generate samples of local protein contentration at active gene
% loci and selected control locations
%
% ARGUMENTS
% project: master ID variable
%
% rawPath: Full or relative path to PreProcessed folder (pr
%               equivalent)
%
% proteinChannel: Integer corresponding to protein channel
%
% OPTIONS
% dropboxFolder: Path to data folder where you wish to save
%                pipeline-generated data sets and figures. If this
%                var is not specified, output will be saved one level
%                above git repo in the folder structure
% zeissFlag: Integer specifying a pixel correction (1 pixel needed for
%               data taken on Zeiss780)
% segmentNuclei: choose whether to segment or re-segment nuclei. follow
% with true or false (default)
%
% OUTPUT: nucleus_struct_protein: compiled data set with protein samples

function nucleus_struct_protein = main02_sample_local_protein(project,DropboxFolder,varargin)

warning('off', 'MATLAB:MKDIR:DirectoryExists');


addpath('./utilities')
ROIRadiusSpot = .2; % radius (um) of region used to query and compare TF concentrations
minSampleSepUm = 1.5; %um
mf_samp_rad = 0.8; % distance (um) from nucleus center to include in sample
minEdgeSepUm = .5; %um
shouldSegmentNuclei = 0;
maskingMethod = 'gradientOtsuHulls';
% PSF info for 3D sampling
use_psf_fit_dims = false; % if true, use fits from PSF fitting
xy_sigma_um = 0.25;% um
z_sigma_um = 0.6; % um
ignoreQC = false;
min_rad_um = 2; % set min and max acceptable area for nucleus segmentation
max_rad_um = 4; %this needs to be 6um for nc12. 4um for nc14
sm_kernel_um = 1; % size of gaussian smoothing kernel
nb_size_um = 10; % determine size of neighborhood to use
pt_snippet_size_um = 1.5; % set snippet to be 3um in size
display_figures = false;
askToOverwrite = true;

% rawPath = 'E:\LocalEnrichment\Data\PreProcessedData\';
pth = getDorsalFolders;
rawPath = [pth, filesep, 'PreProcessedData\'];

proteinChannel = 1;


[~, DataPath, ~] =   header_function(DropboxFolder, project);

for i = 1:(numel(varargin)-1)
    if i ~= numel(varargin)
        if ~ischar(varargin{i+1})
            eval([varargin{i} '=varargin{i+1};']);
        end
    end
end


if display_figures
    
    figure('Units', 'normalized', 'Position', [0.6441 0.5922 0.3175 0.3033]);
    tiledlayout('flow', 'TileSpacing', 'none', 'Padding', 'none')
    clrmp = 'hsv';
    colormap(clrmp);
%     nPlots = 3;
%     ax = {};
%     for i = 1:nPlots
%         ax{i} = nexttile;
%     end
%     drawnow;
    
end


% Load trace data
load([DataPath '/nucleus_struct.mat'],'nucleus_struct')
load([DataPath '/psf_dims.mat'],'psf_dims')
load([DataPath '/set_key.mat'],'set_key')
snipPath = [DataPath '/qc_images/'];
refPath = [DataPath '/refFrames/'];
mkdir(refPath)
mkdir(snipPath)

% get MCP channel
% options = [1 2];
% mcp_channel = options(options~=proteinChannel);
nSets = size(set_key, 1);
for set = 1:nSets
    source = cell2mat(set_key{set, 2});
    [~,~,~,~, ~, ~, ~, ~,Channel1,Channel2,~,Channel3, ~, ~, ~] = readMovieDatabase(source);
    Channels = {Channel1{1},Channel2{1}, Channel3{1}};
    protein_channels(set) = find(contains(Channels, 'input', 'IgnoreCase', true));
    mcp_channels(set)= find(contains(Channels, 'spot', 'IgnoreCase', true)| contains(Channels, 'mcp', 'IgnoreCase', true) | contains(Channels, 'pcp', 'IgnoreCase', true)) ;
    %
    % movieFiles{set} = [rawPath, src, filesep, src, '_movieMat.mat'];
    %     if exist(movieFiles{set}, 'file')
    %         movieMat = loadMovieMat(movieFiles{set});
    %         mcpMovie{set} = movieMat(:, :, :, :, mcp_channels{set});
    %         proteinMovie{set} = movieMat(:, :, :, :, protein_channels{set});
    %         clear movieMat;
    %     end
    % end
end





write_snip_flag = false;
%%%%%%%%%%%%%%%%%%%%%%%%%
% remove all frames that do not contain a segmented particle or that

threeD_flag = nucleus_struct(1).threeD_flag;
% threeD_flag = false;
% remove entries with no particle
nucleus_struct = nucleus_struct(~isnan([nucleus_struct.ParticleID]));

fnames = fieldnames(nucleus_struct);
for i = 1:numel(nucleus_struct)
    nucleus_struct(i).time_orig = nucleus_struct(i).time;
    fluo = nucleus_struct(i).fluo;
    nan_ft = ~isnan(fluo);
    for spotIndex = 1:numel(fnames)
        vec = nucleus_struct(i).(fnames{spotIndex});
        if numel(vec) == numel(fluo)
            nucleus_struct(i).(fnames{spotIndex}) = vec(nan_ft);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%
% generate indexing vectors
frameReference = [nucleus_struct.frames];
nucleusXReference = [nucleus_struct.xPos];
nucleusYReference = [nucleus_struct.yPos];

spotXReference = [nucleus_struct.xPosParticle];
spotYReference = [nucleus_struct.yPosParticle];
spotZReference = [nucleus_struct.zPosParticle];
if threeD_flag
    spotXReference3D = [nucleus_struct.xPosParticle3D];
    spotYReference3D = [nucleus_struct.yPosParticle3D];
    spotZReference3D = [nucleus_struct.zPosParticle3D];
else
    spotXReference3D = [nucleus_struct.xPosParticle];
    spotYReference3D = [nucleus_struct.yPosParticle];
    spotZReference3D = [nucleus_struct.zPosParticle];
end
% spot_z_ref = spot_z_ref3D;

setReference = [];
masterIndexReference = [];
linearIndexReference = [];
subIndexReference = [];
proteinReference = [];
proteinQualityControlReference = [];
for i = 1:numel(nucleus_struct)
    ParticleID = nucleus_struct(i).ParticleID;
    proteinReference = [proteinReference repelem(ParticleID, numel(nucleus_struct(i).frames))];
    proteinQualityControlReference = [proteinQualityControlReference repelem(nucleus_struct(i).qc_flag, numel(nucleus_struct(i).frames))];
    setReference = [setReference repelem(nucleus_struct(i).setID,numel(nucleus_struct(i).frames))];
    masterIndexReference = [masterIndexReference repelem(round(nucleus_struct(i).ncID*1e5),numel(nucleus_struct(i).frames))];
    linearIndexReference = [linearIndexReference repelem(i,numel(nucleus_struct(i).frames))];
    subIndexReference = [subIndexReference 1:numel(nucleus_struct(i).frames)];
end

if ignoreQC
    proteinQualityControlReference = true(size(proteinQualityControlReference));
end
%%%%%%%%%%%%%%%%%%%%
%%% Set size parameters
set_vec = [nucleus_struct.setID];
set_index = unique(set_vec);

PixelSize_um = nucleus_struct(1).PixelSize;
zStep_um = nucleus_struct(1).zStep;

%convert size parameters from um to pixels
neighborhoodSize_px = round(nb_size_um./ PixelSize_um);
smoothingKernel_px = round(sm_kernel_um ./ PixelSize_um);
minArea_px = round(pi*(min_rad_um ./ PixelSize_um).^2);
maxArea_px = round(pi*(max_rad_um ./ PixelSize_um).^2);
proteinSnippetSize_px = round(pt_snippet_size_um ./ PixelSize_um);
% set min separation between control and locus to 2um
minSampleSeparation_px = round(minSampleSepUm ./ PixelSize_um);
minEdgeSeparation_px = round(minEdgeSepUm ./ PixelSize_um);

% calculate ROI size in pixels for spot and control
roiRadiusSpot_px = round(ROIRadiusSpot ./ PixelSize_um);

% calculate average frame-over-frame particle drift from data
linearDifferenceVector = diff(linearIndexReference);
xDifferenceVector = diff(spotXReference);
yDifferenceVector = diff(spotYReference);
displacementVector = sqrt(xDifferenceVector.^2+yDifferenceVector.^2);
displacementVector = displacementVector(linearDifferenceVector==0);
% sets sigma of movement for virtual spot
driftTolerance = nanmedian(displacementVector)*PixelSize_um;

% set dims for 3D protein sampling
if use_psf_fit_dims
    xy_sigma_px = psf_dims.xy_sigma;
    z_sigma_px = psf_dims.z_sigma;
else
    xy_sigma_px = round(xy_sigma_um/PixelSize_um,1);% pixels
    z_sigma_px = round(z_sigma_um/zStep_um,1); % pixels
end

%%%%%%%%%%%%%%%%%%%%%%%
% Designate fields to be added to nucleus structure
newVectorFields = {'spot_protein_vec_3d','spot_protein_vec', 'serial_null_protein_vec',...
    'serial_null_protein_vec_3d','edge_null_protein_vec','edge_null_protein_vec_3d','mf_null_protein_vec',...
    'spot_mcp_vec','edge_mcp_protein_vec','serial_qc_flag_vec','edge_qc_flag_vec', ...
    'edge_null_x_vec', 'serial_null_x_vec','serial_null_y_vec','edge_null_y_vec', 'edge_null_nc_vec',...
    'spot_edge_dist_vec','serial_null_edge_dist_vec'};


newSnippetFields = {'spot_protein_snips', 'edge_null_protein_snips',...
    'spot_mcp_snips','edge_null_mcp_snips'};

% Initialize fields
for i = 1:numel(nucleus_struct)
    reference = nucleus_struct(i).xPos;
    for spotIndex = 1:numel(newVectorFields)
        nucleus_struct(i).(newVectorFields{spotIndex}) = NaN(size(reference));
    end
end

%%%%%%%%%%%%%%%%%%%%%%
%%% make source key
sourceCell = {1,numel(set_index)};
% pt_id_vec = [nucleus_struct.ParticleID];
for s = 1:numel(set_index)
    ind = find([nucleus_struct.setID]==set_index(s),1);
    source = nucleus_struct(ind).source_path;
    sourceCell{s} = source;
end

%%% Generate reference array for set-frame combos
setFrameArray = unique([setReference' frameReference'],'row');
qualityControlStructure = struct;

%%%%%%%%%%%%%%%%%%%%%%%%%
% Nucleus segmentation
%%%%%%%%%%%%%%%%%%%%%%%%%
xDim = nucleus_struct(1).xDim;
yDim = nucleus_struct(1).yDim;
zDim = nucleus_struct(1).zDim;
[xReferenceMesh,yReferenceMesh,zReferenceMesh] = meshgrid(1:xDim,1:yDim,1:zDim);
% first check to see if segmentation files exist
segment_indices = 1:size(setFrameArray,1);
spotFrameVector = false(1,size(setFrameArray,1));
nucleusFrameVector = false(1,size(setFrameArray,1));
for i = 1:size(setFrameArray,1)
    setID_temp = setFrameArray(i,1);
    frame_temp = setFrameArray(i,2);
    nucleusReferenceFile = [refPath 'nc_ref_frame_set' sprintf('%02d',setID_temp) '_frame' sprintf('%03d',frame_temp) '.mat'];
    nucleusFrameVector(i) = isfile(nucleusReferenceFile);
    spotReferenceFile = [refPath 'spot_roi_frame_set' sprintf('%02d',setID_temp) '_frame' sprintf('%03d',frame_temp) '.mat'];
    spotFrameVector(i) = isfile(spotReferenceFile);
end
if all(spotFrameVector) && shouldSegmentNuclei
    warning('previous segmentation results found')
    y = 1;
    n = 0;
    if askToOverwrite
        overwrite = input('overwrite segmentation results? (y/n)');
    else
        overwrite = y;
    end
    shouldSegmentNuclei = overwrite;
elseif ~all(nucleusFrameVector) && ~shouldSegmentNuclei
    warning('some or all frames missing nucleus segmentation data. Segmenting missing frames only')
    shouldSegmentNuclei = 1;
    segment_indices = find(~nucleusFrameVector);
end

if shouldSegmentNuclei
    
    nSets = size(set_key, 1);
    nuclear_mov = cell(1, nSets);
    
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
    
    segmentNuclei_main02(yDim, xDim, segment_indices, setFrameArray, setReference, frameReference, nucleusXReference, ...
        nucleusYReference, masterIndexReference, spotXReference, spotYReference, set_key, rawPath, proteinChannel, ...
        smoothingKernel_px, neighborhoodSize_px, refPath,...
        display_figures, minArea_px, maxArea_px, DropboxFolder,...
        nuclear_mov)
    
    clear nuclear_mov;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%
% Protein sampling
%%%%%%%%%%%%%%%%%%%%%%%%%
disp('taking protein samples...')
for embryoIndex = 1:size(setFrameArray,1)
    
    tic
    setID = setFrameArray(embryoIndex,1);
    frame = setFrameArray(embryoIndex,2);
    
    
    % load spot and nucleus reference frames
    nucleusReferenceFile = [refPath 'nc_ref_frame_set' sprintf('%02d',setID) '_frame' sprintf('%03d',frame) '.mat'];
    load(nucleusReferenceFile,'nc_ref_frame');
    if display_figures, imagescUpdate(nexttile, nc_ref_frame, []); drawnow; end
    
  
    nucleusDistanceFrame = bwdist(~nc_ref_frame);
    if display_figures, imagescUpdate(nexttile, nucleusDistanceFrame, []); drawnow; end

    %generated by segmentNuclei_main02.
    %this is the distance of each pixel from
    %active locus
    spotReferenceFile = [refPath 'spot_roi_frame_set' sprintf('%02d',setID) '_frame' sprintf('%03d',frame) '.mat'];
    load(spotReferenceFile,'spot_dist_frame');
    if display_figures, imagescUpdate(nexttile, spot_dist_frame, []); drawnow; end
    
    
    spotRoiFrame = bwlabel(spot_dist_frame <= roiRadiusSpot_px);
    if display_figures, imagescUpdate(nexttile, spotRoiFrame, []); drawnow; end
    
    
    % get nucleus
    frameSetFilter = setReference==setID&frameReference==frame;
    nucleusXVector = nucleusXReference(frameSetFilter);
    nucleusYVector = nucleusYReference(frameSetFilter);
    
    
    % indexing vectors
    nc_sub_index_vec = subIndexReference(frameSetFilter);
    nc_lin_index_vec = linearIndexReference(frameSetFilter);
    nucleusMasterVector = masterIndexReference(frameSetFilter);
    
    % particle positions
    proteinQualityControlVector = proteinQualityControlReference(frameSetFilter);
    spot_x_vec = spotXReference(frameSetFilter);
    spot_y_vec = spotYReference(frameSetFilter);
    spot_z_vec = spotZReference(frameSetFilter);
    spot_x_vec3D = spotXReference3D(frameSetFilter);
    spot_y_vec3D = spotYReference3D(frameSetFilter);
    spot_z_vec3D = spotZReference3D(frameSetFilter);
    particle_id_vec = proteinReference(frameSetFilter);
    source = set_key(set_key.setID==setID,:).prefix{1};
    
    
    %determine protein channel
    mcp_channel = mcp_channels(setID);
    proteinChannel = protein_channels(setID);
    % load stacks
    
    mcp_stack = load_stacks(rawPath, source, frame, mcp_channel, xDim, yDim, zDim);
    protein_stack = load_stacks(rawPath, source, frame, proteinChannel, xDim, yDim, zDim);
    
   
    % initialize arrays to store relevant info
    for spotIndex = 1:numel(newVectorFields)
        eval([newVectorFields{spotIndex} ' = NaN(size(spot_x_vec));']);
    end
    for spotIndex = 1:numel(newSnippetFields)
        eval([newSnippetFields{spotIndex} ' = NaN(2*proteinSnippetSize_px+1,2*proteinSnippetSize_px+1,numel(spot_x_vec));']);
    end
   
    
[qc_mat, spot_edge_dist_vec] = main02subfunction(nucleusXVector, nucleusYVector, spot_x_vec, spot_y_vec, spot_z_vec,...
    spot_x_vec3D, spot_y_vec3D, spot_z_vec3D, proteinQualityControlVector, nc_ref_frame, nucleusMasterVector,...
    display_figures, protein_stack, mcp_stack,spotRoiFrame,...
    xReferenceMesh,yReferenceMesh,zReferenceMesh,xy_sigma_px,z_sigma_px,...
    proteinSnippetSize_px,nucleusDistanceFrame, minSampleSeparation_px,...
    roiRadiusSpot_px, nucleus_struct, nc_lin_index_vec,  nc_sub_index_vec, minEdgeSeparation_px,...
    driftTolerance,PixelSize_um, setID, frame, particle_id_vec, neighborhoodSize_px, yDim, xDim, minArea_px, maxArea_px,...
    mf_samp_rad, spot_dist_frame);

       
    qualityControlStructure(embryoIndex).qc_mat = fliplr(qc_mat);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % save snip data
    nNuclei = numel(nucleusMasterVector);
    % initialize struct to store snip data
    snip_data = struct;
    % map data back to nucleus_struct
    for spotIndex = 1:nNuclei
        nc_index = nc_lin_index_vec(spotIndex);
        nc_sub_index = nc_sub_index_vec(spotIndex);
        frame = nucleus_struct(nc_index).frames(nc_sub_index);
        for k = 1:numel(newVectorFields)
            vec = eval(newVectorFields{k});
            nucleus_struct(nc_index).(newVectorFields{k})(nc_sub_index) = vec(spotIndex);
        end
        % store snips
        for k = 1:numel(newSnippetFields)
            snip = eval([newSnippetFields{k} '(:,:,spotIndex)']);
            snip_data.(newSnippetFields{k})(:,:,spotIndex) = snip;
        end
    end
    
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
    
    % report time
    t = round(toc);
    disp([num2str(embryoIndex) ' of ' num2str(size(setFrameArray,1)) ' frames completed (' num2str(t) ' sec)'])
end



disp('saving quality control frames...')
% save quality control data
tic
particle_index = unique([nucleus_struct.ParticleID]);
particle_index = particle_index(~isnan(particle_index));
qc_particles = randsample(particle_index,min([100,numel(particle_index)]),false);
particle_index_full = [];
particle_frames_full = [];
for i = 1:numel(qualityControlStructure)
    qc_mat = qualityControlStructure(i).qc_mat;
    for  spotIndex = 1:numel(qc_mat)
        qc_spot = qc_mat(spotIndex);
        if ~isfield(qc_spot,'ParticleID')
            continue
        end
        ParticleID = qc_spot.ParticleID;
        if isempty(ParticleID) || ~ismember(ParticleID,qc_particles)
            continue
        end
        frame = qc_spot.frame;
        particle_index_full = [particle_index_full ParticleID];
        particle_frames_full = [particle_frames_full frame];
        save_name = [snipPath 'pt' num2str(1e4*ParticleID) '_frame' sprintf('%03d',frame) '.mat'];
        save(save_name,'qc_spot');
    end
end

[particle_index_full, si] = sort(particle_index_full);
particle_frames_full = particle_frames_full(si);

qc_ref_struct.particle_frames_full = particle_frames_full;
qc_ref_struct.particle_index_full = particle_index_full;
toc


% save updated nucleus structure
disp('saving nucleus structure...')
nucleus_struct_protein = nucleus_struct;
save([DataPath 'qc_ref_struct.mat'],'qc_ref_struct')
save([DataPath 'nucleus_struct_protein.mat'],'nucleus_struct_protein','-v7.3')