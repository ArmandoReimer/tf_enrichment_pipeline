% main04_check_control_selection(project)
%
% DESCRIPTION
% Generates figures overlaying segmentation results onto histone channel
% and interface in which user can accept or reject spot assignments
%
% ARGUMENTS
% project: master ID variable

function main03_check_control_selection(project,varargin)
close all
% specify paths
DataPath = ['../dat/' project '/'];

for i = 1:numel(varargin)    
    if strcmpi(varargin{i}, 'DropboxFolder')        
        DataPath = [varargin{i+1} '/ProcessedEnrichmentData/' project '/'];
    end
end

SnipPath = [DataPath 'qc_images/'];



% load data
load([DataPath '/nucleus_struct_protein.mat']);
% snip_files = dir([SnipPath '*.mat']);
% check to see if nucleus structure already contains qc review info
if isfield(nucleus_struct_protein, 'qc_review_vec')
    qc_review_vec = [nucleus_struct_protein.qc_review_vec];
else 
    qc_review_vec = NaN(size([nucleus_struct_protein.xPos]));
end
edge_qc_flag_vec = [nucleus_struct_protein.edge_qc_flag_vec];
% set start frame
all_frames = find(edge_qc_flag_vec~=0&~isnan(edge_qc_flag_vec));
outstanding_frames = find(isnan(qc_review_vec)&(edge_qc_flag_vec~=0&~isnan(edge_qc_flag_vec)));
% generate indexing vectors
frame_index = [nucleus_struct_protein.frames];
set_index = [];
particle_index = [];
for i = 1:numel(nucleus_struct_protein)
    set_index = [set_index repelem(nucleus_struct_protein(i).setID, numel(nucleus_struct_protein(i).frames))];
    particle_index = [particle_index repelem(nucleus_struct_protein(i).ParticleID, numel(nucleus_struct_protein(i).frames))];
end    

% iterate through snip files
exit_flag = 0;
cm = jet(64);

index = find(all_frames==outstanding_frames(1));
while ~exit_flag
    % create sister_struct(i) struct    
    frame = frame_index(all_frames(index));
    ParticleID = particle_index(all_frames(index));
    setID = set_index(all_frames(index));    
    % load snip data
    load([SnipPath 'pt' num2str(1e4*ParticleID) '_frame' sprintf('%03d',frame) '.mat']);
    
    %%%%%%%%%%%%%%%%%%%%%%% load image stack %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    cc = '';
    while ~strcmp(cc,'0')&&~strcmp(cc,'1')      
        edge_dist_snip = qc_spot.edge_dist_snip;
        edge_dist_rescaled = 1 + ceil(edge_dist_snip/max(edge_dist_snip(:)) * 63);
        cm = [[1 1 1] ; cm(2:end,:)];
        edge_dist_rgb = ind2rgb(edge_dist_rescaled,cm);
        
        centroid_dist_snip = qc_spot.centroid_dist_snip;
        centroid_dist_snip(edge_dist_rescaled==1) = 0;
        centroid_dist_rescaled = 1 + ceil(centroid_dist_snip/max(centroid_dist_snip(:)) * 63);               
        centroid_dist_rgb = ind2rgb(centroid_dist_rescaled,cm);
        % get frame center
        x_center = qc_spot.x_center;
        y_center = qc_spot.y_center;
        yDim = ceil(size(edge_dist_snip,1)/2);
        xDim = ceil(size(edge_dist_snip,2)/2);
        
        qc_fig = figure('Position',[0 0 1024 512]);                 
        subplot(1,2,1)
        imshow(imadjust(mat2gray(qc_spot.mcp_snip)),'InitialMagnification','fit');                        
        hold on
        p = imshow(edge_dist_rgb);        
        p.AlphaData = .4;        
        s1 = scatter(qc_spot.xp-x_center+xDim,qc_spot.yp-y_center+yDim,30,'MarkerFaceColor',cm(30,:),'MarkerEdgeAlpha',0);
        s2 = scatter(qc_spot.xc_edge-x_center+xDim,qc_spot.yc_edge-y_center+yDim,30,'MarkerFaceColor',cm(60,:),'MarkerEdgeAlpha',0);
        legend([s1 s2], 'spot', 'control')          
        title('Edge Distance Sample')
        
        subplot(1,2,2)
        imshow(imadjust(mat2gray(qc_spot.mcp_snip)),'InitialMagnification','fit');                        
        hold on
        p = imshow(centroid_dist_rgb);        
        p.AlphaData = .4;        
        s1 = scatter(qc_spot.xp-x_center+xDim,qc_spot.yp-y_center+yDim,30,'MarkerFaceColor',cm(30,:),'MarkerEdgeAlpha',0);
        s2 = scatter(qc_spot.xc_centroid-x_center+xDim,qc_spot.yc_centroid-y_center+yDim,30,'MarkerFaceColor',cm(60,:),'MarkerEdgeAlpha',0);
        legend([s1 s2], 'spot', 'control')  
        title('Centroid Distance Sample')
        
        set(gcf,'Name',['Particle ' num2str(ParticleID) ' Frame ' num2str(frame) ' (' num2str(index) ' of ' num2str(numel(all_frames)) ')'])
        if qc_review_vec(all_frames(index)) == 1
            set(gcf, 'color', 'green')
        elseif qc_review_vec(all_frames(index)) == 0
            set(gcf, 'color', 'red')
        end
        ct=waitforbuttonpress;
        cc=get(qc_fig,'currentcharacter');
        if strcmp(cc,'1')||strcmp(cc,'0')                       
            nucleus_struct(qc_spot.nc_index).qc_review_vec(qc_spot.nc_sub_index) = eval(cc);
            index = min(numel(all_frames),index + 1);
        elseif strcmp(cc,'x')
            exit_flag = 1;
            break        
        elseif strcmp(cc,'n')
            index = max(1,index-1);
            break
        elseif strcmp(cc,'m')
            index = min(numel(all_frames),index+1);
            break
        elseif strcmp(cc,'j')
            index = input('enter desired index: ');
            break
        end       
    end 
    close all
    if exit_flag
        disp('Exiting')
        break
    end
end
% nucleus_struct_protein.qc_review_vec = qc_review_vec;
save([DataPath 'nucleus_struct_protein.mat'],'nucleus_struct_protein','-v7.3')