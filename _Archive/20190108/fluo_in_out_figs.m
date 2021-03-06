clear
close all
% add path to utilities
addpath('../utilities')
% set read and write paths
dropboxFolder = 'E:\Nick\LivemRNA\Dropbox (Personal)\';
[~, DataPath, FigRoot] =   header_function(DropboxFolder, project);
FigPath = [FigRoot '\' project '\fluo_input_output\'];
mkdir(FigPath)

% load data set
load([DataPath 'fluo_in_out.mat'],'fluo_io_struct')
%%%
VoxelSize = fluo_io_struct.voxel_size;
PixelSize = sqrt(VoxelSize/.5);

% make fluo heatmap plots

% extract snips
fluo_q1_snip = fluo_io_struct.fluo_q1_snip;
fluo_q4_snip = fluo_io_struct.fluo_q4_snip;
snip_size = size(fluo_q1_snip,1);

% determine bounds
fluo_lb = round(prctile([fluo_q1_snip(:)'  fluo_q4_snip(:)'],1)*10)/10;
fluo_ub = round(prctile([fluo_q1_snip(:)'  fluo_q4_snip(:)'],99)*10)/10;

% define tick strings
xtick_string = "set(gca,'xtick',1:5:snip_size,'xticklabel',round(((1:5:snip_size)-round(snip_size/2))*PixelSize,2))";
ytick_string = "set(gca,'ytick',1:5:snip_size,'yticklabel',round(((1:5:snip_size)-round(snip_size/2))*PixelSize,2))";

% Make titles
spot_mcp_title_dim = 'snail MCP-mCherry (bottom quintile)';
spot_mcp_title_bright = 'snail MCP-mCherry (top quintile)';
spot_mcp_Clabel = 'snail intensity (AU)';

% define heatmap
magenta = [162 60 150]/256/.65;
gray = magenta/3;
cm_magenta = interp1([0,1],vertcat(gray,magenta),linspace(0,1,128));
% dim spot
dim_fluo_fig = makeHeatmapPlots(fluo_q1_snip, 1, spot_mcp_title_dim, spot_mcp_Clabel,cm_magenta,PixelSize,fluo_lb,fluo_ub);
bright_fluo_fig = makeHeatmapPlots(fluo_q4_snip, 1, spot_mcp_title_bright, spot_mcp_Clabel,cm_magenta,PixelSize,fluo_lb,fluo_ub);
saveas(dim_fluo_fig, [FigPath 'dim_fluo_heatmap.pdf'])
saveas(bright_fluo_fig, [FigPath 'bright_fluo_heatmap.pdf'])

%%% Now protein
% extract snips
protein_q1_snip = fluo_io_struct.protein_q1_snip / VoxelSize;
protein_q4_snip = fluo_io_struct.protein_q4_snip / VoxelSize;

% determine bounds
protein_lb = 0;%prctile([protein_q1_snip(:)'  protein_q4_snip(:)'],1);
protein_ub = 90;%floor(max([protein_q1_snip(:)'  protein_q4_snip(:)'])/10)*10;

% Make titles
spot_protein_title_dim = 'Dorsal-Venus at Active snail locus (bottom quintile)';
spot_protein_title_bright = 'Dorsal-Venus at Active snail locus (top quintile)';
spot_protein_Clabel = 'Venus intensity (AU)';

NL_fav_cm = flipud(brewermap(256,'RdYlBu'));
MT_fav_cm = viridis(128);
% define heatmap
% green_top = [12 118 60]/256/.65;
% green_bottom = green_top/3;
% cm_green = interp1([0,1],vertcat(green_bottom,green_top),linspace(0,1,128));
% dim spot
dim_protein_fig = makeHeatmapPlots(protein_q1_snip, 1, spot_protein_title_dim, spot_protein_Clabel,MT_fav_cm,PixelSize,protein_lb,protein_ub);
bright_protein_fig = makeHeatmapPlots(protein_q4_snip, 1, spot_protein_title_bright, spot_protein_Clabel,MT_fav_cm,PixelSize,protein_lb,protein_ub);

saveas(dim_protein_fig, [FigPath 'dim_protein_heatmap.pdf'])
saveas(bright_protein_fig, [FigPath 'bright_protein_heatmap.pdf'])

%%% Now bar plots
protein_target_mean = fluo_io_struct.protein_target_mean;
protein_target_ste = fluo_io_struct.protein_target_ste;
protein_control_mean = fluo_io_struct.protein_control_mean;
protein_control_ste = fluo_io_struct.protein_control_ste;

% calculate y limits
ymax = 1.1*max(protein_target_mean);
ymin = 1.5*min(protein_control_mean);

target_fig = figure;
hold on
b = bar([1 3 5],protein_target_mean);
b.FaceColor = NL_fav_cm(210,:);
errorbar([1 3 5],protein_target_mean,protein_target_ste,'LineStyle','none','Color','black','LineWidth',1.5,'CapSize',20)
box on
xlabel('quintile')
ylabel('Dorsal enrichment (AU)')
set(gca,'xtick',[1 3 5],'FontSize',14);
ylim([ymin ymax])
saveas(target_fig, [FigPath 'bar_target.pdf'])

control_fig = figure;
hold on
b = bar([1 3 5],protein_control_mean);
b.FaceColor = NL_fav_cm(50,:);
errorbar([1 3 5],protein_control_mean,protein_control_ste,'LineStyle','none','Color','black','LineWidth',1.5,'CapSize',20)
box on
xlabel('quintile')
ylabel('Dorsal enrichment (AU)')
set(gca,'xtick',[1 3 5],'FontSize',14);
ylim([ymin ymax])
saveas(control_fig, [FigPath 'bar_control.pdf'])


