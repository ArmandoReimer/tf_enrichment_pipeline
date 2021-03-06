% Script to analyze protein feature variation as a function of size and
% duration of transcription bursts
clear
close all
addpath('utilities')
% define core ID variables
% project = 'Dl-Ven_snaBAC-mCh';
project = 'Dl-Ven_hbP2P-mCh';
DropboxFolder =  'E:\Nick\LivemRNA\Dropbox\';
[~, DataPath, FigureRoot] =   header_function(DropboxFolder, project);
FigPath = [Figure Root '\' project '\burst_analyses\'];
mkdir(FigPath)
% load data
load([DataPath 'hmm_input_output_results.mat'])
w = 7;
K = 3;
load([DataPath 'hmm_input_output_w' num2str(w) '_K' num2str(K) '.mat'],'hmm_input_output')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% uncontrolled burst duration figures

% burst rises rise
if strcmp(project,'Dl-Ven_hbP2P-mCh')
    burst_range = 5:10;
    burst_sigma = 3;
else
    burst_range = 2:12;
    burst_sigma = 3;
end
min_buffer_len = 5;
max_buffer_len = 30;
% extract relevant arrays
lag_dur_vec = results_struct.lag_dur_vec;
lead_dur_vec = results_struct.lead_dur_vec;
center_time_vec = results_struct.center_time_vec;
prev_burst_dist_vec = results_struct.prev_burst_dist;
lead_sz_vec = results_struct.lead_size_vec;
feature_sign_vec = results_struct.feature_sign_vec;
hmm_array = results_struct.hmm_array;
swap_hmm_array = results_struct.swap_hmm_array;
spot_array = results_struct.spot_array;
swap_array = results_struct.swap_array;
virtual_array = results_struct.virtual_array;
window_size = size(swap_array,2);
% initialize data arrays
burst_rise_dur_hmm_mean = NaN(numel(burst_range),window_size);
burst_rise_dur_swap_hmm_mean = NaN(numel(burst_range),window_size);
burst_rise_dur_virt_mean = NaN(numel(burst_range),window_size);
burst_rise_dur_spot_mean = NaN(numel(burst_range),window_size);
burst_rise_dur_swap_mean = NaN(numel(burst_range),window_size);
for i = 1:numel(burst_range)
    burst_vec = burst_range(i)-burst_sigma:burst_range(i)+burst_sigma;
    burst_ft = feature_sign_vec == 1 & ismember(lag_dur_vec,burst_vec) & ...
        lead_dur_vec>= min_buffer_len & lead_dur_vec < max_buffer_len;
    % calculate averages
    burst_rise_dur_hmm_mean(i,:) = nanmean(hmm_array(burst_ft,:));
    burst_rise_dur_swap_hmm_mean(i,:) = nanmean(swap_hmm_array(burst_ft,:));
    burst_rise_dur_spot_mean(i,:) = nanmean(spot_array(burst_ft,:));
    burst_rise_dur_swap_mean(i,:) = nanmean(swap_array(burst_ft,:));
    burst_rise_dur_virt_mean(i,:) = nanmean(virtual_array(burst_ft,:));
end
% apply mild spatio-temporal averaging
% burst_rise_dur_hmm_mean = imgaussfilt(burst_rise_dur_hmm_mean,1);
% burst_rise_dur_swap_hmm_mean = imgaussfilt(burst_rise_dur_swap_hmm_mean,1);
% burst_rise_dur_spot_mean = imgaussfilt(burst_rise_dur_spot_mean,1);
% burst_rise_dur_swap_mean = imgaussfilt(burst_rise_dur_swap_mean,1);
% burst_rise_dur_virt_mean = imgaussfilt(burst_rise_dur_virt_mean,1);

% initialize data arrays
burst_fall_dur_hmm_mean = NaN(numel(burst_range),window_size);
burst_fall_dur_swap_hmm_mean = NaN(numel(burst_range),window_size);
burst_fall_dur_spot_mean = NaN(numel(burst_range),window_size);
burst_fall_dur_swap_mean = NaN(numel(burst_range),window_size);
burst_fall_dur_virt_mean = NaN(numel(burst_range),window_size);
for i = 1:numel(burst_range)
    burst_vec = burst_range(i)-burst_sigma:burst_range(i)+burst_sigma;
    burst_ft = feature_sign_vec == -1 & ismember(lead_dur_vec,burst_vec) & lag_dur_vec >= min_buffer_len ;
    % calculate averages
    burst_fall_dur_hmm_mean(i,:) = nanmean(hmm_array(burst_ft,:));
    burst_fall_dur_swap_hmm_mean(i,:) = nanmean(swap_hmm_array(burst_ft,:));
    burst_fall_dur_spot_mean(i,:) = nanmean(spot_array(burst_ft,:));
    burst_fall_dur_swap_mean(i,:) = nanmean(swap_array(burst_ft,:));
    burst_fall_dur_virt_mean(i,:) = nanmean(virtual_array(burst_ft,:));
end

% apply mild spatio-temporal averaging
% burst_fall_dur_hmm_mean = imgaussfilt(burst_fall_dur_hmm_mean,1);
% burst_fall_dur_swap_hmm_mean = imgaussfilt(burst_fall_dur_swap_hmm_mean,1);
% burst_fall_dur_spot_mean = imgaussfilt(burst_fall_dur_spot_mean,1);
% burst_fall_dur_swap_mean = imgaussfilt(burst_fall_dur_swap_mean,1);
% burst_fall_dur_virt_mean = imgaussfilt(burst_fall_dur_virt_mean,1);
% window_vec = ((1:window_size) - ceil(window_size/2))/3;
cm_rise = brewermap(128,'Reds');
inc = floor(128/numel(burst_range));

%%%%%%%%%%%%%%%%%%%%%%%%%% RISE HEATMAPS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
burst_rise_dur_hm = figure;
burst_rise_dur_hm.Name = 'target spot burst rise hmm';
hm_cm = flipud(brewermap([],'RdYlBu'));
colormap(hm_cm)
% colormap(jet(128)/1.1);
pcolor(flipud(burst_rise_dur_spot_mean))
xlabel('offset (min)')
set(gca,'xtick',1:3:window_size,'xticklabels',-5:5)
ylabel('burst duration (min)')
set(gca,'yticklabels',fliplr(round(burst_range/3,1)))
h = colorbar;
caxis([-15 20])
ylabel(h, 'Dorsal levels (au)')
set(gca,'FontSize', 12);
saveas(burst_rise_dur_hm, [FigPath 'burst_dur_rise_hm_target.tif'])

swap_rise_dur_hm = figure;
swap_rise_dur_hm.Name = 'swap spot burst rise hmm';
colormap(hm_cm)
% colormap(jet(128)/1.1);
pcolor(flipud(burst_rise_dur_swap_mean))
xlabel('offset (min)')
set(gca,'xtick',1:3:window_size,'xticklabels',-5:5)
ylabel('burst duration (min)')
set(gca,'yticklabels',fliplr(round(burst_range/3,1)))
h = colorbar;
ylabel(h, 'Dorsal levels (au)')
set(gca,'FontSize', 12);
caxis([-15 20])
saveas(swap_rise_dur_hm, [FigPath 'burst_dur_rise_hm_swap.tif'])

virt_rise_dur_hm = figure;
virt_rise_dur_hm.Name = 'virtual spot burst rise hmm';
colormap(hm_cm)
% colormap(jet(128)/1.1);
pcolor(flipud(burst_rise_dur_virt_mean))
xlabel('offset (min)')
set(gca,'xtick',1:3:window_size,'xticklabels',-5:5)
ylabel('burst duration (min)')
set(gca,'yticklabels',fliplr(round(burst_range/3,1)))
h = colorbar;
caxis([-15 20])
ylabel(h, 'Dorsal levels (au)')
set(gca,'FontSize', 12);
saveas(virt_rise_dur_hm, [FigPath 'burst_dur_rise_hm_virtual.tif'])

hmm_rise_dur_hm = figure;
hmm_rise_dur_hm.Name = 'target spot burst rise hmm';
colormap(hm_cm)
pcolor(flipud(burst_rise_dur_hmm_mean))
xlabel('offset (min)')
set(gca,'xtick',1:3:window_size,'xticklabels',-5:5)
ylabel('burst duration (min)')
set(gca,'yticklabels',fliplr(round(burst_range/3,1)))
h = colorbar;
caxis([.3 1.2])
ylabel(h, 'sna activity (au)')
set(gca,'FontSize', 12);
saveas(hmm_rise_dur_hm, [FigPath 'burst_dur_rise_hm_hmm.tif'])

swap_hmm_rise_dur_hm = figure;
swap_hmm_rise_dur_hm.Name = 'swap spot burst rise hmm';
colormap(hm_cm)
pcolor(flipud(burst_rise_dur_swap_hmm_mean))
xlabel('offset (min)')
set(gca,'xtick',1:3:window_size,'xticklabels',-5:5)
ylabel('burst duration (min)')
set(gca,'yticklabels',fliplr(round(burst_range/3,1)))
h = colorbar;
caxis([.3 1.2])
ylabel(h, 'sna activity (au)')
set(gca,'FontSize', 12);
saveas(swap_hmm_rise_dur_hm, [FigPath 'swap_burst_dur_rise_hm_hmm.tif'])

%%%%%%%%%%%%%%%%%%%%%%%%%% FALL HEATMAPS %%%%%%%%%%%%%%%%%%%%%%%%%%
window_vec = ((1:window_size) - ceil(window_size/2))/3;
cm_fall = brewermap(128,'PuBu');
inc = floor(128/numel(burst_range));

burst_fall_dur_hm = figure;
burst_fall_dur_hm.Name = 'target spot burst fall protein';
colormap(hm_cm)
% colormap(jet(128)/1.1);
pcolor(flipud(burst_fall_dur_spot_mean))
xlabel('offset (min)')
set(gca,'xtick',1:3:window_size,'xticklabels',-5:5)
ylabel('burst duration (min)')
set(gca,'yticklabels',fliplr(round(burst_range/3,1)))
h = colorbar;
caxis([-12 12])
ylabel(h, 'Dorsal levels (au)')
set(gca,'FontSize', 12);
saveas(burst_fall_dur_hm, [FigPath 'burst_dur_fall_hm_target.tif'])

swap_fall_dur_hm = figure;
swap_fall_dur_hm.Name = 'swap spot burst fall protein';
colormap(hm_cm)
pcolor(flipud(burst_fall_dur_swap_mean))
xlabel('offset (min)')
set(gca,'xtick',1:3:window_size,'xticklabels',-5:5)
ylabel('burst duration (min)')
set(gca,'yticklabels',fliplr(round(burst_range/3,1)))
h = colorbar;
caxis([-12 12])
ylabel(h, 'Dorsal levels (au)')
set(gca,'FontSize', 12);
saveas(swap_fall_dur_hm, [FigPath 'burst_dur_fall_hm_swap.tif'])

virt_fall_dur_hm = figure;
virt_fall_dur_hm.Name = 'virtual spot burst fall protein';
colormap(hm_cm)
pcolor(flipud(burst_fall_dur_virt_mean))
xlabel('offset (min)')
set(gca,'xtick',1:3:window_size,'xticklabels',-5:5)
ylabel('burst duration (min)')
set(gca,'yticklabels',fliplr(round(burst_range/3,1)))
h = colorbar;
caxis([-12 12])
ylabel(h, 'Dorsal levels (au)')
set(gca,'FontSize', 12);
saveas(virt_fall_dur_hm, [FigPath 'burst_dur_fall_hm_virtual.tif'])

hmm_fall_dur_hm = figure;
hmm_fall_dur_hm.Name = 'target burst fall hmm';
colormap(hm_cm)
pcolor(flipud(burst_fall_dur_hmm_mean))
xlabel('offset (min)')
set(gca,'xtick',1:3:window_size,'xticklabels',-5:5)
ylabel('burst duration (min)')
set(gca,'yticklabels',fliplr(round(burst_range/3,1)))
h = colorbar;
caxis([.3 1.3])
ylabel(h, 'sna activity (au)')
set(gca,'FontSize', 12);
saveas(hmm_fall_dur_hm, [FigPath 'burst_dur_fall_hm_hmm.tif'])

swap_hmm_fall_dur_hm = figure;
swap_hmm_fall_dur_hm.Name = 'swap burst fall hmm';
colormap(hm_cm)
pcolor(flipud(burst_fall_dur_swap_hmm_mean))
xlabel('offset (min)')
set(gca,'xtick',1:3:window_size,'xticklabels',-5:5)
ylabel('burst duration (min)')
set(gca,'yticklabels',fliplr(round(burst_range/3,1)))
h = colorbar;
caxis([.3 1.3])
ylabel(h, 'sna activity (au)')
set(gca,'FontSize', 12);
saveas(swap_hmm_fall_dur_hm, [FigPath 'burst_dur_fall_hm_swap_hmm.tif'])



%%
%%%%%%%%%%%%%%%%%%%%%%%%%% RISE WATERFALLS %%%%%%%%%%%%%%%%%%%%%%%%%%
burst_rise_dur_wt = figure;
index_vec = 1:numel(burst_range);
hold on
for i = 1:numel(burst_range)
    temp = burst_rise_dur_spot_mean;
    temp(index_vec~=i,:) = NaN;
    w = waterfall(temp-nanmin(temp(:)),repmat((1:numel(burst_range))',1,window_size));
    w.FaceColor = cm_rise(1+(i-1)*inc,:);
    w.FaceAlpha = .8;
    w.EdgeColor = 'black';
end
xlabel('offset (min)')
set(gca,'xtick',1:3:window_size,'xticklabels',-5:5)    
zlabel('Dorsal levels (au)')    
view(-15,20)
grid on
saveas(burst_rise_dur_wt, [FigPath 'burst_dur_rise_waterfall_target.tif'])

hmm_rise_dur_wt = figure;
index_vec = 1:numel(burst_range);
hold on
for i = 1:numel(burst_range)
    temp = burst_rise_dur_hmm_mean;
    temp(index_vec~=i,:) = NaN;
    w = waterfall(temp-nanmin(temp(:)),repmat((1:numel(burst_range))',1,window_size));
    w.FaceColor = cm_rise(1+(i-1)*inc,:);
    w.FaceAlpha = .8;
    w.EdgeColor = 'black';
end
xlabel('offset (min)')
set(gca,'xtick',1:3:window_size,'xticklabels',-5:5)    
zlabel('sna activity (au)')    
view(-15,20)
grid on
saveas(hmm_rise_dur_wt, [FigPath 'burst_dur_rise_waterfall_hmm.tif'])




%%%%%%%%%%%%%%%%%%%%%%%%%% FALL WATERFALLS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
burst_fall_dur_wt = figure;
index_vec = 1:numel(burst_range);
hold on
for i = 1:numel(burst_range)
    temp = burst_fall_dur_spot_mean;
    temp(index_vec~=i,:) = NaN;
    w = waterfall(temp-nanmin(temp(:)),repmat((1:numel(burst_range))',1,window_size));
    w.FaceColor = cm_fall(1+(i-1)*inc,:);
    w.FaceAlpha = .8;
    w.EdgeColor = 'black';
end
xlabel('offset (min)')
set(gca,'xtick',1:3:window_size,'xticklabels',-5:5)    
zlabel('Dorsal levels (au)')    
view(15,20)
grid on
saveas(burst_fall_dur_wt, [FigPath 'burst_dur_fall_waterfall_target.tif'])


hmm_fall_dur_wt = figure;
index_vec = 1:numel(burst_range);
hold on
for i = 1:numel(burst_range)
    temp = burst_fall_dur_hmm_mean;
    temp(index_vec~=i,:) = NaN;
    w = waterfall(temp-nanmin(temp(:)),repmat((1:numel(burst_range))',1,window_size));
    w.FaceColor = cm_fall(1+(i-1)*inc,:);
    w.FaceAlpha = .8;
    w.EdgeColor = 'black';
end
xlabel('offset (min)')
set(gca,'xtick',1:3:window_size,'xticklabels',-5:5)    
zlabel('sna activity (au)')    
view(15,20)
grid on
saveas(hmm_fall_dur_wt, [FigPath 'burst_dur_fall_waterfall_hmm.tif'])

%%%% Burst example heatmap

% index = 129;
% burst_fig = figure;
% z_vec = hmm_input_output(129).z_vec > 1;
% spot_vec = hmm_input_output(129).spot_protein;
% time_vec = hmm_input_output(129).time;

% yyaxis left
% pp = plot(time_vec/60,spot_vec,'LineWidth',1.5,'Color',hm_cm(end-10,:));
% ax = gca;
% ax.YColor = 'black';
% ylabel('Dorsal concentration (au)')
% 
% yyaxis right 
% hold on
% x = [time_vec/60 time_vec/60];
% y = [z_vec' z_vec'];
% s = area(x,y,'FaceColor',[0 0 0],'EdgeAlpha',0,'FaceAlpha',.3);
% ylim([0 1.35])
% ax = gca;
% ax.YColor = 'black';
% ylabel('activity state')
% p = plot(0,0);
% xlabel('minutes')
% xlim([10 40])
% legend([pp s],'protein level','transcriptional activity');
% StandardFigure(p,gca)
% saveas(burst_fig, [figPath 'burst_example.tif'])