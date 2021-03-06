% script to delve into details of spatio-temporal enrichment dynamics
clear
close all
addpath('utilities')
% define core ID variables
project = 'Dl-Ven_snaBAC-mCh';
% project = 'Dl-Ven_hbP2P-mCh';
dropboxFolder =  'E:\Nick\LivemRNA\Dropbox\';
dataPath = [dropboxFolder 'ProcessedEnrichmentData\' project '\'];
figPath = [dropboxFolder 'LocalEnrichmentFigures\' project '\enrichment_disaggregation\'];
mkdir(figPath)
% load data
load([dataPath 'nucleus_struct_protein.mat'])
%%
tic
%%% Examine distribution of enrichment event sizes
DistLim = .8;
PixelSize = nucleus_struct_protein(1).PixelSize;
spot_protein_vec = [nucleus_struct_protein.spot_protein_vec];
null_protein_vec = [nucleus_struct_protein.mf_null_protein_vec];
enrichment_vec = (spot_protein_vec - null_protein_vec);
time_vec = [nucleus_struct_protein.time];
dist_vec = [nucleus_struct_protein.spot_edge_dist_vec]*PixelSize;
dist_filter = dist_vec >= DistLim;

% generate simple histogram for enrichment values
aggregate_hist = figure;
h = histogram(spot_protein_vec(dist_filter)-null_protein_vec(dist_filter),...
        'Normalization','probability','EdgeAlpha',0);
grid on
xlabel('Dl enrichment at{\it sna} locus (au)')
ylabel('share')
xlim([-200 600])
saveas(aggregate_hist,[figPath 'aggregate_event_hist.tif'])

% break down by time and by protein concentration
t_window = 5;
time_bins = (10:2:60)-5;
mf_bins = linspace(prctile(null_protein_vec,5),prctile(null_protein_vec,95),numel(time_bins));
mf_window = median(diff(mf_bins)) * 5;

% initialize arrays
e_bins = linspace(-200,600,numel(time_bins));
mf_dist_array = zeros(numel(time_bins), numel(e_bins)-1);
time_dist_array = zeros(numel(time_bins), numel(e_bins)-1);
mf_vs_time_array = NaN(numel(time_bins),numel(mf_bins));
ct_array = zeros(size(mf_vs_time_array));
for t = 1:numel(time_bins)
    % generate time counts
    t_center = time_bins(t);
    t_ft = time_vec/60 >= t_center - t_window & time_vec/60 < t_center + t_window;
    t_ct = histcounts(enrichment_vec(dist_filter & t_ft),e_bins,'Normalization','probability');
    time_dist_array(t,:) = t_ct;
    % generate mf counts
    mf_center = mf_bins(t);
    mf_ft = null_protein_vec >= mf_center - mf_window & null_protein_vec < mf_center + mf_window;
    mf_ct = histcounts(enrichment_vec(dist_filter & mf_ft),e_bins,'Normalization','probability');
    mf_dist_array(t,:) = mf_ct;    
    for m = 1:numel(mf_bins)
        mf_center = mf_bins(m);
        mf_ft = null_protein_vec >= mf_center - mf_window & null_protein_vec < mf_center + mf_window;
        % get average for this mf_v_time 
        mf_vs_time_array(t,m) = nanmean(enrichment_vec(dist_filter & mf_ft & t_ft));
        ct_array(t,m) = sum(dist_filter & mf_ft & t_ft);
    end
end


% make figures
% time heatmap 
time_dist_hm = figure;
hm_cm1 = flipud(brewermap([],'RdYlBu'));
colormap(hm_cm1)
imagesc(time_dist_array);
xlabel('Dl enrichment at{\it sna} locus (au)')
set(gca,'xtick',1:5:26,'xticklabels',round(e_bins(1:5:26)/5)*5)
ylabel('time (min)')
set(gca,'ytick',(1:5:26),'yticklabels',time_bins((1:5:26)))
h = colorbar;
ylabel(h, 'share')
set(gca,'FontSize', 12);
caxis([0 .175])
saveas(time_dist_hm,[figPath 'time_event_heatmap.tif'])

mf_dist_hm = figure;
hm_cm1 = flipud(brewermap([],'RdYlBu'));
colormap(hm_cm1)
imagesc(mf_dist_array);
xlabel('Dl enrichment at{\it sna} locus (au)')
set(gca,'xtick',1:5:26,'xticklabels',round(e_bins(1:5:26)/5)*5)
ylabel('background Dl concentration (au)')
set(gca,'ytick',1:5:26,'yticklabels',round(mf_bins(1:5:26)/5)*5)
h = colorbar;
ylabel(h, 'share')
set(gca,'FontSize', 12);
caxis([0 .175])
saveas(mf_dist_hm,[figPath 'mf_event_heatmap.tif'])

pt_array = mf_vs_time_array;

mf_vs_time_hm = figure;
hm_cm2 = flipud(brewermap([],'RdBu'));
colormap(hm_cm1)
p = imagesc(pt_array');
p.AlphaData = .2+(ct_array'>25);
xlabel('time (min)')
ylabel('background Dl concentration (au)')
set(gca,'ytick',1:5:26,'yticklabels',round(mf_bins(1:5:26)/5)*5)
set(gca,'xtick',(1:5:26),'xticklabels',time_bins((1:5:26)))
h = colorbar;
ylabel(h, 'Dl enrichment at{\it sna} locus')
set(gca,'FontSize', 12);
caxis([60 180])
saveas(mf_vs_time_hm,[figPath 'mf_vs_time_enrichment_heatmap.tif'])
toc