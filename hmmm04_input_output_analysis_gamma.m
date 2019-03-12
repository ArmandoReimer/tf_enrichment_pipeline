% Script to probe relationship between input protein concentration and
% output transcriptional response
clear 
close all
% define ID variables
K = 3;
w = 7;
project = 'Dl_Venus_snaBAC_mCherry_Leica_hp';
dropboxFolder =  'E:\Nick\Dropbox (Garcia Lab)\';
dataPath = [dropboxFolder '\ProcessedEnrichmentData\' project '\'];
figPath = [dropboxFolder '\LocalEnrichmentFigures\' project '\hmm_input_output_K' num2str(K) '_w' num2str(w) '\'];
mkdir(figPath)

nTraces = 50; % number of individual traces to select for plotting
window_size = 10; % number of lags over which to track protein/fluo dynamics
nBoots = 100;
n_ref_hist_bins = 10;
min_time = 10;
make_trace_plots = 0;
% extract protein, gene, fluorophore info
underscores = strfind(project,'_');
protein_name = project(1:underscores(1)-1);
protein_fluor = project(underscores(1)+1:underscores(2)-1);
gene_name = project(underscores(2)+1:underscores(3)-1);
if numel(underscores) == 3
    ind = numel(project);
else
    ind = underscores(4)-1;
end
gene_fluor = project(underscores(3)+1:end);


% load data set
load([dataPath 'hmm_input_output_w' num2str(w) '_K' num2str(K) '.mat'])

% first make figures to ensure that hmmm results have been properly
% concatenated with protein data
if make_trace_plots
    for i = 1:numel(master_struct)
        subProject = master_struct(i).project;
        subID = master_struct(i).ID;
        qcPath = [figPath '\' subID '_qc_' subProject '\'];
        mkdir(qcPath);
        hmm_input_output = master_struct(i).hmm_input_output;
        s_index = 1:numel(hmm_input_output);
        rng(123);
        plot_indices = randsample(s_index,min([20,numel(s_index)]),false);
        for j = 1:numel(plot_indices)
            % MCP channel checks
            mcp_check = hmm_input_output(plot_indices(j)).mcp_check;
            fluo_check = hmm_input_output(plot_indices(j)).fluo_check;
            fluo = hmm_input_output(plot_indices(j)).fluo;
            time = hmm_input_output(plot_indices(j)).time;
            r_vec = sum(hmm_input_output(plot_indices(j)).r_mat,2);
            % make figure
            qc_fig = figure('Visible','off');
            hold on
            plot(time,fluo / nanmean(fluo))
            plot(time,fluo_check / nanmean(fluo_check));
            plot(time,mcp_check / nanmean(mcp_check))
            plot(time,r_vec / nanmean(r_vec))
            legend('fluo (HMM)', 'fluo (data)','raw mcp','activity state (HMM)')
            xlabel('time')
            ylabel([gene_name ' activity (au)'])
            saveas(qc_fig,[qcPath 'mcp_check_nc_' sprintf('%03d',plot_indices(j)) '.png'])

            % Protein Channel checks
            spot_protein = hmm_input_output(plot_indices(j)).spot_protein;
            null_protein = hmm_input_output(plot_indices(j)).null_protein;
            mf_protein = hmm_input_output(plot_indices(j)).mf_protein;
            % make figure
            qc_fig = figure('Visible','off');
            hold on
            plot(time,spot_protein)
            plot(time,null_protein)
            plot(time,mf_protein)
            legend('protein (spot)', 'protein (control spot)','protein (mf control)')
            xlabel('time')
            ylabel([protein_name ' - ' protein_fluor ' (au)'])
            saveas(qc_fig,[qcPath 'protein_check_nc_' sprintf('%03d',plot_indices(j)) '.png'])
        end
    end

    % Make single trace input-output plots
    for i = 1:numel(master_struct)
        subProject = master_struct(i).project;
        subID = master_struct(i).ID;
        tracePath = [figPath '\' subID '_single_trace_' subProject '\'];
        mkdir(tracePath);
        hmm_input_output = master_struct(i).hmm_input_output;
        s_index = 1:numel(hmm_input_output);
        rng(321);
        plot_indices = randsample(s_index,min([nTraces,numel(s_index)]),false);
        for j = 1:numel(plot_indices)
            % MCP channel checks
            time = hmm_input_output(plot_indices(j)).time;
            r_vec = sum(hmm_input_output(plot_indices(j)).r_mat,2);
            spot_protein = hmm_input_output(plot_indices(j)).spot_protein;
            null_protein = hmm_input_output(plot_indices(j)).null_protein;
            delta_protein = spot_protein - null_protein;
            % make figure
            trace_fig = figure('Visible','off');
            hold on

            yyaxis left
            plot(time,r_vec)
            ylabel(['instantaneous ' gene_name ' activity (au)'])

            yyaxis right
            plot(time,delta_protein);
            ylabel(['absolute ' protein_name ' enrichment(au)'])

            ax = gca;
            ax.YAxis(1).Color = 'black';
            ax.YAxis(2).Color = 'black';

            legend('transcriptional activity', 'local protein concentration')
            xlabel('time')
            ylabel([gene_name ' activity (au)'])
            saveas(trace_fig,[tracePath 'input_output_nc' sprintf('%03d',plot_indices(j)) '.png'])               
        end
    end
end

close all
% ignore possibility of second project for noq 
hmm_input_output = master_struct(1).hmm_input_output;
Tres = nanmedian(diff([hmm_input_output.time]));

% Generate effective transcription rate vector and de-trend protein data
occ_vec = nanmean(vertcat(hmm_input_output.z_mat));
r_mean = nanmean(vertcat(hmm_input_output.r));
if K==3
    r_mean = [r_mean(1) (r_mean(2).*occ_vec(2)+r_mean(3).*occ_vec(3))/(occ_vec(2)+occ_vec(3))];
end
for i = 1:numel(hmm_input_output)
    % generate discrete state and continuous production vectors
    z_mat = hmm_input_output(i).z_mat;
    r = hmm_input_output(i).r;    
    if K == 3
        z_mat(:,2) = sum(z_mat(:,2:3),2);
        z_mat = z_mat(:,1:2);
    end
    r_mat = z_mat .* r_mean;
    hmm_input_output(i).r_vec = sum(r_mat,2);
    hmm_input_output(i).z_mat2 = z_mat;
    [~,z_vec] = max(z_mat,[],2);
    hmm_input_output(i).z_vec = z_vec;       
end
% obtain time-averaged fluo trend
t_vec = round([hmm_input_output.time]/60);
f_vec = [hmm_input_output.fluo];
r_vec = vertcat(hmm_input_output.r_vec);

time_index = 1:50;
fluo_trend = NaN(size(time_index));
r_trend = NaN(size(time_index));
for t = 1:numel(time_index)
    fluo_trend(t) = nanmean(f_vec(t_vec==time_index(t)));
    r_trend(t) = nanmean(r_vec(t_vec==time_index(t)));
end
filter_size = 5;
% detrend input and output time series
for i = 1:numel(hmm_input_output)
    % detrend spot protein
    spot_vec = hmm_input_output(i).spot_protein;
    frame_vec = 1:numel(spot_vec);
    [p,~,mu] = polyfit(frame_vec(~isnan(spot_vec)),spot_vec(~isnan(spot_vec)),2);
    spot_trend = polyval(p,(1:numel(spot_vec)),[],mu);
    hmm_input_output(i).spot_protein_dt = spot_vec - spot_trend;
    hmm_input_output(i).spot_protein_dt_smooth = imgaussfilt(hmm_input_output(i).spot_protein_dt,filter_size);
    % detrend control 
    null_vec = hmm_input_output(i).null_protein;
    [p,~,mu] = polyfit(frame_vec(~isnan(null_vec)),null_vec(~isnan(null_vec)),2);
    null_trend = polyval(p,(1:numel(null_vec)),[],mu);
    hmm_input_output(i).null_protein_dt = null_vec - null_trend;
    % add smoothed version of control 
    hmm_input_output(i).null_protein_sm = imgaussfilt(hmm_input_output(i).null_protein_all,filter_size);
    hmm_input_output(i).delta_protein_sm = spot_vec - hmm_input_output(i).null_protein_sm;
    % detrend fluo
    fluo = hmm_input_output(i).fluo;
    activity = hmm_input_output(i).r_vec;
    time = round(hmm_input_output(i).time/60);
    for j = 1:numel(fluo)
        fluo(j) = fluo(j) - fluo_trend(time(j)==time_index);
        activity(j) = activity(j) - r_trend(time(j)==time_index);
    end
    hmm_input_output(i).fluo_dt = fluo;
    hmm_input_output(i).r_vec_dt = activity;
end


%%% Simple(ish) things first...
prctile_vec = 0:10:100;
spot_prctile_vec = NaN(size(prctile_vec));
null_prctile_vec = NaN(size(prctile_vec));
spot_protein = [hmm_input_output.spot_protein];
null_protein_smooth = [hmm_input_output.null_protein_sm];
fluo_vec = [hmm_input_output.fluo];
activity_vec = vertcat(hmm_input_output.r_vec)';

for i = 1:numel(prctile_vec)
    spot_prctile_vec(i) = prctile(spot_protein,prctile_vec(i));
    null_prctile_vec(i) = prctile(null_protein_smooth,prctile_vec(i));
end
[bin_counts, ~,~,binSpot,binNull] = histcounts2(spot_protein,null_protein_smooth,spot_prctile_vec,null_prctile_vec);
norm_counts = bin_counts ./ repmat(sum(bin_counts),numel(prctile_vec)-1,1);

resample_weights =  repmat(sum(norm_counts,2),1,numel(prctile_vec)-1) ./ norm_counts;

% now draw fluo, delta, and mf samples for each delta bin
null_boot_array = NaN(nBoots,numel(prctile_vec)-1);
spot_boot_array = NaN(nBoots,numel(prctile_vec)-1);
fluo_boot_array = NaN(nBoots,numel(prctile_vec)-1);
activity_boot_array = NaN(nBoots,numel(prctile_vec)-1);

for i = 1:numel(prctile_vec)-1
    bin_indices = find(binSpot==i);
    resample_wt_vec = resample_weights(:,i);
    sample_weights = resample_wt_vec(binNull(bin_indices));
    for n = 1:nBoots
        s_ids = randsample(bin_indices,numel(bin_indices),true,sample_weights);
        null_boot_array(n,i) = nanmean(null_protein_smooth(s_ids));
        spot_boot_array(n,i) = nanmean(spot_protein(s_ids));
        fluo_boot_array(n,i) = nanmean(fluo_vec(s_ids));
        activity_boot_array(n,i) = nanmean(activity_vec(s_ids));
    end
end

fluo_mean = nanmean(fluo_boot_array);
fluo_ste= nanstd(fluo_boot_array);

activity_mean = nanmean(activity_boot_array);
activity_ste= nanstd(activity_boot_array);

mf_mean = nanmean(null_boot_array);
mf_ste= nanstd(null_boot_array);

spot_mean = nanmean(spot_boot_array);
spot_ste= nanstd(spot_boot_array);

aspirational_fig = figure('Visible','off');
hold on
errorbar(spot_mean,activity_mean,activity_ste);
ylabel([gene_name ' activity']);
yyaxis right
plot(spot_mean,spot_mean)
hold on
plot(spot_mean,mf_mean)
ylabel([protein_name ' concentration (au)'])

legend([gene_name ' activity'], [protein_name ' (local)'], [protein_name ' (background)'],'Location','southeast')
xlabel(['local ' protein_name ' concentration (au)'])
saveas(aspirational_fig, [figPath 'aspirational_fig.png'])

%%% Examine sna activity in vicinity of anamalous high and low pt points

spot_protein_dt_smooth = [hmm_input_output.spot_protein_dt_smooth];
null_protein_smooth = [hmm_input_output.null_protein_sm];

% calculate necessary sample weights
bulk_prctile_index = 0:5:100;
bulk_pt_prctile_vec = NaN(size(bulk_prctile_index));
for i = 1:numel(bulk_prctile_index)
    bulk_pt_prctile_vec(i) = prctile(spot_protein_dt_smooth,bulk_prctile_index(i));
end

nSamples = nansum(binSpot==1);

% initialize arrays
fluo_cell = cell(1,3);
local_protein_cell = cell(1,3);
null_protein_cell = cell(1,3);
activity_cell = cell(1,3);
ind_ref_cell = cell(1,3);

% pull time series samples from vicinity of protein events
name_cell = {'low','mid','high'};
ref_vec = -window_size:window_size;%1:(2*window_size+1);
for i = 1:numel(hmm_input_output)
    time_vec = hmm_input_output(i).time/60;
    frame_vec = 1:numel(time_vec);
    fluo_vec = hmm_input_output(i).fluo;
    activity_vec = hmm_input_output(i).r_vec;
    % find protein peaks and troughs
    pt_sp_dt_sm = hmm_input_output(i).spot_protein_dt_smooth;   
    pt_sp = hmm_input_output(i).spot_protein;   
    pt_nn_sm = hmm_input_output(i).null_protein_sm;   
    [~,~,nullBins] = histcounts(pt_sp_dt_sm,bulk_pt_prctile_vec);
    nan_filter = ~isnan(pt_nn_sm)&~isnan(pt_sp_dt_sm)&frame_vec>3&frame_vec<numel(frame_vec)-1;
    % find points ot interest
    high_ids = find(nan_filter&nullBins==numel(bulk_pt_prctile_vec)-1);
    low_ids = find(nan_filter&nullBins==1);
    mid_ids = find(nan_filter&nullBins>5&nullBins<16);

    % sample proiten and activity traces
    for j = 1:numel(name_cell)
        eval(['ids =' name_cell{j} '_ids;'])
        fluo_temp = NaN(numel(ids),2*window_size+1);
        act_temp = NaN(numel(ids),2*window_size+1);
        loc_temp = NaN(numel(ids),2*window_size+1);
        bkg_temp = NaN(numel(ids),2*window_size+1);
        ind_temp = NaN(1,numel(ids));
        for k = 1:numel(ids)
            raw_ind = ref_vec+ids(k);
            ft_vec1 = raw_ind > 0 & raw_ind <= numel(fluo_vec);
            ft_vec2 = ref_vec(ft_vec1) + ids(k);
            % record
            fluo_temp(k,ft_vec1) = fluo_vec(ft_vec2);
            act_temp(k,ft_vec1) = activity_vec(ft_vec2);
            loc_temp(k,ft_vec1) = pt_sp_dt_sm(ft_vec2);
            bkg_temp(k,ft_vec1) = pt_nn_sm(ft_vec2);
            ind_temp(k) = i;
            if isnan(bkg_temp(k,window_size+1))
                error('wtf')
            end
        end
        % add to main cell structures
        fluo_cell{j} = vertcat(fluo_cell{j}, fluo_temp);
        activity_cell{j} = vertcat(activity_cell{j}, act_temp);
        local_protein_cell{j} = vertcat(local_protein_cell{j}, loc_temp);
        null_protein_cell{j} = vertcat(null_protein_cell{j}, bkg_temp);
        ind_ref_cell{j} =[ind_ref_cell{j} ind_temp];
    end
end
tic
% now construct control sets from mid traces that mimic temporal mf profile
% of high and low traces
trend_id_vec = [1 2 3];
control_id_vec = [2 2 2];
sampling_struct = struct;
for i = 1:numel(trend_id_vec)
    % get background pt levels
    nn_mat_trend = null_protein_cell{trend_id_vec(i)};
    nn_mat_control = null_protein_cell{control_id_vec(i)};
    ind_vec_trend = ind_ref_cell{trend_id_vec(i)};
    ind_vec_control = ind_ref_cell{control_id_vec(i)};
    % other control arrays
    sp_mat_control = local_protein_cell{control_id_vec(i)};    
    fluo_mat_control = fluo_cell{control_id_vec(i)};        
    act_mat_control = activity_cell{control_id_vec(i)};
    % initialize targeted control arrays
    nn_spec_control = NaN(size(nn_mat_trend));
    sp_spec_control = NaN(size(nn_mat_trend));
    fluo_spec_control = NaN(size(nn_mat_trend));
    act_spec_control = NaN(size(nn_mat_trend));
    % ok. now iterate through trend set and find closest match to each
    % trace in control
    ctrl_ids = NaN(1,size(nn_mat_trend,1));
    for j = 1:size(nn_mat_trend,1)
        bkg_vec = nn_mat_trend(j,:);
        bkg_ind = ind_vec_trend(j);
        bkg_mat = repmat(bkg_vec,size(nn_mat_control,1),1);
        nan_ft = ~isnan(bkg_vec);
        overlap_vec = sum(~isnan(bkg_mat)==~isnan(nn_mat_control),2);
        diff_vec = nanmean(abs(nn_mat_control-bkg_vec),2);
        if trend_id_vec(i) == control_id_vec(i)
            diff_vec(j) = Inf;
        end     
        diff_vec(overlap_vec~=sum(nan_ft)) = Inf;
        diff_vec(ind_vec_control==bkg_ind) = Inf;
        [~, mi] = min(diff_vec);
        ctrl_ids(j) = mi;
        nn_spec_control(j,nan_ft) = nn_mat_control(mi,nan_ft);
        sp_spec_control(j,nan_ft) = sp_mat_control(mi,nan_ft);
        fluo_spec_control(j,nan_ft) = fluo_mat_control(mi,nan_ft);
        act_spec_control(j,nan_ft) = act_mat_control(mi,nan_ft);
    end
    % save relevant arrayes in data structure
    sampling_struct(i).trend_id = trend_id_vec(i);
    sampling_struct(i).trend_name = name_cell{trend_id_vec(i)};
    sampling_struct(i).ctrl_id = control_id_vec(i);
    sampling_struct(i).ctrl_name = name_cell{control_id_vec(i)};
    sampling_struct(i).ctrl_id_vec = ctrl_ids;
    % save trend data fields
    sampling_struct(i).trend_null_protein = null_protein_cell{trend_id_vec(i)};
    sampling_struct(i).trend_spot_protein = local_protein_cell{trend_id_vec(i)};
    sampling_struct(i).trend_fluo = fluo_cell{trend_id_vec(i)};
    sampling_struct(i).trend_activity = activity_cell{trend_id_vec(i)};
    % controls
    sampling_struct(i).ctrl_null_protein = nn_spec_control;
    sampling_struct(i).ctrl_spot_protein = sp_spec_control;
    sampling_struct(i).ctrl_fluo = fluo_spec_control;
    sampling_struct(i).ctrl_activity = act_spec_control;
end
toc    
%%
