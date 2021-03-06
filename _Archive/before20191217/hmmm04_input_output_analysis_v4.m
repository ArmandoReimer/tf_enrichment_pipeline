% Script to probe relationship between input protein concentration and
% output transcriptional response
clear 
close all
% define ID variables
K = 3;
w = 7;
project = 'Dl_Venus_snaBAC_mCherry_Leica_hp';
% dropboxFolder =  'E:\Nick\Dropbox (Garcia Lab)\';
dropboxFolder = 'C:\Users\nlamm\Dropbox (Garcia Lab)\';
dataPath = [dropboxFolder '\ProcessedEnrichmentData\' project '\'];
figPath = [dropboxFolder '\LocalEnrichmentFigures\' project '\hmm_input_output_K' num2str(K) '_w' num2str(w) '\'];
mkdir(figPath)

nTraces = 50; % number of individual traces to select for plotting
window_size = 10; % number of lags and leads over which to track protein/fluo dynamics
nBoots = 100;
min_time = 5; % minimum # minutes into trace to take
make_trace_plots = 0;
dTime = 2; % size of temporal deviation permitted btw trend and control
feature_filter_size = 3;
bkg_filter_size = 9;
maxSamples = 4000;
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

% detrend input and output time series
for i = 1:numel(hmm_input_output)
    % detrend spot protein
    spot_vec = hmm_input_output(i).spot_protein;
    frame_vec = 1:numel(spot_vec);
    [p,~,mu] = polyfit(frame_vec(~isnan(spot_vec)),spot_vec(~isnan(spot_vec)),2);
    spot_trend = polyval(p,(1:numel(spot_vec)),[],mu);
    hmm_input_output(i).spot_protein_dt = hmm_input_output(i).spot_protein - spot_trend;
    hmm_input_output(i).spot_protein_smooth = imgaussfilt(hmm_input_output(i).spot_protein_all,feature_filter_size,'FilterDomain','spatial');
    hmm_input_output(i).spot_protein_smooth(isnan(spot_vec)) = NaN;    
    % add smoothed version of control 
    hmm_input_output(i).null_protein_smooth = imgaussfilt(hmm_input_output(i).null_protein_all,feature_filter_size,'FilterDomain','spatial');
    hmm_input_output(i).null_protein_bkg = imgaussfilt(hmm_input_output(i).null_protein_all,bkg_filter_size,'FilterDomain','spatial');
    hmm_input_output(i).null_protein_smooth(isnan(spot_vec)) = NaN;
    hmm_input_output(i).null_protein_bkg(isnan(spot_vec)) = NaN;
    hmm_input_output(i).delta_protein_smooth = hmm_input_output(i).spot_protein_smooth - hmm_input_output(i).null_protein_bkg;
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
null_protein_smooth = [hmm_input_output.null_protein_smooth];
fluo_vec = [hmm_input_output.fluo];
activity_vec = vertcat(hmm_input_output.r_vec)';

for i = 1:numel(prctile_vec)
    spot_prctile_vec(i) = prctile(spot_protein,prctile_vec(i));
    null_prctile_vec(i) = prctile(null_protein_smooth,prctile_vec(i));
end
[bin_counts, ~,~,binSpot,binNull] = histcounts2(spot_protein,null_protein_smooth,spot_prctile_vec,null_prctile_vec);
norm_counts = bin_counts ./ repmat(sum(bin_counts),numel(prctile_vec)-1,1);
% calculate weights needed to mimic generall population
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Examine sna activity in vicinity of anamalous high and low pt points

delta_protein_smooth = [hmm_input_output.delta_protein_smooth];
null_protein_smooth = [hmm_input_output.null_protein_bkg];

% calculate necessary sample 
bulk_prctile_index = 0:5:100;

bulk_pt_prctile_vec = NaN(size(bulk_prctile_index));
for i = 1:numel(bulk_prctile_index)
    bulk_pt_prctile_vec(i) = prctile(delta_protein_smooth,bulk_prctile_index(i));
end

% initialize arrays
fluo_cell = cell(1,3);
local_protein_cell = cell(1,3);
null_protein_cell = cell(1,3);
activity_cell = cell(1,3);
ind_ref_cell = cell(1,3);
time_cell = cell(1,3);

% pull time series samples from vicinity of protein events
name_cell = {'low','mid','high'};
ref_vec = -window_size:window_size;%1:(2*window_size+1);
for i = 1:numel(hmm_input_output)
    % core features
    time_vec = hmm_input_output(i).time/60;   
    frame_vec = 1:numel(time_vec);
    fluo_vec = hmm_input_output(i).fluo;
    activity_vec = hmm_input_output(i).r_vec;
    dt_sm = hmm_input_output(i).delta_protein_smooth;      
    pt_nn_sm = hmm_input_output(i).null_protein_bkg; 
    % find protein peaks and troughs
    [~,high_ids,blip_widths] = findpeaks(dt_sm,frame_vec,'MinPeakWidth',60*.5/Tres,'MinPeakProminence',.3); % NL: eye-balled atm
    [~,low_ids,dip_widths] = findpeaks(-dt_sm,frame_vec,'MinPeakWidth',60*.5/Tres,'MinPeakProminence',.3); % NL: eye-balled atm
    
    nan_filter = ~isnan(pt_nn_sm)&frame_vec>2&frame_vec<numel(frame_vec)-1;
    nan_ids = find(nan_filter);
    high_ids = high_ids(ismember(high_ids,nan_ids));
    low_ids = low_ids(ismember(low_ids,nan_ids));
    
    mid_ids = find(nan_filter&~ismember(frame_vec,[high_ids high_ids-1 high_ids+1])&...
        ~ismember(frame_vec,[low_ids low_ids-1 low_ids+1]));

    % sample proiten and activity traces
    for j = 1:numel(name_cell)
        eval(['ids =' name_cell{j} '_ids;'])
        fluo_temp = NaN(numel(ids),2*window_size+1);
        act_temp = NaN(numel(ids),2*window_size+1);
        loc_temp = NaN(numel(ids),2*window_size+1);
        bkg_temp = NaN(numel(ids),2*window_size+1);
        ind_temp = NaN(1,numel(ids));
        time_temp = NaN(1,numel(ids));
        for k = 1:numel(ids)
            raw_ind = ref_vec+ids(k);
            ft_vec1 = raw_ind > 0 & raw_ind <= numel(fluo_vec);
            ft_vec2 = raw_ind(ft_vec1);
            % record
            fluo_temp(k,ft_vec1) = fluo_vec(ft_vec2);
            act_temp(k,ft_vec1) = activity_vec(ft_vec2);
            loc_temp(k,ft_vec1) = dt_sm(ft_vec2);
            bkg_temp(k,ft_vec1) = pt_nn_sm(ft_vec2);
            ind_temp(k) = i;
            time_temp(k) = time_vec(ids(k));
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
        time_cell{j} =[time_cell{j} time_temp];
    end
end
tic
% now construct control sets from mid traces that mimic temporal mf profile
% of high and low traces
trend_id_vec = [1 2 3];
control_id_vec = [2 2 2];
sampling_struct = struct;
for i = 1:numel(trend_id_vec)
    % get dore matching vars
    nn_mat_trend = null_protein_cell{trend_id_vec(i)};
    time_vec_trend = time_cell{trend_id_vec(i)};
    ind_vec_trend = ind_ref_cell{trend_id_vec(i)};
    
    nn_mat_control = null_protein_cell{control_id_vec(i)};
    time_vec_control = time_cell{control_id_vec(i)};    
    ind_vec_control = ind_ref_cell{control_id_vec(i)};
    % other control arrays
    loc_mat_control = local_protein_cell{control_id_vec(i)};    
    fluo_mat_control = fluo_cell{control_id_vec(i)};        
    act_mat_control = activity_cell{control_id_vec(i)};
    % initialize targeted control arrays
    nRows = min([maxSamples size(nn_mat_trend,1)]);
    nn_matched_control = NaN(nRows,size(nn_mat_trend,2));
    sp_matched_control = NaN(nRows,size(nn_mat_trend,2));
    fluo_matched_control = NaN(nRows,size(nn_mat_trend,2));
    act_matched_control = NaN(nRows,size(nn_mat_trend,2));
    % ok. now iterate through trend set and find closest match to each
    % trace in control
    ctrl_ids = NaN(1,nRows);
    trend_ids = NaN(1,nRows);
    index_vec = 1:size(nn_mat_trend,1);
    shuffled_indices = randsample(index_vec,numel(index_vec),false);
    parfor j = 1:nRows
        rowInd = shuffled_indices(j);
        trend_ids(j) = rowInd;
        trend_bkg_vec = nn_mat_trend(rowInd,:);
        trend_ind = ind_vec_trend(rowInd);
        trend_time = time_vec_trend(rowInd);
        % apply selection criteria
        trend_bkg_mat = repmat(trend_bkg_vec,size(nn_mat_control,1),1);
        nn_bkg_control_temp = nn_mat_control;
        nn_bkg_control_temp(isnan(trend_bkg_mat)) = NaN;
        nan_ft = ~isnan(trend_bkg_vec);
        overlap_vec = sum(~isnan(trend_bkg_mat)==~isnan(nn_bkg_control_temp),2);
        t_diff_vec = abs(time_vec_control - trend_time);        
        pt_diff_vec = nanmean(abs(nn_bkg_control_temp-trend_bkg_vec),2);
        
        % must have same pattern of missing and prese
        pt_diff_vec(overlap_vec~=sum(nan_ft)) = Inf;
        % cannot be from same trace
        pt_diff_vec(ind_vec_control==trend_ind) = Inf;
        % must fall during comparable period of time
        pt_diff_vec(t_diff_vec>dTime) = Inf;
        if sum(~isinf(pt_diff_vec))==0
            warning('no matches found under constraints')
            continue
        end
        % find closest match
        [~, mi] = min(pt_diff_vec);
        ctrl_ids(j) = mi;
        nn_blank = nn_mat_control(mi,:);
        nn_blank(~nan_ft) = NaN;
        nn_matched_control(j,:) = nn_blank;
        
        sp_blank = loc_mat_control(mi,:);
        sp_blank(~nan_ft) = NaN;
        sp_matched_control(j,:) = sp_blank;
        
        fluo_blank = fluo_mat_control(mi,:);
        fluo_blank(~nan_ft) = NaN;
        fluo_matched_control(j,:) = fluo_blank;
        
        act_blank = act_mat_control(mi,:);
        act_blank(~nan_ft) = NaN;
        act_matched_control(j,:) = act_blank;
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
    sampling_struct(i).trend_id_vec = trend_ids;    
    % controls
    sampling_struct(i).ctrl_null_protein = nn_matched_control;
    sampling_struct(i).ctrl_spot_protein = sp_matched_control;
    sampling_struct(i).ctrl_fluo = fluo_matched_control;
    sampling_struct(i).ctrl_activity = act_matched_control;
end
toc    
%% Calculate bootstrap estimates of Mean and SE
calc_fields = {'null_protein','spot_protein','fluo','activity'};
id_types = {'trend','ctrl'};
    
trend_activity_se_array = NaN(2*window_size+1,numel(name_cell));

for i = 1:numel(sampling_struct)
    qc_ctrl_id_vec = sampling_struct(i).ctrl_id_vec;
    qc_trend_id_vec = sampling_struct(i).trend_id_vec;
    qc_trend_indices = qc_trend_id_vec(~isnan(qc_ctrl_id_vec));
    for j = 1:numel(calc_fields)
        f_str = calc_fields{j};
        for k = 1:numel(id_types)
            id_str = id_types{k};
            data_str = [id_str '_' f_str];
            data_mat = sampling_struct(i).(data_str);
            if strcmpi(id_str,'trend')
                data_mat = data_mat(qc_trend_indices,:);
            end
            boot_mat = NaN(nBoots,2*window_size+1);
            index_vec = 1:size(data_mat,1);
            for n = 1:nBoots
                s_ids = randsample(index_vec,numel(index_vec),true);
                boot_mat(n,:) = nanmean(data_mat(s_ids,:));
            end
            sampling_struct(i).([data_str '_mean']) = nanmean(boot_mat);
            sampling_struct(i).([data_str '_se']) = nanstd(boot_mat);
        end
    end
end
                       
time_axis = Tres*ref_vec / 60;
cm = jet(128);
red1 = cm(115,:);
red2 = cm(128,:);
blue1 = cm(35,:);
blue2 = cm(20,:);

for i = 1:numel(name_cell)
   
    qc_fig = figure;
    hold on
    errorbar(time_axis,sampling_struct(i).trend_null_protein_mean,...
        sampling_struct(i).trend_null_protein_se,'Color',blue1,'LineWidth',1.2,'CapSize',0)
    plot(time_axis,sampling_struct(i).ctrl_null_protein_mean,'--','Color',blue2)
%     errorbar(time_axis,sampling_struct(i).ctrl_null_protein_mean,...
%         sampling_struct(i).ctrl_null_protein_se,'--','Color',blue2,'CapSize',0)
    ylabel([protein_name ' concentration (au)'])

    xlabel('offset')
    ax = gca;
    ax.YColor = 'black';
    legend('background protein (control)','background protein (trend)');%,'fluorescence (control)','fluorescence (trend)')
    grid on

    saveas(qc_fig,[figPath name_cell{i} '_sampling_qc.png'])
end

close all
% Now make figures examining input-output relationships
for i = 1:numel(name_cell)
    % figure to check consistency of method

    in_out_fig = figure;
    hold on
    yyaxis left
    plot(time_axis,sampling_struct(i).trend_spot_protein_mean,'Color',blue2,'LineWidth',1)
    plot(time_axis,sampling_struct(i).ctrl_spot_protein_mean,'-','Color',[blue1 .3])
    ylabel([protein_name ' concentration (au)'])
    ax = gca;
    ax.YColor = blue2;
    
    
    yyaxis right
    errorbar(time_axis,sampling_struct(i).trend_activity_mean,...
        sampling_struct(i).trend_activity_se,'-','Color','black','LineWidth',1,'CapSize',0)
    plot(time_axis,sampling_struct(i).ctrl_activity_mean,'-','Color',[0 0 0 .3])
    ylabel([protein_name ' concentration (au)'])
    ylabel([gene_name ' activity (au)'])

    xlabel('offset (minutes)')
    ax = gca;
    ax.YColor = 'black';
    legend('protein (trend)','protein (control)','production rate (trend)',...
        'production rate(control)','Location','southwest')
    grid on

    saveas(in_out_fig,[figPath name_cell{i} '_input_output_act.png'])
end