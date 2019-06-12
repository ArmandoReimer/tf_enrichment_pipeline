clear
close all

% define core ID variables
project = 'Dl-Ven_snaBAC-mCh';
% project = 'Dl-Ven_hbP2P-mCh';
dropboxFolder =  'E:\Nick\Dropbox (Garcia Lab)\';
dataPath = [dropboxFolder 'ProcessedEnrichmentData\' project '/'];
w = 7;
K = 3;
% window analysis params
window_size = 15; 
% load input-output data set
load([dataPath 'hmm_input_output_w' num2str(w) '_K' num2str(K) '.mat'],'hmm_input_output')
for i = 1:numel(hmm_input_output)    
    z_vec = hmm_input_output(i).z_vec' > 1;
    r_vec = hmm_input_output(i).r_vec';
    hmm_input_output(i).z_vec = z_vec;
    z_prob_vec = sum(hmm_input_output(i).z_mat(:,2:3),2);
    zd = [0 diff(z_vec)];
    change_points = find(zd~=0);
    % calculate duration of bursts
    dur_vec_lag = diff([change_points NaN]);
    dur_vec_lead = diff([NaN change_points]);
    % calculate intensity of bursts
    sz_vec_lag = NaN(size(dur_vec_lag));
    sz_vec_lead = NaN(size(dur_vec_lag));
    for j = 1:numel(change_points)-1
        r_mean = nanmean(r_vec(change_points(j):change_points(j+1)-1));
        sz_vec_lag(j) = r_mean;
        sz_vec_lead(j+1) = r_mean;
    end         
    % generate full-length vectors
    z_dur_lag_vec = NaN(size(z_vec));
    z_dur_lead_vec = NaN(size(z_vec));
    z_dur_lag_vec(change_points) = dur_vec_lag;  
    z_dur_lead_vec(change_points) = dur_vec_lead;  
    sz_lag_vec = NaN(size(z_vec));
    sz_lead_vec = NaN(size(z_vec));
    sz_lag_vec(change_points) = sz_vec_lag;  
    sz_lead_vec(change_points) = sz_vec_lead;  
    % record
    hmm_input_output(i).z_dur_lag_vec = z_dur_lag_vec;  
    hmm_input_output(i).z_dur_lead_vec = z_dur_lead_vec;
    hmm_input_output(i).sz_lag_vec = sz_lag_vec;  
    hmm_input_output(i).sz_lead_vec = sz_lead_vec;
    hmm_input_output(i).z_diff_vec = zd;
    hmm_input_output(i).z_prob_vec = z_prob_vec';
    % detrend data
    spot_pt_vec = hmm_input_output(i).spot_protein;%(~dt_filter_gap);
    time_vec = hmm_input_output(i).time;%(~dt_filter_gap);
    p = polyfit(time_vec,spot_pt_vec,2);
    spot_pt_trend = polyval(p,time_vec);
    hmm_input_output(i).spot_protein_dt = spot_pt_vec - spot_pt_trend; 
    % virtual spot protein
    virtual_pt_vec = hmm_input_output(i).serial_protein;%(~dt_filter_gap);    
    p = polyfit(time_vec,virtual_pt_vec,2);
    virtual_pt_trend = polyval(p,time_vec);
    hmm_input_output(i).serial_protein_dt = virtual_pt_vec - virtual_pt_trend;
    % swap spot protein
    swap_pt_vec = hmm_input_output(i).swap_spot_protein;%(~dt_filter_gap);    
    p = polyfit(time_vec,swap_pt_vec,2);
    swap_pt_trend = polyval(p,time_vec);
    hmm_input_output(i).swap_spot_protein_dt = swap_pt_vec - swap_pt_trend;
end

% generate master set of vectors for feature classification
gap_filter_vec = [hmm_input_output.dt_filter_gap];
time_vec = round([hmm_input_output.time]/60);
z_dur_lag_list = [hmm_input_output.z_dur_lag_vec];
z_dur_lead_list = [hmm_input_output.z_dur_lead_vec];
sz_lag_list = [hmm_input_output.sz_lag_vec];
sz_lead_list = [hmm_input_output.sz_lead_vec];
z_diff_list = [hmm_input_output.z_diff_vec];
% initialize results structure
results_struct = struct;
% arrays for fitting linear offsets
window_vec = -window_size:window_size;
fit_array = [ones(numel(window_vec),1) window_vec'];
   
% calculate expected number of features for pre-allocation
n_entries = sum(z_diff_list~=0&~gap_filter_vec);
% initialize arrays to store time series snip
fluo_array = NaN(n_entries,2*window_size+1);
hmm_array = NaN(n_entries,2*window_size+1);
spot_array = NaN(n_entries,2*window_size+1);
swap_array = NaN(n_entries,2*window_size+1);
virtual_array = NaN(n_entries,2*window_size+1);        
% initialize arrays to store id variables
particle_id_vec = NaN(1,n_entries);
center_time_vec = NaN(1,n_entries);
lag_dur_vec = NaN(1,n_entries);
lead_dur_vec = NaN(1,n_entries);
lag_size_vec = NaN(1,n_entries);
lead_size_vec = NaN(1,n_entries);
feature_sign_vec = NaN(1,n_entries);

% iterate through structure
iter = 1;                    
for j = 1:numel(hmm_input_output)
    % core ID vectors
    ParticleID = hmm_input_output(j).ParticleID;
    time = hmm_input_output(j).time;
    % feature classification vectors
    gap_filter = hmm_input_output(j).dt_filter_gap;      
    z_dur_lag_vec = hmm_input_output(j).z_dur_lag_vec;            
    z_dur_lead_vec = hmm_input_output(j).z_dur_lead_vec;  
    sz_lag_vec = hmm_input_output(j).sz_lag_vec;            
    sz_lead_vec = hmm_input_output(j).sz_lead_vec;
    z_diff_vec = hmm_input_output(j).z_diff_vec;  
    % activity
    fluo = hmm_input_output(j).fluo;
    r_vec = hmm_input_output(j).r_vec';
    % protein fields                        
    spot_protein = hmm_input_output(j).spot_protein_dt;
    swap_spot_protein = hmm_input_output(j).swap_spot_protein_dt;
    virtual_protein = hmm_input_output(j).serial_protein_dt;        
    % apply filter             
    spot_protein(gap_filter) = NaN;
    swap_spot_protein(gap_filter) = NaN;
    virtual_protein(gap_filter) = NaN;
    % find features
    id_list = find(z_diff_vec~=0&~gap_filter);
    for id = id_list
        full_range = id - window_size:id+window_size;
        true_range = full_range(full_range>0&full_range<=numel(virtual_protein));
        % record
        ft1 = ismember(full_range,true_range);
        % qc check
        if sum(~isnan(spot_protein(true_range))) >= window_size && sum(~isnan(swap_spot_protein(true_range)))...
                >= window_size && sum(~isnan(virtual_protein(true_range))) >= window_size
            % extract raw fragments
            spot_fragment = spot_protein(true_range);
            swap_fragment = swap_spot_protein(true_range);
            virtual_fragment = virtual_protein(true_range);
            fluo_fragment = fluo(true_range);
            hmm_fragment = r_vec(true_range);
            % fit linear offsets            
            fit_sub_array = fit_array(ft1,:);
%             spot_fit = fit_sub_array(~isnan(spot_fragment),:) \ spot_fragment(~isnan(spot_fragment))';
%             swap_fit = fit_sub_array(~isnan(swap_fragment),:) \ swap_fragment(~isnan(swap_fragment))';
%             virtual_fit = fit_sub_array(~isnan(virtual_fragment),:) \ virtual_fragment(~isnan(virtual_fragment))';                
            fluo_fit = fit_sub_array(~isnan(fluo_fragment),:) \ fluo_fragment(~isnan(fluo_fragment))';
%             % save time snips               
            spot_array(iter,ft1) = spot_fragment;
            swap_array(iter,ft1) = swap_fragment;
            virtual_array(iter,ft1) = virtual_fragment;
            fluo_array(iter,ft1) = fluo_fragment - fluo_fit(1) - fluo_fit(2)*window_vec(ft1);
            hmm_array(iter,ft1) = hmm_fragment;
            % save other info
            particle_id_vec(iter) = ParticleID;
            center_time_vec(iter) = time(id);
            lag_dur_vec(iter) = z_dur_lag_vec(id);
            lead_dur_vec(iter) = z_dur_lead_vec(id);
            lag_size_vec(iter) = sz_lag_vec(id);
            lead_size_vec(iter) = sz_lead_vec(id);
            feature_sign_vec(iter) = sign(z_diff_vec(id));
        end
        % increment
        iter = iter + 1;  
        if mod(iter,100) == 0
            disp(iter)
        end
    end    
end        
% record data
results_struct.spot_array = spot_array;
results_struct.swap_array = swap_array;
results_struct.virtual_array = virtual_array;
results_struct.fluo_array = fluo_array;
results_struct.hmm_array = hmm_array;      
particle_id_vec(iter) = ParticleID;
results_struct.center_time_vec = center_time_vec;
results_struct.lag_dur_vec = lag_dur_vec;
results_struct.lead_dur_vec = lead_dur_vec;
results_struct.lag_size_vec = lag_size_vec;
results_struct.lead_size_vec = lead_size_vec;
results_struct.feature_sign_vec = feature_sign_vec;
% save
save([dataPath 'hmm_input_output_results.mat'],'results_struct')