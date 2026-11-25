function [outputs] = CHO_v27(signal, srate, param)
% CHO_v27() - Return detected bounding boxes of oscillations, 1/f fitting output,
%                and the std of the input signal
%
%           * Each bounding box has following attributes:
%               center_fp: center frequency
%               center_tp: middle time point of the bounding box
%               peak_val: peak_amplitude within the bounding box
%               min_F: lower range of the detected frequency
%               max_F: upper range of the detected frequency
%               start: onset of the oscillation
%               stop: off of the oscillation
%
%           * 1/f fitting output has following attributes:
%               init_ap_fit: 1/f fitted data points
%               initial_ap_fit_coef: fitting coefficients [b(1) b(2) b(3)] for model
%               function as below:
%               modelfun = @(b,x) b(1) - log10(b(2) + x(:,1).^b(3));
%
% Usage:
%       >>  param.minimum_cycles = 2;
%       >>  param.ovlp_threshold = 0.5;
%       >>  param.frequency_range = 1:60;  % optional, default 1:40 Hz
%       >>  cho_output = CHO_v27(data, srate, param);
%
% NOTE:
%       * none
%
% Required inputs:
%       data        = Single-channel data vector (required)
%       srate      = Sampling rate
%       param.plot = [0 or 1] if the value is 1, CHO plots raw
%                   signal with detected bounding boxes, a time-frequency map with
%                   ininital bounding boxes, a time-frequency map with final bounding
%                   boxes, and the power spectrum with 1/f fitting.
%       param.minimum_cycles = [a double value] if a bounding box has
%                   smaller number of cycles than 'minimum_cycles', CHO reject the
%                   bounding box.
%       param.ovlp_threshold = [0 < a double value < 1] if two bounding
%                   boxes overlap more than 'ovlp_threshold' in time domain and the center
%                   frequencies are next each other (difference <=1) CHO merges
%                   those two bounding boxes.
%
%
%
% Known problems:
%
%
%

% 06-03-23 v24 released
% 01-22-25 v25: wavelet v6 envelope and correlation-based overlap merge
% 06-14-26 v26: safer defaults, frequency-range handling, NaN defense,
%                Nyquist clipping, and richer output metadata.
% 06-15-26 v27: high-frequency ACF validation uses period tolerance instead
%                of replacing center frequency with integer-lag ACF frequency.
%                This avoids artificial sparse peaks such as 85/102/128/171 Hz
%                at 512 Hz sampling.
% 06-16-26 v27 update: detection wavelet cycles can be frequency-dependent.
%                Default is 'linear_3_10', matching the birdsong TF map:
%                low frequencies use 3 cycles and high frequencies use 10 cycles.
outputs = [];

if nargin < 3 || isempty(param)
    param = struct();
end
param = cho_set_default(param, 'plot', 0);
param = cho_set_default(param, 'three_d_plot', 0);
param = cho_set_default(param, 'acf_plots', 0);
param = cho_set_default(param, 'AC_PK_Threshold', 1);
param = cho_set_default(param, 'conceding_threshold_for_peaks', 10);
param = cho_set_default(param, 'minimum_cycles', 2);
param = cho_set_default(param, 'ovlp_threshold', 0.50);
param = cho_set_default(param, 'frequency_range', 1:40);
param = cho_set_default(param, 'frequency_step', 1);
param = cho_set_default(param, 'wavelet_cycles', 'linear_3_10');
param = cho_set_default(param, 'osc_threshold_std', 2);
param = cho_set_default(param, 'n_ap_fit_chunks', 4);
param = cho_set_default(param, 'acf_validation_mode', 'frequency_bin');
param = cho_set_default(param, 'acf_period_tolerance_pct', 25);
param = cho_set_default(param, 'acf_preserve_wavelet_center', 0);
param = cho_set_default(param, 'acf_lag_interpolation', 1);
param = cho_set_default(param, 'merge_mode', 'correlation'); % correlation | fast_tf_overlap | none
param = cho_set_default(param, 'max_boxes_for_correlation_merge', 250);
param = cho_set_default(param, 'fast_merge_time_overlap_frac', 0.05);
param = cho_set_default(param, 'fast_merge_freq_gap_hz', max(1, param.frequency_step));
param = cho_set_default(param, 'return_maps', 0);
param = cho_set_default(param, 'birdsong_tfmap', 0);
param = cho_set_default(param, 'plot_birdsong', 0);
param = cho_set_default(param, 'birdsong_freq_range', []);
param = cho_set_default(param, 'birdsong_n_freq', 80);
param = cho_set_default(param, 'birdsong_baseline_idx', []);
param = cho_set_default(param, 'birdsong_clim', [-6 12]);
param = cho_set_default(param, 'birdsong_colormap', 'turbo');

signal = double(signal(:));
if isempty(signal)
    outputs = make_empty_outputs(param, srate, 'empty_signal');
    return;
end
if any(~isfinite(signal))
    signal = fillmissing(signal, 'linear', 'EndValues', 'nearest');
end
if any(~isfinite(signal)) || std(signal) < eps
    outputs = make_empty_outputs(param, srate, 'invalid_or_flat_signal');
    return;
end

minimum_cycles = param.minimum_cycles;
AC_PK_Threshold = param.AC_PK_Threshold;
frequency_range = param.frequency_range;

%%
timex = 1:length(signal);
f_steps = param.frequency_step;
freqs = normalize_frequency_range(frequency_range, f_steps, srate);
if isempty(freqs)
    outputs = make_empty_outputs(param, srate, 'empty_frequency_range_after_nyquist_clip');
    return;
end
for ifreq = 1:length(freqs)
    freq = freqs(ifreq);
%         [ filtered_data(:,ifreq) ] = fft_bandpass_filtering(signal, srate, [freq-f_steps freq+f_steps]);
    this_wavelet_cycles = get_wavelet_cycles_for_frequency(freqs, ifreq, param);
    [~, filtered_data(:,ifreq)] = wavelet_single_freq_filtering_v6(signal, srate, freq, this_wavelet_cycles);
    env_data(:,ifreq) = log10(abs(filtered_data(:,ifreq)));
    med_env(ifreq) = median(env_data(:,ifreq));
end

param.freqs = freqs;
param.wavelet_cycles_by_freq = arrayfun(@(ii) get_wavelet_cycles_for_frequency(freqs, ii, param), 1:length(freqs));

%% select lowest one over f fit
clear med_psd
q4_points = round(linspace(1,timex(end),param.n_ap_fit_chunks+1),0);
q4_start = q4_points(1:end-1);
q4_end = q4_points(2:end);
X = freqs';
for iq = 1:length(q4_start)
    med_psd(:,iq) = median(env_data(q4_start(iq):q4_end(iq),:),1);
    Y = med_psd(:,iq);
    
    [oof_outputs] = oof_fitting_v2(Y,X);
    ap_fits(:,iq) = oof_outputs.init_ap_fit;
    %     kai_coef(:,iq) = oof_outputs.initial_ap_fit_coef(3);
    %     oof_outputs
end

% m_ap_fits = median(ap_fits);
[minval minidx] = min(mean(ap_fits,1));
% [max_kai_val max_kai_idx] = max(kai_coef);

%%
% figure();
% plot(ap_fits);
% figure();
% plot(med_psd);

%% remove pink noise
ap_fit = ap_fits(:,minidx);
ap_mat = repmat(ap_fit,1,length(timex));
denoised_env_data = env_data - ap_mat';
% denoised_env_data = denoised_env_data';
std_residue = std(denoised_env_data(:),'omitnan');

%% 3D ap fit
if param.three_d_plot == 1
    %%
    timesec = linspace(0,1.5,length(timex));
    f=figure();
    B = repmat(ap_fit,1,length(timesec));
    ap_s = surf(timesec,freqs,B);
    ap_s.EdgeColor = 'none';
    colormap bone;
    %         colormap(f, flipud(colormap(f)))
    %         colorbar
    hold on;
    freezeColors
    
    tf_s = surf(timesec,freqs,(env_data)');
    colormap parula;
    tf_s.EdgeColor = 'none';
    set(gca, 'xDir','normal');
    set(gca, 'yDir','reverse');
    xlabel('time');
    ylabel('frequency');
    zlabel('amplitude');
    xlim([0 1.5]);
    ylim([freqs(1) freqs(end)]);
    zlim([min(ap_fit)-.5 max(ap_fit)+1]);
    caxis([min(ap_fit)-.5 max(ap_fit)+.5]);
    
    set(gca,'fontsize',20);
    view(45,30);
    set(gcf,'color','w');
    
    %%
    f=figure();
    tf_s = surf(timesec,freqs,(denoised_env_data)');
    colormap parula;
    tf_s.EdgeColor = 'none';
    set(gca, 'xDir','normal');
    set(gca, 'yDir','reverse');
    xlabel('time');
    ylabel('frequency');
    zlabel('amplitude');
    xlim([0 1.5]);
    ylim([freqs(1) freqs(end)]);
    zlim([min(ap_fit)-.5 max(ap_fit)+1]);
    caxis([min(ap_fit)-.5 max(ap_fit)+.5]);
    
    set(gca,'fontsize',20);
    view(45,30);
    set(gcf,'color','w');
    
end

%% masking osc data
osc_threshold = param.osc_threshold_std * std_residue;
masked_env_data = zeros(size(env_data));
osc_data = zeros(size(denoised_env_data));
oevent_map = zeros(size(denoised_env_data));
for ifp = 1:length(freqs)
    sig_idxs = find(denoised_env_data(:,ifp)>osc_threshold);
    oevent_map(sig_idxs,ifp) = 1;
    osc_data(sig_idxs,ifp) = filtered_data(sig_idxs,ifp) .* oevent_map(sig_idxs,ifp);
    masked_env_data(sig_idxs,ifp) = denoised_env_data(sig_idxs,ifp) .* oevent_map(sig_idxs,ifp);
end

%%
% figure();
% imagesc(masked_env_data);

%%
freq_pnts = find(sum(masked_env_data,1)>0);
n_freq_pnts = length(freq_pnts);
bounding_boxes = [];
cycles_boxes = [];
cnt = 1;
% noise_rssq = rssq((filtered_data(:,1)));
% noise_amp = 10^max(ap_fits(:,iq));
noise_signal = signal;
for ifreq = 1:n_freq_pnts
    freq_data = masked_env_data(:,freq_pnts(ifreq));
    center_freq = freqs(freq_pnts(ifreq));
    ccl_data = bwlabel(freq_data);
    n_clusters = length(setdiff(unique(ccl_data),0));
    for icl = 1:n_clusters
        cl_idxs = find(ccl_data == icl);
        bounding_boxes(cnt).center_fp = center_freq;
        bounding_boxes(cnt).center_tp = round((cl_idxs(1)+cl_idxs(end))/2);
        bounding_boxes(cnt).peak_val = max(freq_data(cl_idxs(1):cl_idxs(end)));
        bounding_boxes(cnt).minF = max(freqs(1), center_freq-f_steps);
        bounding_boxes(cnt).maxF = min(freqs(end), center_freq+f_steps);
        bounding_boxes(cnt).start = cl_idxs(1);
        bounding_boxes(cnt).stop = cl_idxs(end);
        
        cycles_boxes(cnt) = center_freq * (bounding_boxes(cnt).stop - bounding_boxes(cnt).start)/srate;
        cnt = cnt + 1;
    end
end

initial_bounding_boxes = bounding_boxes;

%% remove bounding boxes that have less than 2 cycle

final_cycles_boxes = cycles_boxes;
bounding_boxes(cycles_boxes < minimum_cycles) = [];
final_cycles_boxes(cycles_boxes < minimum_cycles) = [];

bounding_boxes_second_step = bounding_boxes;

%% remove non autocorrelated boxes
n_boxes = length(bounding_boxes);
pxx_judge = [];
auto_corr_judge = [];
for ibox = 1:n_boxes
    
    time_window = bounding_boxes(ibox).start:bounding_boxes(ibox).stop;
    bounded_signal = signal(time_window);
    bounded_signal = detrend(bounded_signal);
    %     %%
    %     figure();
    %     plot(bounded_signal);
    %     hold on;
    %     plot(detrend(bounded_signal));
    
    
    [acf,lags,bounds] = autocorr(bounded_signal,'numlags',length(time_window)-1,'NumSTD',AC_PK_Threshold);
    %     [pxx_acf, acf_f] = fft_spectrum(acf,srate);
    
    %%
    pos_center_lags = [];
    %     neg_center_lags = [];
    
    pos_data = zeros(size(acf));
    pos_idxs = find(acf>bounds(1));
    pos_data(pos_idxs) = 1;
    pos_cluster_data = bwlabel(pos_data);
    n_clusters = length(setdiff(unique(pos_cluster_data),0));
    
    
    for icl = 1:n_clusters
        if icl == 1
            pos_center_lags(icl) = lags(1);
        else
            cl_idxs = find(pos_cluster_data == icl);
            [~, pk_rel_idx] = max(acf(cl_idxs));
            peak_idx = cl_idxs(pk_rel_idx);
            if param.acf_lag_interpolation
                pos_center_lags(icl) = interpolate_acf_peak_lag(lags, acf, peak_idx);
            else
                pos_center_lags(icl) = lags(peak_idx);
            end
        end
    end
    
    %     neg_data = zeros(size(acf));
    %     neg_idxs = find(acf<bounds(2));
    %     neg_data(neg_idxs) = 1;
    %     neg_cluster_data = bwlabel(neg_data);
    %     n_clusters = length(setdiff(unique(neg_cluster_data),0));
    %     for icl = 1:n_clusters
    %
    %         cl_idxs = find(neg_cluster_data == icl);
    %         center_idx = round((cl_idxs(1) + cl_idxs(end))/2);
    %         neg_center_lags(icl) = lags(center_idx);
    %
    %     end
    %
    if ~isempty(pos_center_lags)
        pos_center_lags(1) = [];
    else
        pos_center_lags = [];
    end
    
    
    %%
    pos_ac_LOCS = pos_center_lags;
    neg_ac_LOCS = [];
    %     neg_ac_LOCS = neg_center_lags;
    
    if length(pos_ac_LOCS)>1  % cycle is bigger than 2
        pos_invterval_simliarity_btw_pos_pks = abs(pos_ac_LOCS(1) - pos_ac_LOCS(2))/pos_ac_LOCS(1)*100;
        %         neg_invterval_simliarity_btw_pos_pks = abs(neg_ac_LOCS(1) - neg_ac_LOCS(2))/neg_ac_LOCS(1)*100;
        pos_herz = srate/pos_ac_LOCS(1);
        acf_freq_hz = srate/pos_ac_LOCS(1);
        expected_period_samples = srate / bounding_boxes(ibox).center_fp;
        acf_period_error_pct = abs(pos_ac_LOCS(1) - expected_period_samples) / expected_period_samples * 100;
        %         neg_herz = round(srate/neg_ac_LOCS(1));
        
        
        %         if abs(bounding_boxes(ibox).center_fp - herz) < 1
        if pos_herz >= freqs(1) && pos_herz <= freqs(end)
            bounding_boxes(ibox).acf_fp = acf_freq_hz;
            bounding_boxes(ibox).acf_period_error_pct = acf_period_error_pct;
            bounding_boxes(ibox).acf_interval_similarity_pct = pos_invterval_simliarity_btw_pos_pks;

            if strcmpi(param.acf_validation_mode, 'period_tolerance')
                acf_freq_ok = acf_period_error_pct <= param.acf_period_tolerance_pct;
            else
                acf_freq_ok = bounding_boxes(ibox).minF <= pos_herz && bounding_boxes(ibox).maxF >= pos_herz;
            end

            if acf_freq_ok && ...
                    ( abs(pos_invterval_simliarity_btw_pos_pks - 100) < param.conceding_threshold_for_peaks )
                
                if ~param.acf_preserve_wavelet_center
                    bounding_boxes(ibox).center_fp = pos_herz;
                    nearest_freq_idx = nearest_frequency_index(freqs, pos_herz);
                    bounding_boxes(ibox).peak_val = max(masked_env_data( ...
                        bounding_boxes(ibox).start:bounding_boxes(ibox).stop, nearest_freq_idx));
                end
                
                auto_corr_judge(ibox) = 1;
            else
                auto_corr_judge(ibox) = 0;
            end
        else
            auto_corr_judge(ibox) = 0;
        end
    else
        pos_herz = 0;
        acf_freq_hz = NaN;
        acf_period_error_pct = NaN;
        pos_ac_LOCS = 0;
        neg_ac_LOCS = [];
        auto_corr_judge(ibox) = 0;
        pos_invterval_simliarity_btw_pos_pks = 0;
        neg_invterval_simliarity_btw_pos_pks = 0;
        bounding_boxes(ibox).acf_fp = acf_freq_hz;
        bounding_boxes(ibox).acf_period_error_pct = acf_period_error_pct;
        bounding_boxes(ibox).acf_interval_similarity_pct = pos_invterval_simliarity_btw_pos_pks;
    end
    
    %%
    if param.acf_plots == 1
        n_plots = 2;
        figure();
        subplot(n_plots,1,1);
        plot(signal);
        vline([time_window(1) time_window(end)]);
        
        subplot(n_plots,1,2);
        plot(lags,acf);
        hline(bounds);
        vline(pos_ac_LOCS);
        vline(neg_ac_LOCS);
        text(0,0.5,[num2str(pos_herz) 'Hz detected'],'color','red');
        text(0,0.3,['Frequency bounds: ' num2str(bounding_boxes(ibox).minF) '-' num2str(bounding_boxes(ibox).maxF) 'Hz'],'color','red');
        text(0,0.1,['Similarity btw 2 peaks ' num2str(pos_invterval_simliarity_btw_pos_pks) '%'],'color','red');
        text(0,-0.1,['Neg. Similarity btw 2 peaks ' num2str(neg_invterval_simliarity_btw_pos_pks) '%'],'color','red');
        
        if auto_corr_judge(ibox)
            title('True oscillation!');
            
        else
            title('False oscillation!');
        end
    end
end
before_auto_corr_check = bounding_boxes;
auto_corr_judge = logical(auto_corr_judge);
bounding_boxes(~auto_corr_judge) = [];
final_cycles_boxes(~auto_corr_judge) = [];
% bounding_boxes(~(auto_corr_judge & pxx_judge)) = [];
% final_cycles_boxes(~(auto_corr_judge & pxx_judge)) = [];

%% merge overlapping boxes
n_boxes = length(bounding_boxes);
end_criteria = 0;
final_bounding_boxes = bounding_boxes;
% find overlapping boxes
overlap_list = [];
cnt = 1;
box_clusters = [];
for ibox = 1:n_boxes
    box_clusters(ibox).member = [ibox];
end

if strcmpi(param.merge_mode, 'none') || isempty(final_bounding_boxes)
    cluster_label_for_boxes = 1:length(bounding_boxes);
elseif strcmpi(param.merge_mode, 'fast_tf_overlap') || n_boxes > param.max_boxes_for_correlation_merge
    [final_bounding_boxes, cluster_label_for_boxes] = merge_boxes_fast_tf_overlap(bounding_boxes, param);
else
    while(~end_criteria)
        n_boxes = length(final_bounding_boxes);
        inner_flag = 0;

        for ibox = 1:n_boxes
            for jbox = 1:n_boxes
                if jbox > ibox
                    %                 if is_boxes_overlapp(final_bounding_boxes(ibox),final_bounding_boxes(jbox),ovlp_threshold)
                    if is_boxes_overlapped_signal(final_bounding_boxes(ibox),final_bounding_boxes(jbox),signal,srate,param)

                        overlap_list(cnt,:) = [ibox,jbox];

                        new_box = merge_two_boxes(final_bounding_boxes(ibox), final_bounding_boxes(jbox));

                        final_bounding_boxes(ibox) = new_box;
                        final_bounding_boxes(jbox) = [];

                        box_clusters(ibox).member = unique([box_clusters(ibox).member box_clusters(jbox).member]);
                        box_clusters(jbox) = [];

                        inner_flag = 1;
                        cnt = cnt + 1;
                        break;
                    else

                    end
                end
            end
            if inner_flag == 1
                break;
            end
        end

        if inner_flag == 0
            end_criteria = 1;
        end

    end

    n_clusters = length(box_clusters);
    cluster_label_for_boxes = 1:length(bounding_boxes);
    for icls = 1:n_clusters
        cluster_label_for_boxes(box_clusters(icls).member) = icls;
    end
end


%% visualization
if param.plot > 0
    n_tp = length(signal);
    
    linewidth_val = 2;
    timex = 1:n_tp;
    n_plots = 4;
    figure('position',[10 10 400 1000]);
    subplot(n_plots,1,1);
    plot(timex,(signal));
    hold on;
    %     plot(timex,zscore(hg_env_signal));
    p_bounding_boxes1 = final_bounding_boxes;
    for ibox = 1:length(p_bounding_boxes1)
        %     figure();
        min_val = min(signal(p_bounding_boxes1(ibox).start:p_bounding_boxes1(ibox).stop));
        max_val = max(signal(p_bounding_boxes1(ibox).start:p_bounding_boxes1(ibox).stop));
        line([p_bounding_boxes1(ibox).start p_bounding_boxes1(ibox).stop],[max_val max_val],'color','red');
        line([p_bounding_boxes1(ibox).start p_bounding_boxes1(ibox).start],[min_val max_val],'color','red');
        line([p_bounding_boxes1(ibox).stop p_bounding_boxes1(ibox).stop],[min_val max_val],'color','red');
        line([p_bounding_boxes1(ibox).start p_bounding_boxes1(ibox).stop],[min_val min_val],'color','red');
        text((p_bounding_boxes1(ibox).start+p_bounding_boxes1(ibox).stop)/2,0, [num2str(p_bounding_boxes1(ibox).center_fp) 'Hz'],'color','red');
        fprintf('\n%dHz\n', p_bounding_boxes1(ibox).center_fp);
    end
    xsteps = n_tp/10;
    xticks = 0:xsteps:n_tp;
    xlim([1 n_tp])
    set(gca,'xtick',xticks);
    set(gca,'xticklabel',round(xticks/srate,2));
    %     set(gca,'fontsize',20);
    %     set(gca,'linewidth',linewidth_val);
    xlabel('time');
    %     ylim([-200 200]);
    
    subplot(n_plots,1,2);
    if isnan(std_residue)
        imagesc(timex,freqs,denoised_env_data',[0 1]);
    else
        imagesc(timex,freqs,denoised_env_data',[0 3*std_residue]);
    end
    % imagesc(timex,freqs,BW');
    % imagesc(timex,freqs,ovent_map');
    set(gca, 'YDir','normal');
    ylim([max(0.5, freqs(1)-0.5) freqs(end)]);
    %     hline([3 8 14 30]);
    hold on;
    p_bounding_boxes2 = initial_bounding_boxes;
    
    if length(p_bounding_boxes2) < 20
        for ibox = 1:length(p_bounding_boxes2)
            plot(p_bounding_boxes2(ibox).center_tp,p_bounding_boxes2(ibox).center_fp,'r.','markersize',10);
            hold on;
            line([p_bounding_boxes2(ibox).start p_bounding_boxes2(ibox).stop],[p_bounding_boxes2(ibox).maxF p_bounding_boxes2(ibox).maxF],'color','red');
            line([p_bounding_boxes2(ibox).start p_bounding_boxes2(ibox).start],[p_bounding_boxes2(ibox).minF p_bounding_boxes2(ibox).maxF],'color','red');
            line([p_bounding_boxes2(ibox).stop p_bounding_boxes2(ibox).stop],[p_bounding_boxes2(ibox).minF p_bounding_boxes2(ibox).maxF],'color','red');
            line([p_bounding_boxes2(ibox).start p_bounding_boxes2(ibox).stop],[p_bounding_boxes2(ibox).minF p_bounding_boxes2(ibox).minF],'color','red');
        end
    else
        text(10,freqs(end/2),'More than 20 bounding boxes are initially generated.','color','red');
    end
    
    
    %     set(gca,'fontsize',20);
    %     set(gca,'linewidth',linewidth_val);
    xlabel('time');
    %         xticks = 0:100:n_tp;
    set(gca,'xtick',xticks);
    set(gca,'xticklabel',round(xticks/srate,2));
    
    
    subplot(n_plots,1,3);
    %     imagesc(timex,freqs,osc_data');
%     imagesc(timex,freqs,denoised_env_data',[0 3*std_residue]);

    if isnan(std_residue)
        imagesc(timex,freqs,denoised_env_data');
    else
        imagesc(timex,freqs,denoised_env_data',[0 3*std_residue]);
    end

    % imagesc(timex,freqs,BW');
    set(gca, 'YDir','normal');
    ylim([max(0.5, freqs(1)-0.5) freqs(end)]);
    %     hline([3 8 14 30]);
    hold on;
    p_bounding_boxes3 = final_bounding_boxes;
    for ibox = 1:length(p_bounding_boxes3)
        %     figure();
        plot(p_bounding_boxes3(ibox).center_tp,p_bounding_boxes3(ibox).center_fp,'r.','markersize',10);
        hold on;
        line([p_bounding_boxes3(ibox).start p_bounding_boxes3(ibox).stop],[p_bounding_boxes3(ibox).maxF p_bounding_boxes3(ibox).maxF],'color','red');
        line([p_bounding_boxes3(ibox).start p_bounding_boxes3(ibox).start],[p_bounding_boxes3(ibox).minF p_bounding_boxes3(ibox).maxF],'color','red');
        line([p_bounding_boxes3(ibox).stop p_bounding_boxes3(ibox).stop],[p_bounding_boxes3(ibox).minF p_bounding_boxes3(ibox).maxF],'color','red');
        line([p_bounding_boxes3(ibox).start p_bounding_boxes3(ibox).stop],[p_bounding_boxes3(ibox).minF p_bounding_boxes3(ibox).minF],'color','red');
    end
    
    %     set(gca,'fontsize',20);
    %     set(gca,'linewidth',linewidth_val);
    xlabel('time');
    set(gca,'xtick',xticks);
    set(gca,'xticklabel',round(xticks/srate,2));
    
    subplot(n_plots,1,4);
    plot(freqs,med_env);
    hold on;
    plot(freqs,ap_fit);
    xlabel('Hz');
    
end


%%
% eegplot([signal osc_comps]','srate',srate,'winlength',1);

%%
outputs = oof_outputs;
outputs.bounding_boxes = final_bounding_boxes;
outputs.signal_std = std(signal);
outputs.freqs = freqs;
outputs.param = param;
outputs.status = "ok";
outputs.n_initial_boxes = length(initial_bounding_boxes);
outputs.n_after_cycle_check_boxes = length(bounding_boxes_second_step);
outputs.n_before_autocorr_check_boxes = length(before_auto_corr_check);
outputs.n_after_autocorr_check_boxes = length(bounding_boxes);
outputs.n_final_boxes = length(final_bounding_boxes);
outputs.initial_bounding_boxes = initial_bounding_boxes;
outputs.before_auto_corr_check = before_auto_corr_check;
outputs.bounding_boxes_second_step = bounding_boxes_second_step;
outputs.ap_fit = ap_fit;
outputs.med_env = med_env;
outputs.std_residue = std_residue;
outputs.osc_threshold = osc_threshold;
if param.return_maps
    outputs.filtered_data = filtered_data;
    outputs.env_data = env_data;
    outputs.denoised_env_data = denoised_env_data;
    outputs.masked_env_data = masked_env_data;
    outputs.oevent_map = oevent_map;
end
if param.birdsong_tfmap || param.plot_birdsong
    if isempty(param.birdsong_freq_range)
        param.birdsong_freq_range = [freqs(1) freqs(end)];
    end
    TF = cho_birdsong_tfmap(signal, srate, param.birdsong_freq_range, ...
        param.birdsong_n_freq, param.birdsong_baseline_idx);
    outputs.birdsong_tfmap = TF;
    if param.plot_birdsong
        outputs.birdsong_figure = plot_cho_birdsong_tfmap(TF, final_bounding_boxes, ...
            srate, param);
    end
end

% outputs.osc_comps = osc_comps;

end


function [merged_boxes, cluster_label_for_boxes] = merge_boxes_fast_tf_overlap(boxes, param)
if isempty(boxes)
    merged_boxes = boxes;
    cluster_label_for_boxes = [];
    return;
end

n = numel(boxes);
parent = 1:n;
[~, order] = sortrows([[boxes.start]' [boxes.center_fp]']);

for oi = 1:n
    i = order(oi);
    for oj = oi+1:n
        j = order(oj);
        if boxes(j).start > boxes(i).stop
            break;
        end
        if should_fast_merge_boxes(boxes(i), boxes(j), param)
            parent = union_parent(parent, i, j);
        end
    end
end

for i = 1:n
    parent(i) = find_parent(parent, i);
end
[roots, ~, cluster_label_for_boxes] = unique(parent);
merged_boxes = boxes(roots);
for k = 1:numel(roots)
    members = find(cluster_label_for_boxes == k);
    merged_boxes(k) = merge_box_members(boxes(members));
end
end

function tf = should_fast_merge_boxes(box1, box2, param)
time_overlap = min(box1.stop, box2.stop) - max(box1.start, box2.start) + 1;
if time_overlap <= 0
    tf = false;
    return;
end
dur1 = max(1, box1.stop - box1.start + 1);
dur2 = max(1, box2.stop - box2.start + 1);
time_frac = time_overlap / min(dur1, dur2);

freq_gap = max([box1.minF box2.minF]) - min([box1.maxF box2.maxF]);
freq_touch = freq_gap <= param.fast_merge_freq_gap_hz;
center_close = abs(box1.center_fp - box2.center_fp) <= max(param.fast_merge_freq_gap_hz, param.frequency_step);

tf = time_frac >= param.fast_merge_time_overlap_frac && (freq_touch || center_close);
end

function box = merge_box_members(member_boxes)
box = member_boxes(1);
box.start = min([member_boxes.start]);
box.stop = max([member_boxes.stop]);
box.minF = min([member_boxes.minF]);
box.maxF = max([member_boxes.maxF]);
[box.peak_val, idx] = max([member_boxes.peak_val]);
box.center_fp = member_boxes(idx).center_fp;
box.center_tp = member_boxes(idx).center_tp;
if isfield(member_boxes, 'acf_fp')
    acf_vals = [member_boxes.acf_fp];
    box.acf_fp = acf_vals(idx);
end
if isfield(member_boxes, 'acf_period_error_pct')
    err_vals = [member_boxes.acf_period_error_pct];
    box.acf_period_error_pct = err_vals(idx);
end
if isfield(member_boxes, 'acf_interval_similarity_pct')
    sim_vals = [member_boxes.acf_interval_similarity_pct];
    box.acf_interval_similarity_pct = sim_vals(idx);
end
end

function box = merge_two_boxes(box1, box2)
box = box1;
box.start = min(box1.start, box2.start);
box.stop = max(box1.stop, box2.stop);
box.minF = min(box1.minF, box2.minF);
box.maxF = max(box1.maxF, box2.maxF);
two_fp_candidates = [box1.center_fp box2.center_fp];
two_tp_candidates = [box1.center_tp box2.center_tp];
two_pkvals = [box1.peak_val box2.peak_val];
[box.peak_val, maxi] = max(two_pkvals);
box.center_fp = two_fp_candidates(maxi);
box.center_tp = two_tp_candidates(maxi);
if isfield(box1, 'acf_fp') && isfield(box2, 'acf_fp')
    acf_vals = [box1.acf_fp box2.acf_fp];
    box.acf_fp = acf_vals(maxi);
end
if isfield(box1, 'acf_period_error_pct') && isfield(box2, 'acf_period_error_pct')
    err_vals = [box1.acf_period_error_pct box2.acf_period_error_pct];
    box.acf_period_error_pct = err_vals(maxi);
end
if isfield(box1, 'acf_interval_similarity_pct') && isfield(box2, 'acf_interval_similarity_pct')
    sim_vals = [box1.acf_interval_similarity_pct box2.acf_interval_similarity_pct];
    box.acf_interval_similarity_pct = sim_vals(maxi);
end
end

function parent = union_parent(parent, a, b)
ra = find_parent(parent, a);
rb = find_parent(parent, b);
if ra ~= rb
    parent(rb) = ra;
end
end

function r = find_parent(parent, a)
r = a;
while parent(r) ~= r
    r = parent(r);
end
end


function [outputs] = is_boxes_overlapped_signal(box1,box2,signal,srate,param)
vol1 = (box1.maxF-box1.minF) * (box1.stop-box1.start);
vol2 = (box2.maxF-box2.minF) * (box2.stop-box2.start);
%%
overlapping_threshold = param.ovlp_threshold;
%%
time1 = box1.start:box1.stop;
time2 = box2.start:box2.stop;
freq1 = box1.minF:box1.maxF;
freq2 = box2.minF:box2.maxF;
center_freq1 = box1.center_fp;
center_freq2 = box2.center_fp;
[time_intersect] = intersect(time1,time2);
[freq_intersect] = intersect(freq1,freq2);
vol_intersect = length(time_intersect)*length(freq_intersect);

% if abs(center_freq1 - center_freq2) > 2
%     outputs = 0;
%     return;
% end

if length(time_intersect) == 0
    outputs = 0;
    return;
end

%%
% overlapped_data1 = fft_bandpass_filtering(signal,srate,[box1.minF box1.maxF]);
cycles1 = get_wavelet_cycles_for_center_frequency(center_freq1, param);
cycles2 = get_wavelet_cycles_for_center_frequency(center_freq2, param);
overlapped_data1 = wavelet_single_freq_filtering_v4(signal,srate,center_freq1,cycles1);
overlapped_data2 = wavelet_single_freq_filtering_v4(signal,srate,center_freq2,cycles2);
windowing_mask = tukeywin(length(signal),0.25);
overlapped_data1 = overlapped_data1 .* windowing_mask;
overlapped_data2 = overlapped_data2 .* windowing_mask;
% angle_data1 = angle(hilbert(overlapped_data1));
% angle_data2 = angle(hilbert(overlapped_data2));
% cos_data1 = cos(angle_data1);
% cos_data2 = cos(angle_data2);

[coef, pval] = corr(overlapped_data1(time_intersect(1):time_intersect(end)),overlapped_data2(time_intersect(1):time_intersect(end)));


% if coef > overlapping_threshold || (vol_intersect>vol1*overlapping_threshold) || (vol_intersect>vol2*overlapping_threshold)
if pval < 0.05 && coef > 0 || (vol_intersect>vol1*overlapping_threshold) || (vol_intersect>vol2*overlapping_threshold)
    outputs = 1;
else
    outputs = 0;
end

end

% 
% function [outputs] = is_boxes_overlapped_signal(box1,box2,signal,srate,param)
% % vol1 = (box1.maxF-box1.minF) * (box1.stop-box1.start);
% % vol2 = (box2.maxF-box2.minF) * (box2.stop-box2.start);
% 
% vol1 = (box1.stop-box1.start);
% vol2 = (box2.stop-box2.start);
% %%
% overlapping_threshold = param.ovlp_threshold;
% %%
% time1 = box1.start:box1.stop;
% time2 = box2.start:box2.stop;
% freq1 = box1.minF:box1.maxF;
% freq2 = box2.minF:box2.maxF;
% center_freq1 = box1.center_fp;
% center_freq2 = box2.center_fp;
% [time_intersect] = intersect(time1,time2);
% % [freq_intersect] = intersect(freq1,freq2);
% % vol_intersect = length(time_intersect)*length(freq_intersect);
% vol_intersect = length(time_intersect);
% 
% 
% if max([center_freq1, center_freq2])/min([center_freq1, center_freq2]) > 2
%     outputs = 0;
%     return;
% end
% 
% % if abs(center_freq1 - center_freq2) > 2
% %     outputs = 0;
% %     return;
% % end
% 
% % if length(time_intersect) < 2
% %     outputs = 0;
% % else
% %     if vol_intersect/vol1 > overlapping_threshold || vol_intersect/vol2 > overlapping_threshold
% %         outputs = 1;
% %     else
% %         outputs = 0;
% %     end    
% % end
% 
% if length(time_intersect) < 2
%     outputs = 0;
% else
%     if vol_intersect/vol1 > overlapping_threshold || vol_intersect/vol2 > overlapping_threshold
%         outputs = 1;
%     else
%         outputs = 0;
%     end    
% end
% 
% 
% end

function param = cho_set_default(param, field_name, default_value)
if ~isfield(param, field_name) || isempty(param.(field_name))
    param.(field_name) = default_value;
end
end

function lag_est = interpolate_acf_peak_lag(lags, acf, peak_idx)
lag_est = lags(peak_idx);
if peak_idx <= 1 || peak_idx >= numel(acf)
    return;
end
y1 = acf(peak_idx - 1);
y2 = acf(peak_idx);
y3 = acf(peak_idx + 1);
den = y1 - 2*y2 + y3;
if ~isfinite(den) || abs(den) < eps
    return;
end
delta = 0.5 * (y1 - y3) / den;
if ~isfinite(delta) || abs(delta) > 1
    return;
end
lag_step = median(diff(lags), 'omitnan');
if ~isfinite(lag_step) || lag_step <= 0
    lag_step = 1;
end
lag_est = lags(peak_idx) + delta * lag_step;
end

function freqs = normalize_frequency_range(frequency_range, f_step, srate)
frequency_range = double(frequency_range(:)');
frequency_range = frequency_range(isfinite(frequency_range) & frequency_range > 0);
if isempty(frequency_range)
    freqs = [];
    return;
end
if numel(frequency_range) == 2 && diff(frequency_range) > f_step
    freqs = frequency_range(1):f_step:frequency_range(2);
else
    freqs = unique(frequency_range, 'stable');
end
nyq_guard = max(0, floor(srate/2) - f_step);
freqs = freqs(freqs >= f_step & freqs <= nyq_guard);
freqs = round(freqs, 6);
end

function cycles = get_wavelet_cycles_for_frequency(freqs, ifreq, param)
cycles_param = param.wavelet_cycles;
if isnumeric(cycles_param)
    cycles_param = double(cycles_param(:)');
    if isempty(cycles_param)
        cycles_vec = linspace(3, 10, numel(freqs));
        cycles = cycles_vec(ifreq);
    elseif isscalar(cycles_param)
        cycles = cycles_param;
    elseif numel(cycles_param) == numel(freqs)
        cycles = cycles_param(ifreq);
    else
        x = linspace(1, numel(freqs), numel(cycles_param));
        cycles = interp1(x, cycles_param, ifreq, 'linear', 'extrap');
    end
elseif ischar(cycles_param) || isstring(cycles_param)
    mode = lower(string(cycles_param));
    if mode == "linear" || mode == "linear_3_10" || mode == "birdsong"
        cycles_vec = linspace(3, 10, numel(freqs));
        cycles = cycles_vec(ifreq);
    else
        error('Unknown wavelet_cycles mode: %s', mode);
    end
else
    cycles_vec = linspace(3, 10, numel(freqs));
    cycles = cycles_vec(ifreq);
end
cycles = max(0.5, double(cycles));
end

function cycles = get_wavelet_cycles_for_center_frequency(center_freq, param)
if isfield(param, 'freqs') && isfield(param, 'wavelet_cycles_by_freq') && ...
        ~isempty(param.freqs) && ~isempty(param.wavelet_cycles_by_freq)
    [~, ix] = min(abs(param.freqs - center_freq));
    cycles = param.wavelet_cycles_by_freq(ix);
else
    cycles = get_wavelet_cycles_for_frequency(center_freq, 1, param);
end
end

function idx = nearest_frequency_index(freqs, target_freq)
[~, idx] = min(abs(freqs - target_freq));
end

function outputs = make_empty_outputs(param, srate, status)
outputs = struct();
outputs.status = string(status);
outputs.bounding_boxes = [];
outputs.signal_std = NaN;
outputs.freqs = [];
outputs.param = param;
outputs.srate = srate;
outputs.n_initial_boxes = 0;
outputs.n_after_cycle_check_boxes = 0;
outputs.n_before_autocorr_check_boxes = 0;
outputs.n_after_autocorr_check_boxes = 0;
outputs.n_final_boxes = 0;
outputs.initial_bounding_boxes = [];
outputs.before_auto_corr_check = [];
outputs.bounding_boxes_second_step = [];
outputs.ap_fit = [];
outputs.med_env = [];
outputs.std_residue = NaN;
outputs.osc_threshold = NaN;
end

function TF = cho_birdsong_tfmap(x, fs, freq_range, n_freq, baseline_idx)
% Birdsong-style log-frequency Morlet TF map, adapted from
% seizure_tfmap_birdsongStyle.m.
x = double(x(:)');
n_time = numel(x);
time = (0:n_time-1) ./ fs;

freq_range = double(freq_range(:)');
freq_range = freq_range(isfinite(freq_range) & freq_range > 0);
if isempty(freq_range)
    freq_range = [1 min(200, floor(fs/2)-1)];
elseif numel(freq_range) == 1
    freq_range = [freq_range freq_range];
else
    freq_range = [min(freq_range) max(freq_range)];
end
freq_range(2) = min(freq_range(2), floor(fs/2)-1);
freq_range(1) = max(freq_range(1), 0.5);
if freq_range(2) <= freq_range(1)
    freq_range(2) = min(floor(fs/2)-1, freq_range(1)+1);
end

n_freq = max(2, round(n_freq));
freqs_tf = logspace(log10(freq_range(1)), log10(freq_range(2)), n_freq);
n_cycles = linspace(3, 10, n_freq);

power = zeros(n_freq, n_time);
for fi = 1:n_freq
    f = freqs_tf(fi);
    nc = n_cycles(fi);
    t_wave = -2:1/fs:2;
    sigma = nc / (2*pi*f);
    wavelet = exp(2*1i*pi*f*t_wave) .* exp(-t_wave.^2 ./ (2*sigma^2));
    wavelet = wavelet ./ sqrt(sum(abs(wavelet).^2));
    conv_result = conv(x, wavelet, 'same');
    power(fi, :) = abs(conv_result).^2;
end

if isempty(baseline_idx)
    baseline_idx = 1:n_time;
else
    baseline_idx = baseline_idx(:)';
    baseline_idx = baseline_idx(isfinite(baseline_idx));
    baseline_idx = round(baseline_idx);
    baseline_idx = baseline_idx(baseline_idx >= 1 & baseline_idx <= n_time);
    if isempty(baseline_idx)
        baseline_idx = 1:n_time;
    end
end

baseline_power = mean(power(:, baseline_idx), 2, 'omitnan');
baseline_power(~isfinite(baseline_power) | baseline_power <= eps) = eps;
power_db = 10 .* log10(power ./ baseline_power);

TF = struct();
TF.power = power_db;
TF.freqs = freqs_tf;
TF.time = time;
TF.baseline_idx = baseline_idx;
TF.freq_range = freq_range;
TF.n_cycles = n_cycles;
end

function fig = plot_cho_birdsong_tfmap(TF, bounding_boxes, fs, param)
fig = figure('Color','w', 'Position',[100 100 1050 560]);
imagesc(TF.time, TF.freqs, TF.power);
set(gca, 'YDir','normal', 'YScale','log', 'FontSize', 12);
xlabel('Time (sec)');
ylabel('Frequency (Hz)');
title('CHO birdsong-style time-frequency map', 'FontSize', 16, 'FontWeight','bold');
try
    colormap(gca, feval(param.birdsong_colormap, 256));
catch
    colormap(gca, turbo(256));
end
colorbar;
if ~isempty(param.birdsong_clim) && numel(param.birdsong_clim) == 2
    clim(param.birdsong_clim);
end
hold on;
for ibox = 1:numel(bounding_boxes)
    x1 = (bounding_boxes(ibox).start - 1) ./ fs;
    x2 = (bounding_boxes(ibox).stop - 1) ./ fs;
    y1 = bounding_boxes(ibox).minF;
    y2 = bounding_boxes(ibox).maxF;
    rectangle('Position',[x1 y1 max(x2-x1, eps) max(y2-y1, eps)], ...
        'EdgeColor',[1 0.1 0.1], 'LineWidth', 1.5);
    text((x1+x2)/2, bounding_boxes(ibox).center_fp, ...
        sprintf('%g Hz', bounding_boxes(ibox).center_fp), ...
        'Color','w', 'FontSize', 9, 'FontWeight','bold', ...
        'HorizontalAlignment','center', 'VerticalAlignment','middle');
end
end
