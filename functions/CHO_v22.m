function [outputs] = CHO_v22(signal, srate, param)
% CHO_v16() - Return detected bounding boxes of oscillations, 1/f fitting output,
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
%       >>  [cho_output] = CHO_v21(data, srate, param);    
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
% Required toolboxes: 
%   Econometric Toolbox for autocorr()
%   
%

% 06-03-23 released
% 12-18-23 zero value issue in 'peak_val' of bounding box fixed (reported by Eleonora Marcantoni <2832601M@student.gla.ac.uk>) 
%


outputs = [];


param.three_d_plot = 0; % plotting time-frequency map --- debuging purpose
param.acf_plots = 0; % plotting auto-correlation results --- debuging purpose
param.AC_PK_Threshold = 1; % 'NumSTD' parameter for autocorr() in Matlab );
param.rounding = 1;

minimum_cycles = param.minimum_cycles;
AC_PK_Threshold = param.AC_PK_Threshold;
ovlp_threshold = param.ovlp_threshold;
frequency_vector = param.frequency_vector;

%%
timex = 1:length(signal);
freqs = frequency_vector;
f_steps = freqs(2)-freqs(1);
for ifreq = 1:length(freqs)
    freq = freqs(ifreq);
    [ filtered_data(:,ifreq) ] = fft_bandpass_filtering(signal, srate, [freq-f_steps freq+f_steps]);
    %     [ filtered_data(:,ifreq) ] = wavelet_single_freq_filtering_v2(signal, srate, ifreq, cycles);
    env_data(:,ifreq) = log10(abs(hilbert(filtered_data(:,ifreq))));
    med_env(ifreq) = median(env_data(:,ifreq));
end

%% select lowest one over f fit
clear med_psd
q4_points = round(linspace(1,timex(end),5),0);
q4_start = q4_points(1:end-1);
q4_end = q4_points(2:end);
X = freqs';
for iq = 1:length(q4_start)
    med_psd(:,iq) = median(env_data(q4_start(iq):q4_end(iq),:),1);   
    Y = med_psd(:,iq);
    
    [oof_outputs] = oof_fitting_v2(Y,X);
    ap_fits(:,iq) = oof_outputs.init_ap_fit;
end

% m_ap_fits = median(ap_fits);
[minval minidx] = min(ap_fits(1,:));


%% remove pink noise
ap_fit = ap_fits(:,minidx);
ap_mat = repmat(ap_fit,1,length(timex));
denoised_env_data = env_data - ap_mat';
% denoised_env_data = denoised_env_data';
std_residue = std(denoised_env_data(:));

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
        ylim([1 40]);
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
        ylim([1 40]);
        zlim([min(ap_fit)-.5 max(ap_fit)+1]);
        caxis([min(ap_fit)-.5 max(ap_fit)+.5]);
    
        set(gca,'fontsize',20);
        view(45,30);
        set(gcf,'color','w');
        
end

%% masking osc data
osc_threshold = 2*std_residue;
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
    ccl_data = bwlabel(freq_data);
    n_clusters = length(setdiff(unique(ccl_data),0));
    for icl = 1:n_clusters
        cl_idxs = find(ccl_data == icl);
        bounding_boxes(cnt).center_fp = freq_pnts(ifreq);
        bounding_boxes(cnt).center_tp = round((cl_idxs(1)+cl_idxs(end))/2);
        bounding_boxes(cnt).peak_val = max(freq_data(cl_idxs(1):cl_idxs(end)));
        bounding_boxes(cnt).minF = freq_pnts(ifreq)-1;
        bounding_boxes(cnt).maxF = freq_pnts(ifreq)+1;
        bounding_boxes(cnt).start = cl_idxs(1);
        bounding_boxes(cnt).stop = cl_idxs(end);
        
        cycles_boxes(cnt) = freq_pnts(ifreq) * (bounding_boxes(cnt).stop - bounding_boxes(cnt).start)/srate;
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
    [acf,lags,bounds] = autocorr(bounded_signal,'numlags',length(time_window)-1,'NumSTD',AC_PK_Threshold);
    [pxx_acf, acf_f] = fft_spectrum(acf,srate);
    
    %%
    pos_center_lags = [];
    neg_center_lags = [];
    
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
            center_idx = round((cl_idxs(1) + cl_idxs(end))/2);
            pos_center_lags(icl) = lags(center_idx);
        end
    end
    
    neg_data = zeros(size(acf));
    neg_idxs = find(acf<bounds(2));
    neg_data(neg_idxs) = 1;
    neg_cluster_data = bwlabel(neg_data);
    n_clusters = length(setdiff(unique(neg_cluster_data),0));
    for icl = 1:n_clusters
        
        cl_idxs = find(neg_cluster_data == icl);
        center_idx = round((cl_idxs(1) + cl_idxs(end))/2);
        neg_center_lags(icl) = lags(center_idx);
        
    end
    
    pos_center_lags(1) = [];
    
   
    %%
    pos_ac_LOCS = pos_center_lags;
    neg_ac_LOCS = neg_center_lags;
    
    if length(pos_ac_LOCS)>1 && length(neg_ac_LOCS)>1 % cycle is bigger than 2
        pos_invterval_simliarity_btw_pos_pks = abs(pos_ac_LOCS(1) - pos_ac_LOCS(2))/pos_ac_LOCS(1)*100;
        neg_invterval_simliarity_btw_pos_pks = abs(neg_ac_LOCS(1) - neg_ac_LOCS(2))/pos_ac_LOCS(1)*100;
        pos_herz = round(srate/pos_ac_LOCS(1));
        neg_herz = round(srate/neg_ac_LOCS(1));
        %         if abs(bounding_boxes(ibox).center_fp - herz) < 1
        if pos_herz >= freqs(1) && pos_herz <= freqs(end)
            if bounding_boxes(ibox).minF <= pos_herz && bounding_boxes(ibox).maxF >= pos_herz && ( abs(pos_invterval_simliarity_btw_pos_pks - 100) < 30 || abs(neg_invterval_simliarity_btw_pos_pks - 100) < 30 )
                bounding_boxes(ibox).center_fp = pos_herz;
                tmp = masked_env_data(bounding_boxes(ibox).start:bounding_boxes(ibox).stop,bounding_boxes(ibox).minF:bounding_boxes(ibox).maxF);
                bounding_boxes(ibox).peak_val = max(tmp(:));
                                
                auto_corr_judge(ibox) = 1;
            else
                auto_corr_judge(ibox) = 0;
            end
        else
            auto_corr_judge(ibox) = 0;
        end
    else
        pos_herz = 0;
        pos_ac_LOCS = 0;
        auto_corr_judge(ibox) = 0;
        pos_invterval_simliarity_btw_pos_pks = 0;
        neg_invterval_simliarity_btw_pos_pks = 0;
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
bounding_boxes(~(auto_corr_judge)) = [];
final_cycles_boxes(~(auto_corr_judge)) = [];
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


while(~end_criteria)
    n_boxes = length(final_bounding_boxes);
    inner_flag = 0;
    
    for ibox = 1:n_boxes
        for jbox = 1:n_boxes
            if jbox > ibox
                %                 if is_boxes_overlapp(final_bounding_boxes(ibox),final_bounding_boxes(jbox),ovlp_threshold)
                if is_boxes_overlapped_signal(final_bounding_boxes(ibox),final_bounding_boxes(jbox),signal,srate,ovlp_threshold)
                    
                    overlap_list(cnt,:) = [ibox,jbox];
                    
                    new_box.start = min(final_bounding_boxes(ibox).start,final_bounding_boxes(jbox).start);
                    new_box.stop = max(final_bounding_boxes(ibox).stop,final_bounding_boxes(jbox).stop);
                    new_box.minF = min(final_bounding_boxes(ibox).minF,final_bounding_boxes(jbox).minF);
                    new_box.maxF = max(final_bounding_boxes(ibox).maxF,final_bounding_boxes(jbox).maxF);
                    two_fp_candidates = [final_bounding_boxes(ibox).center_fp final_bounding_boxes(jbox).center_fp];
                    two_tp_candidates = [final_bounding_boxes(ibox).center_tp final_bounding_boxes(jbox).center_tp];
                    two_pkvals = [final_bounding_boxes(ibox).peak_val final_bounding_boxes(jbox).peak_val];
                    [maxv maxi] = max(two_pkvals);
                    new_box.center_fp = two_fp_candidates(maxi);
                    new_box.center_tp = two_tp_candidates(maxi);
                    new_box.peak_val = two_pkvals(maxi);
                    
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

%%
n_clusters = length(box_clusters);
cluster_label_for_boxes = 1:length(bounding_boxes);
for icls = 1:n_clusters
    cluster_label_for_boxes(box_clusters(icls).member) = icls;
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
    imagesc(timex,freqs,denoised_env_data',[0 3*std_residue]);
    % imagesc(timex,freqs,BW');
    % imagesc(timex,freqs,ovent_map');
    set(gca, 'YDir','normal');
    ylim([.5 40]);
    %     hline([3 8 14 30]);
    hold on;
    p_bounding_boxes2 = initial_bounding_boxes;
    for ibox = 1:length(p_bounding_boxes2)
        plot(p_bounding_boxes2(ibox).center_tp,p_bounding_boxes2(ibox).center_fp,'r.','markersize',10);
        hold on;
        line([p_bounding_boxes2(ibox).start p_bounding_boxes2(ibox).stop],[p_bounding_boxes2(ibox).maxF p_bounding_boxes2(ibox).maxF],'color','red');
        line([p_bounding_boxes2(ibox).start p_bounding_boxes2(ibox).start],[p_bounding_boxes2(ibox).minF p_bounding_boxes2(ibox).maxF],'color','red');
        line([p_bounding_boxes2(ibox).stop p_bounding_boxes2(ibox).stop],[p_bounding_boxes2(ibox).minF p_bounding_boxes2(ibox).maxF],'color','red');
        line([p_bounding_boxes2(ibox).start p_bounding_boxes2(ibox).stop],[p_bounding_boxes2(ibox).minF p_bounding_boxes2(ibox).minF],'color','red');
    end
    
 
    %     set(gca,'fontsize',20);
    %     set(gca,'linewidth',linewidth_val);
    xlabel('time');
    %         xticks = 0:100:n_tp;
    set(gca,'xtick',xticks);
    set(gca,'xticklabel',round(xticks/srate,2));
    
    
    subplot(n_plots,1,3);
%     imagesc(timex,freqs,osc_data');
    imagesc(timex,freqs,denoised_env_data',[0 3*std_residue]);
    % imagesc(timex,freqs,BW');
    set(gca, 'YDir','normal');
    ylim([.5 40]);
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
outputs= oof_outputs;
outputs.bounding_boxes = final_bounding_boxes;
outputs.signal_std = std(signal);

% outputs.osc_comps = osc_comps;

end


function [outputs] = is_boxes_overlapped_signal(box1,box2,signal,srate,ovlp_threshold)
vol1 = (box1.maxF-box1.minF) * (box1.stop-box1.start);
vol2 = (box2.maxF-box2.minF) * (box2.stop-box2.start);
%%
overlapping_threshold = ovlp_threshold;
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

if abs(center_freq1 - center_freq2) > 2
    outputs = 0;
    return;
end

if length(time_intersect) == 0
    outputs = 0;
    return;
end

%%
overlapped_data1 = fft_bandpass_filtering(signal,srate,[box1.minF box1.maxF]);
overlapped_data2 = fft_bandpass_filtering(signal,srate,[box2.minF box2.maxF]);
windowing_mask = tukeywin(length(signal),0.25);
overlapped_data1 = overlapped_data1 .* windowing_mask;
overlapped_data2 = overlapped_data2 .* windowing_mask;
% angle_data1 = angle(hilbert(overlapped_data1));
% angle_data2 = angle(hilbert(overlapped_data2));
% cos_data1 = cos(angle_data1);
% cos_data2 = cos(angle_data2);

[coef, pval] = corr(overlapped_data1(time_intersect(1):time_intersect(end)),overlapped_data2(time_intersect(1):time_intersect(end)));


% if coef > overlapping_threshold || (vol_intersect>vol1*overlapping_threshold) || (vol_intersect>vol2*overlapping_threshold)
if pval < 0.05 && coef > 0
    outputs = 1;
else
    outputs = 0;
end

end

