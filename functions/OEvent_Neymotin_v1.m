function [outputs] = OEvent_Neymotin_v1(signal, srate, plot_on_off)

outputs = [];

minimum_cycles = 2;
%%
% [psd,f] = pwelch(signal,[],srate/2,1:srate/2,srate);
% 
% figure();
% plot(f,log10(psd));
% xlim([0 40]);

%%
% cycles = 1;
% freq = 6;
% [ filtered_data ] = wavelet_single_freq_filtering_v2(signal, srate, freq, cycles);

%%
timex = 1:length(signal);
f_steps = 1;
freqs = f_steps:f_steps:40;
pos_norm_env_data = zeros(length(timex),length(freqs));
for ifreq = 1:length(freqs)
    freq = freqs(ifreq);
    [ filtered_data(:,ifreq) ] = fft_bandpass_filtering(signal, srate, [freq-f_steps freq+f_steps]);
%     [ filtered_data(:,ifreq) ] = wavelet_single_freq_filtering_v2(signal, srate, ifreq, cycles);
    env_data(:,ifreq) = abs(hilbert(filtered_data(:,ifreq)));
    med_env(ifreq) = median(env_data(:,ifreq));
end

%% remove pink noise

for itp = 1:length(timex)
    denoised_env_data(itp,:) = env_data(itp,:)' - med_env';
%     denoised_env_data(itp,:) = env_data(itp,:)' - ap_fit;
end

std_residue = std(denoised_env_data(:));

%% masking osc data
masked_env_data = [];
osc_data = zeros(size(denoised_env_data));
for ifp = 1:length(freqs)
    sig_idxs = find(denoised_env_data(:,ifp)>1*std_residue);
    oevent_map(sig_idxs,ifp) = 1;
    osc_data(sig_idxs,ifp) = filtered_data(sig_idxs,ifp) .* oevent_map(sig_idxs,ifp);
    masked_env_data(sig_idxs,ifp) = denoised_env_data(sig_idxs,ifp) .* oevent_map(sig_idxs,ifp);
end

%% find local power peaks
BW = imregionalmax(masked_env_data);
[pks_row pks_col] = find(BW > 0);

%% find boxes
n_peaks = length(pks_row);
n_freqs = length(freqs);
n_times = length(timex);
box_threshold = 1*std_residue;
for ipks = 1:n_peaks
    tp = pks_row(ipks);
    fp = pks_col(ipks);
    peak_val = denoised_env_data(tp,fp);
    
%     box_threshold = peak_val/2;
    
    bounding_boxes(ipks).center_fp = fp;
    bounding_boxes(ipks).center_tp = tp;
    bounding_boxes(ipks).peak_val = peak_val;
    % find minF
    for ifp = fp:-1:1
%         if denoised_env_data(tp,ifp) < peak_val/2 || denoised_env_data(tp,ifp) < 0
        if denoised_env_data(tp,ifp) < peak_val/2 || denoised_env_data(tp,ifp) < 0
            break;
        end
    end
    bounding_boxes(ipks).minF = ifp;
    
    % find maxF
    for ifp = fp:n_freqs
        if denoised_env_data(tp,ifp) < peak_val/2 || denoised_env_data(tp,ifp) < 0
            break;
        end
    end
    bounding_boxes(ipks).maxF = ifp;
    
    % find start
    for itp = tp:-1:1
        if denoised_env_data(itp,fp) < peak_val/2 || denoised_env_data(itp,fp) < 0
            break;
        end
    end
    bounding_boxes(ipks).start = itp;
    
    % find maxF
    for itp = tp:n_times
        if denoised_env_data(itp,fp) < peak_val/2 || denoised_env_data(itp,fp) < 0
            break;
        end
    end
    bounding_boxes(ipks).stop = itp;
    
    cycles_boxes(ipks) = fp * (bounding_boxes(ipks).stop - bounding_boxes(ipks).start)/srate;
end
initial_bounding_boxes = bounding_boxes;

%% remove bounding boxes that have less than 2 cycle

final_cycles_boxes = cycles_boxes;
bounding_boxes(cycles_boxes < minimum_cycles) = [];
final_cycles_boxes(cycles_boxes < minimum_cycles) = [];


%% merge overlapping boxes (50%)
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
                if is_boxes_overlapp(final_bounding_boxes(ibox),final_bounding_boxes(jbox))
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

%%
n_boxes = length(final_bounding_boxes);
osc_comps = [];
for ibox = 1:n_boxes
    osc_comps(:,ibox) = fft_bandpass_filtering(signal,srate,[final_bounding_boxes(ibox).minF final_bounding_boxes(ibox).maxF]);
    L = final_bounding_boxes(ibox).stop - final_bounding_boxes(ibox).start;
    windowing_mask = tukeywin(L,0.25);
    zero_data = zeros(length(signal),1);
    zero_data(final_bounding_boxes(ibox).start:final_bounding_boxes(ibox).stop-1) = windowing_mask;
    osc_comps(:,ibox) = osc_comps(:,ibox) .* zero_data;
end

% periodic_signal = sum(osc_comps,2);
% aperiodic_signal = signal - periodic_signal;

%% visualization

if plot_on_off > 0
    n_plots = 4;
    figure('position',[10 10 400 1000]);
    subplot(n_plots,1,1);
    plot(timex,signal);
    hold on;
    for ibox = 1:length(final_bounding_boxes)
        %     figure();
        min_val = min(signal(final_bounding_boxes(ibox).start:final_bounding_boxes(ibox).stop));
        max_val = max(signal(final_bounding_boxes(ibox).start:final_bounding_boxes(ibox).stop));
        line([final_bounding_boxes(ibox).start final_bounding_boxes(ibox).stop],[max_val max_val],'color','red');
        line([final_bounding_boxes(ibox).start final_bounding_boxes(ibox).start],[min_val max_val],'color','red');
        line([final_bounding_boxes(ibox).stop final_bounding_boxes(ibox).stop],[min_val max_val],'color','red');
        line([final_bounding_boxes(ibox).start final_bounding_boxes(ibox).stop],[min_val min_val],'color','red');
        text((final_bounding_boxes(ibox).start+final_bounding_boxes(ibox).stop)/2,0, [num2str(final_bounding_boxes(ibox).center_fp) 'Hz'],'color','red');
        fprintf('\n%dHz\n', final_bounding_boxes(ibox).center_fp);
    end
    
    
    subplot(n_plots,1,2);
    imagesc(timex,freqs,denoised_env_data',[std_residue 6*std_residue]);
    % imagesc(timex,freqs,BW');
    % imagesc(timex,freqs,ovent_map');
    set(gca, 'YDir','normal');
    ylim([.5 40]);
%     hline([3 8 14 30]);
    hold on;
    
    for ibox = 1:length(initial_bounding_boxes)
        plot(pks_row(ibox),pks_col(ibox),'r.','markersize',10);
        hold on;
        line([initial_bounding_boxes(ibox).start initial_bounding_boxes(ibox).stop],[initial_bounding_boxes(ibox).maxF initial_bounding_boxes(ibox).maxF],'color','red');
        line([initial_bounding_boxes(ibox).start initial_bounding_boxes(ibox).start],[initial_bounding_boxes(ibox).minF initial_bounding_boxes(ibox).maxF],'color','red');
        line([initial_bounding_boxes(ibox).stop initial_bounding_boxes(ibox).stop],[initial_bounding_boxes(ibox).minF initial_bounding_boxes(ibox).maxF],'color','red');
        line([initial_bounding_boxes(ibox).start initial_bounding_boxes(ibox).stop],[initial_bounding_boxes(ibox).minF initial_bounding_boxes(ibox).minF],'color','red');
    end
    
    
    subplot(n_plots,1,3);
    imagesc(timex,freqs,osc_data');
    % imagesc(timex,freqs,BW');
    set(gca, 'YDir','normal');
    ylim([.5 40]);
%     hline([3 8 14 30]);
    hold on;
    
    for ibox = 1:length(final_bounding_boxes)
        %     figure();
        plot(final_bounding_boxes(ibox).center_tp,final_bounding_boxes(ibox).center_fp,'r.','markersize',10);
        hold on;
        line([final_bounding_boxes(ibox).start final_bounding_boxes(ibox).stop],[final_bounding_boxes(ibox).maxF final_bounding_boxes(ibox).maxF],'color','red');
        line([final_bounding_boxes(ibox).start final_bounding_boxes(ibox).start],[final_bounding_boxes(ibox).minF final_bounding_boxes(ibox).maxF],'color','red');
        line([final_bounding_boxes(ibox).stop final_bounding_boxes(ibox).stop],[final_bounding_boxes(ibox).minF final_bounding_boxes(ibox).maxF],'color','red');
        line([final_bounding_boxes(ibox).start final_bounding_boxes(ibox).stop],[final_bounding_boxes(ibox).minF final_bounding_boxes(ibox).minF],'color','red');
    end
    
    subplot(n_plots,1,4);
    plot(freqs,log10(med_env));
%     plot(freqs,ap_fit);
    
%     eegplot([signal osc_comps periodic_signal aperiodic_signal]','srate',srate,'winlength',1);
end


%%
% eegplot([signal osc_comps]','srate',srate,'winlength',1);

%%
outputs.bounding_boxes = final_bounding_boxes;
outputs.osc_comps = osc_comps;

end

function [outputs] = is_boxes_overlapp(box1,box2)
vol1 = (box1.maxF-box1.minF) * (box1.stop-box1.start);
vol2 = (box2.maxF-box2.minF) * (box2.stop-box2.start);
%%
overlapping_threshold = 0.5;
%%
time1 = box1.start:box1.stop;
time2 = box2.start:box2.stop;
freq1 = box1.minF:box1.maxF;
freq2 = box2.minF:box2.maxF;
[time_intersect] = intersect(time1,time2);
[freq_intersect] = intersect(freq1,freq2);

vol_intersect = length(time_intersect)*length(freq_intersect);
if (vol_intersect>vol1*overlapping_threshold) || (vol_intersect>vol2*overlapping_threshold)
    outputs = 1;
else
    outputs = 0;
end

end