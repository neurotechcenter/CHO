clear;
clc;
close all;


srate = 400;
atimes = linspace(-1500,1500,srate*3);
mtimes = linspace(-1500,1500,srate*3);


%%
warning('off');
band_labels = {'\delta (<3Hz)','\theta (3-6Hz)','\alpha (7-14Hz)','\beta (15-30Hz)','low \gamma (31-40Hz)'};
band_limits = [3 6 14 30 40];
band_limits_lp = [1 3 6 14 30];
band_limits_up = [3 6 14 30 40];
n_bands = length(band_labels);

n_subjects = 7;

%%

custom_colormap = jet(10);

figure();
colormap(custom_colormap);
colorbar;
set(gca, 'visible', 'off');
%%
for isub =  1:n_subjects
%     load(['s' num2str(isub) '_auditory_ds_v7.mat']);
    load(['y' num2str(isub) '_cho_method_2std_EEG_v1.mat']); 
    
    n_trials = size(method_outputs,2);
    n_ch = size(method_outputs,1);
%     ch_pk_freq_prob = [];
    ch_power = [];
    
    for ch = 1:n_ch
        for tr = 1:n_trials
            n_boxes = length(method_outputs(ch,tr).bounding_boxes);
            
            for iband = 1:length(band_labels)
                ch_power{iband}(ch,tr) = 0;
            end
            
            for ibox = 1:n_boxes
                center_fp = method_outputs(ch,tr).bounding_boxes(ibox).center_fp;
                start_time = method_outputs(ch,tr).bounding_boxes(ibox).start;
                stop_time = method_outputs(ch,tr).bounding_boxes(ibox).stop;
                power = method_outputs(ch,tr).bounding_boxes(ibox).peak_val;
                noise_power = max(method_outputs(ch,tr).fooofed_ap_fit);

                for iband = 1:n_bands
                    if center_fp < band_limits_up(iband) && band_limits_lp(iband) < center_fp
%                         ch_duration_prob{iband}(ch,tr) = ch_duration_prob{iband}(ch,tr) + ch_prob;
                        ch_power{iband}(ch,tr) = max(ch_power{iband}(ch,tr), power);
                        break;
                    end
                end
                
            end
            
        end
    end
    
        %%
    max_bandpower = [];
    norm_ch_power = [];
    for iband = 1:n_bands
%         max_bandpower(iband) = max(ch_power{iband}(:));
%         norm_ch_power{iband} = ch_power{iband}/max_bandpower(iband);
        
        
        median_bandpower(iband) = median(ch_power{iband}(:));
        std_bandpower(iband) = std(ch_power{iband}(:));
        
        if (median_bandpower(iband) + 2*std_bandpower(iband)) == 0
            norm_ch_power{iband} = ch_power{iband};
        else
            norm_ch_power{iband} = ch_power{iband}/(median_bandpower(iband) + 2*std_bandpower(iband));
        end
    end
    
    
   
    
    %%

    h=figure('position',[10 10 1000 300]);
    for iband = 1:5
%         subaxis(1,5,iband,'SpacingVert',0,'MR',0); 
%         subplot(1,5,iband);
        subaxis(1,5,iband,'SpacingHorizontal',0,'MR',0.1); 
%         ch_vect = mean(ch_duration{iband},2,'omitnan')/(srate*1.5)*10;
        ch_vect = mean(norm_ch_power{iband},2,'omitnan');
        elocs = readlocs('elec_loc_64.xyz');
        topoplot(ch_vect, elocs, 'maplimits' , [0 1] ,'plotchans' ,[1:42, 45:63]);
%         topoplot(ch_vect, elocs, 'maplimits' , [0 1]);
        % title('Auditory response at 107ms');
        ylim([-1 1]);
%         title(band_labels{iband});
    end
%     exportgraphics(h,['osc_topo_sub_' num2str(isub) '_CHO_EEG_v1.jpg'],'Resolution',600);
    
    
    %%
    ch_bins = 1:length(elocs);
    
    ch_labels = [];
    for ch = 1:length(elocs)
        ch_labels{ch,1} = (elocs(ch).labels);
    end
    
    
    %%
    
    eeg_ch = [];
    for iband = 1:length(band_labels)
        eeg_ch{iband} = [];
    end
    
    for ch = 1:length(elocs)
        for tr = 1:n_trials
            n_boxes = length(method_outputs(ch,tr).bounding_boxes);
           
            
            for ibox = 1:n_boxes
                center_fp = method_outputs(ch,tr).bounding_boxes(ibox).center_fp;
                start_time = method_outputs(ch,tr).bounding_boxes(ibox).start;
                stop_time = method_outputs(ch,tr).bounding_boxes(ibox).stop;
                power = method_outputs(ch,tr).bounding_boxes(ibox).peak_val;
                noise_power = max(method_outputs(ch,tr).fooofed_ap_fit);

                for iband = 1:n_bands
                    if center_fp < band_limits_up(iband) && band_limits_lp(iband) < center_fp
%                         ch_duration_prob{iband}(ch,tr) = ch_duration_prob{iband}(ch,tr) + ch_prob;
%                         ch_power{iband}(ch,tr) = max(ch_power{iband}(ch,tr), power);
                        eeg_ch{iband} = [eeg_ch{iband}; ch];
                    end
                end
                
            end
            
        end
    end

      %%

    figure();
    h = [];
    for iband = 1:length(band_labels) 
        h(iband) = subaxis(1,5,iband,'SpacingHorizontal',0.04,'MR',0.1);
%         h(iband) = subplot(1,5,iband);        
        B = eeg_ch{iband};
        [counts,centers] = hist(B,ch_bins');
        bar_data = counts(:)./n_trials;
        bh = barh(centers,bar_data);  

        bh.FaceColor = [.8 .8 .99];
        bh.EdgeColor = [.8 .8 .99];

        set(gca,'yTick',ch_bins);
        ylim([ch_bins(1)-1 ch_bins(end)+1]);
        set(gca,'yTickLabel',[]);
        if max(counts) == 0
            xlim([0 1]);
        else
            xlim([0 max(bar_data)+max(bar_data)*0.1]);
        end
        
        title(band_labels{iband});
        
        subj_counts{iband}(isub,:) = bar_data(:);
    end
    set(h(1),'yTickLabel',ch_labels);
    xlabel('occurence per electrode');
    
end

%%
figure('position',[10 10 1000 400]);
h = [];
for iband = 1:length(band_labels)
    h(iband) = subaxis(1,5,iband,'SpacingHorizontal',0.02,'MR',0.1);
    %         h(iband) = subplot(1,5,iband);
    bar_data_med = median(subj_counts{iband},1);
    bar_data_mad = mad(subj_counts{iband})/n_subjects^0.5;

    bh = barh(centers,bar_data_med);
    
    hold on;
    for ibin = 1:length(ch_bins)        
        line([bar_data_med(ibin) bar_data_med(ibin)+bar_data_mad(ibin)],[ibin ibin]);
    end
    
    bh.FaceColor = [.8 .8 .99];
    bh.EdgeColor = [.8 .8 .99];
    
    set(gca,'yTick',ch_bins);
    ylim([ch_bins(1)-1 ch_bins(end)+1]);
    set(gca,'yTickLabel',[]);
%     if max(counts) == 0
%         xlim([0 1]);
%     else
%         xlim([0 max(bar_data_med)+max(bar_data_med)*0.5]);
%     end
    
    title(band_labels{iband});
%     xlim([0 200]);
end
set(h(1),'yTickLabel',ch_labels);
xlabel('occurence per electrode');

%%
load subj_counts_fooof_EEG.mat
figure('position',[10 10 1000 800]);
h = [];

for iband = 1:length(band_labels)
    h(iband) = subaxis(1,5,iband,'SpacingHorizontal',0.02,'MR',0.1);
    %         h(iband) = subplot(1,5,iband);
%     bar_data_med1 = mean(subj_counts{iband},1);
%     bar_data_mad1 = std(subj_counts{iband})/length(subjects)^0.5;
%     
%     bar_data_med2 = mean(subj_counts2{iband},1);
%     bar_data_mad2 = std(subj_counts2{iband})/length(subjects)^0.5;
% 
    bar_data_med1 = median(subj_counts{iband},1);
    bar_data_mad1 = mad(subj_counts{iband});
    
    bar_data_med2 = median(subj_counts2{iband},1);
    bar_data_mad2 = mad(subj_counts2{iband});
    
    bh = barh(centers,[bar_data_med1; bar_data_med2]);
    
    
    for ibin = 1:length(ch_bins)
        if subj_counts2{iband}(:,ibin) == 0 | isnan(subj_counts2{iband}(:,ibin))
            p(ibin) = 1;
            rankh(ibin) = 0;
%             stats(ibin) = nan;
        else
            [p(ibin),rankh(ibin),stats(ibin)] = ranksum(subj_counts{iband}(:,ibin),subj_counts2{iband}(:,ibin));
        end        
    end  
    

    hold on;
    for ibin = 1:length(ch_bins)        
        line([bar_data_med1(ibin) bar_data_med1(ibin)+bar_data_mad1(ibin)],[ibin-0.15 ibin-0.15]);
    end
    
    for ibin = 1:length(ch_bins)        
        line([bar_data_med2(ibin) bar_data_med2(ibin)+bar_data_mad2(ibin)],[ibin+0.15 ibin+0.15]);
    end
    
    for ibin = 1:length(ch_bins)        
        if p(ibin) < 0.05/(length(ch_bins)*1)
            text(bar_data_med2(ibin)+bar_data_mad2(ibin) ,ibin, '*', 'fontsize', 20);
        end
    end
    
    bh(1).FaceColor = [.6 .6 .99];
    bh(1).EdgeColor = [.6 .6 .99];
    
    bh(2).FaceColor = [.9 .8 .99];
    bh(2).EdgeColor = [.9 .8 .99];
    
    set(gca,'yTick',ch_bins);
    ylim([ch_bins(1)-1 ch_bins(end)+1]);
    set(gca,'yTickLabel',[]);
%     if max(counts) == 0
%         xlim([0 1]);
%     else
%         xlim([0 max(bar_data_med)+max(bar_data_med)*0.5]);
%     end
    
    title(band_labels{iband});
    xlim([0 1.5]);

    E1 = EntropyEstimationHist(bar_data_med1,centers);
    E2 = EntropyEstimationHist(bar_data_med2,centers);
    
    fprintf('Entropy of CHO = %d and FOOOF = %d\n',E1, E2);
end
set(h(1),'yTickLabel',ch_labels);
xlabel('occurence per sec');
legend('CHO','FOOOF');

E3 = EntropyEstimationHist(ones(1,length(centers))*200,centers);
    
fprintf('Maximal Entropy = %d\n',E3);





%%
FC_elec_idx = [ 5 6 ];
PO_elec_idx = [ 61 62 ];

% FC_elec_idx = [ 1:46 ];
% PO_elec_idx = [ 47:64 ];

figure('position',[10 10 1000 300]);
h = [];
for iband = 1:length(band_labels)
    h(iband) = subaxis(1,5,iband,'SpacingHorizontal',0.02,'MR',0.1);
    F_data1 = subj_counts{iband}(:,FC_elec_idx);
    F_data2 = subj_counts2{iband}(:,FC_elec_idx);
    
    PO_data1 = subj_counts{iband}(:,PO_elec_idx);
    PO_data2 = subj_counts2{iband}(:,PO_elec_idx);
    
%     bar_data = [ median([F_data1(:) PO_data1(:)],1); median([F_data2(:) PO_data2(:)],1)]';
    bar_data = [ median(F_data1(:))/length(FC_elec_idx) median(PO_data1(:))/length(PO_elec_idx); median(F_data2(:))/length(FC_elec_idx) median(PO_data2(:))/length(PO_elec_idx)];
    bh = barh(bar_data);
%     
%     bh(1).FaceColor = [.6 .6 .99];
%     bh(1).EdgeColor = [.6 .6 .99];
%     
%     bh(2).FaceColor = [.9 .8 .99];
%     bh(2).EdgeColor = [.9 .8 .99];
    
    bh(1).FaceColor = [.99 .6 .6 ];
    bh(1).EdgeColor = [.99 .6 .6 ];
    
    bh(2).FaceColor = [.99 .9 .8];
    bh(2).EdgeColor = [.99 .9 .8];
    
    [cho_p(iband),F_rankh(iband),F_stats(iband)] = ranksum(F_data1(:),PO_data1(:));
    [fooof_p(iband),PO_rankh(iband),PO_stats(iband)] = ranksum(F_data2(:),PO_data2(:));
    
% 
    line([median(PO_data1(:))/length(PO_elec_idx) median(PO_data1(:))/length(PO_elec_idx)+mad(PO_data1(:))/length(PO_elec_idx)],[1+0.15 1+0.15]);
    line([median(F_data1(:))/length(FC_elec_idx) median(F_data1(:))/length(FC_elec_idx)+mad(F_data1(:))/length(FC_elec_idx)],[1-0.15 1-0.15]);
    
    line([median(PO_data2(:))/length(PO_elec_idx) median(PO_data2(:))/length(PO_elec_idx)+mad(PO_data2(:))/length(PO_elec_idx)],[2+0.15 2+0.15]);
    line([median(F_data2(:))/length(FC_elec_idx) median(F_data2(:))/length(FC_elec_idx)+mad(F_data2(:))/length(FC_elec_idx)],[2-0.15 2-0.15]);


%    
    if cho_p(iband) < 0.05/(2)
        text(max([median(F_data1(:))/length(FC_elec_idx) median(PO_data1(:))/length(PO_elec_idx)]) ,1, '*', 'fontsize', 20);
    end
%     
    if fooof_p(iband) < 0.05/(2)
        text(max([median(F_data2(:))/length(FC_elec_idx) median(PO_data2(:))/length(PO_elec_idx)]) ,2, '*', 'fontsize', 20);
    end
% %         

    title(band_labels{iband});
    set(gca,'yTickLabel',[]);
    
    xlim([0 0.7]);
end
set(h(1),'yTickLabel',{'CHO','FOOOF'});
xlabel('occurence per sec');
legend('FC2-FC4','O1-O2');




