clear;
clc;
% close all;

if exist('subj_counts_fooof_ecog.mat','file')

else
    error('Figure5_AB_CHO_v1: Run "Figure5_A_Fooof_v1.m" first');
end

srate = 400;
n_subjects = 8;

%%
warning('off');
band_labels = {'\delta (<3Hz)','\theta (3-6Hz)','\alpha (7-14Hz)','\beta (15-30Hz)','low \gamma (31-40Hz)'};
band_limits = [3 6 14 30 40];
band_limits_lp = [1 3 6 14 30];
band_limits_up = [3 6 14 30 40];
n_bands = length(band_labels);

%%

custom_colormap = jet(10);

figure();
colormap(custom_colormap);
colorbar;
set(gca, 'visible', 'off');

%%
load ctx_bins_labels_ecog.mat

%%
for isub = 1:n_subjects
    %     load(['s' num2str(isub) '_auditory_ds_v7.mat']);
    load(['x' num2str(isub) '_cho_method_1std_v1.mat']);
    
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
    
    h=figure('position',[10 10 1800 300]);
    for iband = 1:5
        subaxis(1,5,iband,'SpacingHorizontal',0,'MR',0.1);
        %         ch_vect = mean(ch_duration{iband},2,'omitnan')/(srate*1.5)*10;
        ch_vect = mean(norm_ch_power{iband},2,'omitnan');
        
        load(['x' num2str(isub) '_brain.mat']);
        if mean(tala.trielectrodes(:,1)) < 0
            hemi(isub) = 1;
        else
            hemi(isub) = -1;
        end
        %         load(['s' num2str((isub)) '_brain.mat']);
        n_ch = min([n_ch length(tala.electrodes)]);
        %     ch_vect = zeros(1,length(tala.electrodes));
        topo_param.cmapstruct = cmapstruct;
        topo_param.cortex = cortex;
        topo_param.tala = tala;
        topo_param.viewstruct = viewstruct;
        topo_param.ix = ix;
        topo_param.vcontribs = vcontribs;
        topo_param.right_hemi = hemi(isub);
        topo_param.cmapstruct.enablecolormap = true;
        topo_param.cmapstruct.enablecolorbar = false;
        topo_param.view = 'lateral';
        topo_param.alpha = 1;
        
        ecog_topo_plot(topo_param,ch_vect(1:n_ch));
        hold on;
        elecLoc = tala.trielectrodes;
        
        for ich = 1:n_ch
            textshift = 10 * hemi(isub) * -1;
            cidx = floor(ch_vect(ich) * 10)+1;
            if ch_vect(ich) == 0
                plot3(elecLoc(ich,1)+textshift,elecLoc(ich,2),elecLoc(ich,3),'o','markerfacecolor',[0 0 0]/256,'markersize',7,'markeredgecolor','none');
            elseif cidx > 10
                plot3(elecLoc(ich,1)+textshift,elecLoc(ich,2),elecLoc(ich,3),'o','markerfacecolor',custom_colormap(end,:),'markersize',7,'markeredgecolor','k');
            else
                plot3(elecLoc(ich,1)+textshift,elecLoc(ich,2),elecLoc(ich,3),'o','markerfacecolor',custom_colormap(cidx,:),'markersize',7,'markeredgecolor','k');
            end
            
        end
        %         title(band_labels{iband});
    end
    %     exportgraphics(h,['osc_topo_sub_' num2str(isub) '_CHO_ECOG_v1.jpg'],'Resolution',600);
    
    %%
    
    
    ctx_bins = 1:length(ctx_bins_labels);
    
    ch_labels = [];
    for ch = 1:n_ch
        ch_labels{ch,1} = cell2mat(SecondaryLabel{ch});
    end
    
    elec_bin_numbers = [];
    for ibin = 1:length(ctx_bins_labels)
        lb_idxs = find(strcmp(ctx_bins_labels{ibin},ch_labels));
        elec_bin_numbers(1,ibin) = length(lb_idxs);
    end
    
    %%
    brain_area = [];
    for iband = 1:length(band_labels)
        brain_area{iband} = [];
    end
    for ch = 1:n_ch
        for tr = 1:n_trials
            n_boxes = length(method_outputs(ch,tr).bounding_boxes);           

            
            for ibox = 1:n_boxes
                center_fp = method_outputs(ch,tr).bounding_boxes(ibox).center_fp;
                start_time = method_outputs(ch,tr).bounding_boxes(ibox).start;
                stop_time = method_outputs(ch,tr).bounding_boxes(ibox).stop;
                power = method_outputs(ch,tr).bounding_boxes(ibox).peak_val;
                
                for iband = 1:n_bands
                    if center_fp < band_limits_up(iband) && band_limits_lp(iband) < center_fp
%                         ch_power{iband}(ch,tr) = max(ch_power{iband}(ch,tr), power);
                        
                        ctx_number = find(strcmp(ctx_bins_labels,SecondaryLabel{ch}));
                        
                        brain_area{iband} = [brain_area{iband}; ctx_number];

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
        B = brain_area{iband};
        [counts,centers] = hist(B,ctx_bins');
        bar_data = (counts(:)./elec_bin_numbers(:))./n_trials;
        bh = barh(centers,bar_data);  

        bh.FaceColor = [.8 .8 .99];
        bh.EdgeColor = [.8 .8 .99];

        set(gca,'yTick',ctx_bins);
        ylim([ctx_bins(1)-1 ctx_bins(end)+1]);
        set(gca,'yTickLabel',[]);
        if max(counts) == 0
            xlim([0 1]);
        else
            xlim([0 max(bar_data)+max(bar_data)*0.1]);
        end
        
        title(band_labels{iband});
        
        subj_counts{iband}(isub,:) = bar_data(:);
    end
    set(h(1),'yTickLabel',ctx_bins_labels);
    xlabel('occurence per electrode');
    
   
end

%%

figure();
h = [];
for iband = 1:length(band_labels)
    h(iband) = subaxis(1,5,iband,'SpacingHorizontal',0.04,'MR',0.1);
    imagesc(subj_counts{iband}');
    set(gca,'yTick',ctx_bins);
    ylim([ctx_bins(1)-1 ctx_bins(end)+1]);
    set(gca,'yTickLabel',[]);
    set(gca, 'YDir','normal')
end
set(h(1),'yTickLabel',ctx_bins_labels);

%%
figure('position',[10 10 1000 400]);
h = [];
for iband = 1:length(band_labels)
    h(iband) = subaxis(1,5,iband,'SpacingHorizontal',0.02,'MR',0.1);
    %         h(iband) = subplot(1,5,iband);
    bar_data_med = mean(subj_counts{iband},1);
    bar_data_mad = std(subj_counts{iband})/n_subjects^0.5;

    bh = barh(centers,bar_data_med);
    
    hold on;
    for ibin = 1:length(ctx_bins)        
        line([bar_data_med(ibin)-bar_data_mad(ibin) bar_data_med(ibin)+bar_data_mad(ibin)],[ibin ibin]);
    end
    
    bh.FaceColor = [.8 .8 .99];
    bh.EdgeColor = [.8 .8 .99];
    
    set(gca,'yTick',ctx_bins);
    ylim([ctx_bins(1)-1 ctx_bins(end)+1]);
    set(gca,'yTickLabel',[]);
%     if max(counts) == 0
%         xlim([0 1]);
%     else
%         xlim([0 max(bar_data_med)+max(bar_data_med)*0.5]);
%     end
    
    title(band_labels{iband});
%     xlim([0 200]);
end
set(h(1),'yTickLabel',ctx_bins_labels);
xlabel('occurence per electrode');

%%
load subj_counts_fooof_ecog.mat

figure('position',[10 10 1000 400]);
h = [];
for iband = 1:length(band_labels)
    h(iband) = subaxis(1,5,iband,'SpacingHorizontal',0.02,'MR',0.1);
    %         h(iband) = subplot(1,5,iband);
    bar_data_med = mean(subj_counts2{iband},1);
    bar_data_mad = std(subj_counts2{iband})/n_subjects^0.5;

    bh = barh(centers,bar_data_med);
    
    hold on;
    for ibin = 1:length(ctx_bins)        
        line([bar_data_med(ibin)-bar_data_mad(ibin) bar_data_med(ibin)+bar_data_mad(ibin)],[ibin ibin]);
    end
    
    bh.FaceColor = [.8 .8 .99];
    bh.EdgeColor = [.8 .8 .99];
    
    set(gca,'yTick',ctx_bins);
    ylim([ctx_bins(1)-1 ctx_bins(end)+1]);
    set(gca,'yTickLabel',[]);
%     if max(counts) == 0
%         xlim([0 1]);
%     else
%         xlim([0 max(bar_data_med)+max(bar_data_med)*0.5]);
%     end
    
    title(band_labels{iband});
%     xlim([0 200]);
end
set(h(1),'yTickLabel',ctx_bins_labels);
xlabel('occurence per electrode');

%%
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

    bar_data_med1 = median(subj_counts{iband},1);
    bar_data_mad1 = mad(subj_counts{iband});
    
    bar_data_med2 = median(subj_counts2{iband},1);
    bar_data_mad2 = mad(subj_counts2{iband});
    
    bh = barh(centers,[bar_data_med1; bar_data_med2]);
    
    
    for ibin = 1:length(ctx_bins)
        if subj_counts2{iband}(:,ibin) == 0 | isnan(subj_counts2{iband}(:,ibin))
            p(ibin) = 1;
            rankh(ibin) = 0;
%             stats(ibin) = nan;
        else
            [p(ibin),rankh(ibin),stats(ibin)] = ranksum(subj_counts{iband}(:,ibin),subj_counts2{iband}(:,ibin));
        end        
    end  
    

    hold on;
    for ibin = 1:length(ctx_bins)        
        line([bar_data_med1(ibin) bar_data_med1(ibin)+bar_data_mad1(ibin)],[ibin-0.15 ibin-0.15]);
    end
    
    for ibin = 1:length(ctx_bins)        
        line([bar_data_med2(ibin) bar_data_med2(ibin)+bar_data_mad2(ibin)],[ibin+0.15 ibin+0.15]);
    end
    
    for ibin = 1:length(ctx_bins)        
        if p(ibin) < 0.05/(length(ctx_bins))
            text(bar_data_med2(ibin)+bar_data_mad2(ibin) ,ibin, '*', 'fontsize', 20);
        end
    end
    
    bh(1).FaceColor = [.6 .6 .99];
    bh(1).EdgeColor = [.6 .6 .99];
    
    bh(2).FaceColor = [.9 .8 .99];
    bh(2).EdgeColor = [.9 .8 .99];
    
    set(gca,'yTick',ctx_bins);
    ylim([ctx_bins(1)-1 ctx_bins(end)+1]);
    set(gca,'yTickLabel',[]);
%     if max(counts) == 0
%         xlim([0 1]);
%     else
%         xlim([0 max(bar_data_med)+max(bar_data_med)*0.5]);
%     end
    
    title(band_labels{iband});
    xlim([0 2]);
    
%     E1 = EntropyEstimationHist(bar_data_med1,centers);
%     E2 = EntropyEstimationHist(bar_data_med2,centers);
%     
%     fprintf('Entropy of CHO = %d and FOOOF = %d\n',E1, E2);
end
set(h(1),'yTickLabel',ctx_bins_labels);
xlabel('occurence per sec');
legend('CHO','FOOOF');

% E3 = EntropyEstimationHist(ones(1,length(centers))*200,centers);
%     
% fprintf('Maximal Entropy = %d\n',E3);

%% plot without empty bars

plot_idx = [2 6 11 12 14 15 16 20 21];
new_ctx_bins = 1:length(plot_idx);
new_ctx_bins_labels = ctx_bins_labels(plot_idx);
new_centers = 1:length(plot_idx);

figure('position',[10 10 1000 800]);
h = [];
for iband = 1:length(band_labels)
    h(iband) = subaxis(1,5,iband,'SpacingHorizontal',0.02,'MR',0.1);
    
    new_subj_counts{iband} = subj_counts{iband}(:,plot_idx);
    new_subj_counts2{iband} = subj_counts2{iband}(:,plot_idx);

    bar_data_med1 = median(new_subj_counts{iband},1);
    bar_data_mad1 = mad(new_subj_counts{iband});
    
    bar_data_med2 = median(new_subj_counts2{iband},1);
    bar_data_mad2 = mad(new_subj_counts2{iband});
    
    all_bar_data = [bar_data_med1; bar_data_med2];
    bh = barh(new_centers,all_bar_data);
    
    
    for ibin = 1:length(new_ctx_bins)
        if new_subj_counts2{iband}(:,ibin) == 0 | isnan(new_subj_counts2{iband}(:,ibin))
            p(ibin) = 1;
            rankh(ibin) = 0;
%             stats(ibin) = nan;
        else
            [p(ibin),rankh(ibin),stats(ibin)] = ranksum(new_subj_counts{iband}(:,ibin),new_subj_counts2{iband}(:,ibin));
        end        
    end  
    

    hold on;
    for ibin = 1:length(new_ctx_bins)        
        line([bar_data_med1(ibin) bar_data_med1(ibin)+bar_data_mad1(ibin)],[ibin-0.15 ibin-0.15]);
    end
    
    for ibin = 1:length(new_ctx_bins)        
        line([bar_data_med2(ibin) bar_data_med2(ibin)+bar_data_mad2(ibin)],[ibin+0.15 ibin+0.15]);
    end
    
    for ibin = 1:length(new_ctx_bins)        
        if p(ibin) < 0.05/(length(new_ctx_bins))
            
%             line([bar_data_med2(ibin)+bar_data_mad2(ibin) bar_data_med2(ibin)+bar_data_mad2(ibin)],[ibin-0.1 ibin+0.1],'color','k');
            text(bar_data_med2(ibin)+bar_data_mad2(ibin)+0.01 ,ibin-0.01, '*', 'fontsize', 20);
        end
    end
    
    bh(1).FaceColor = [.6 .6 .99];
    bh(1).EdgeColor = [.6 .6 .99];
    
    bh(2).FaceColor = [.9 .8 .99];
    bh(2).EdgeColor = [.9 .8 .99];
    
    set(gca,'yTick',new_ctx_bins);
    ylim([new_ctx_bins(1)-1 new_ctx_bins(end)+1]);
    set(gca,'yTickLabel',[]);
%     if max(counts) == 0
%         xlim([0 1]);
%     else
%         xlim([0 max(bar_data_med)+max(bar_data_med)*0.5]);
%     end
    
    title(band_labels{iband});
    xlim([0 1.5]);
    
%     E1 = EntropyEstimationHist(bar_data_med1,centers);
%     E2 = EntropyEstimationHist(bar_data_med2,centers);
%     
%     fprintf('Entropy of CHO = %d and FOOOF = %d\n',E1, E2);
end
set(h(1),'yTickLabel',new_ctx_bins_labels);
xlabel('occurence per sec');
legend('CHO','FOOOF');

% E3 = EntropyEstimationHist(ones(1,length(centers))*200,centers);
%     
% fprintf('Maximal Entropy = %d\n',E3);
