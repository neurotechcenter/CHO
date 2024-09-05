clear;
clc;
close all;


srate = 400;

atimes = linspace(-1500,1500,srate*3);
mtimes = linspace(-1500,1500,srate*3);


frequency_range = [1:40];

%%
warning('off');
band_labels = {'\delta (<3Hz)','\theta (3-6Hz)','\alpha (7-14Hz)','\beta (15-30Hz)','low \gamma (31-40Hz)'};
band_limits = [3 6 14 30 40];
band_limits_lp = [1 3 6 14 30];
band_limits_up = [3 6 14 30 40];
n_bands = length(band_labels);

custom_colormap = jet(10);
n_subjects = 7;


%%
for isub =  1:n_subjects
    load(['y' num2str(isub) '_fooof_chpower_eeg_v1.mat']);
    
    
    %%
    max_bandpower = [];
    norm_ch_power = [];
    for iband = 1:n_bands
%         max_bandpower(iband) = max(ch_power{iband}(:));
%         norm_ch_power{iband} = ch_power{iband}/max_bandpower(iband);
        
        median_bandpower(iband) = median(ch_power{iband}(:));
        std_bandpower(iband) = std(ch_power{iband}(:));
        norm_ch_power{iband} = ch_power{iband}/(median_bandpower(iband) + 2*std_bandpower(iband));
    end
    
    
    
    h=figure('position',[10 10 1000 300]);
    for iband = 1:5
        %         subaxis(1,5,iband,'SpacingVert',0,'MR',0,'SpacingHoriz',0,'MR',0);
        %         subplot(1,5,iband);
        subaxis(1,5,iband,'SpacingHorizontal',0,'MR',0.1);
        %         ch_vect = mean(ch_duration{iband},2,'omitnan')/(srate*1.5)*10;
        ch_vect = mean(norm_ch_power{iband},2,'omitnan');
        elocs = readlocs('elec_loc_64.xyz');
        topoplot(ch_vect, elocs, 'maplimits' , [0 1] ,'plotchans' ,[1:42, 45:63]);
        % title('Auditory response at 107ms');
        ylim([-1 1]);
        %         title(band_labels{iband});
    end
%     exportgraphics(h,['osc_topo_sub_' num2str(isub) '_FOOOF_EEG_v1.jpg'],'Resolution',600);
    
    %%
    ch_bins = 1:length(elocs);
    
    ch_labels = [];
    for ch = 1:length(elocs)
        ch_labels{ch,1} = (elocs(ch).labels);
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
        
        subj_counts2{iband}(isub,:) = bar_data(:);
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
    bar_data_med = mean(subj_counts2{iband},1);
    bar_data_mad = std(subj_counts2{iband})/n_subjects^0.5;

    bh = barh(centers,bar_data_med);
    
    hold on;
    for ibin = 1:length(ch_bins)        
        line([bar_data_med(ibin)-bar_data_mad(ibin) bar_data_med(ibin)+bar_data_mad(ibin)],[ibin ibin]);
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
save('./data/subj_counts_fooof_eeg.mat','subj_counts2');
