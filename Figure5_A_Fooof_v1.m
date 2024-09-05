clear;
clc;
close all

srate = 400;
n_subjects = 8;

%%
warning('off');
band_labels = {'\delta (<3Hz)','\theta (3-6Hz)','\alpha (7-14Hz)','\beta (15-30Hz)','low \gamma (31-40Hz)'};
band_limits = [3 6 14 30 40];
band_limits_lp = [1 3 6 14 30];
band_limits_up = [3 6 14 30 40];
n_bands = length(band_labels);

van_gogh_starry_colormap = jet(10);
load ctx_bins_labels_ecog.mat

%%
for isub = 1:8
    load(['x' num2str(isub) '_fooof_chpower_v1.mat']);
    
    
    
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
                plot3(elecLoc(ich,1)+textshift,elecLoc(ich,2),elecLoc(ich,3),'o','markerfacecolor',van_gogh_starry_colormap(end,:),'markersize',7,'markeredgecolor','k');
            else
                plot3(elecLoc(ich,1)+textshift,elecLoc(ich,2),elecLoc(ich,3),'o','markerfacecolor',van_gogh_starry_colormap(cidx,:),'markersize',7,'markeredgecolor','k');
            end
            
        end
%         title(band_labels{iband});
    end
%     exportgraphics(h,['osc_topo_sub_' num2str(isub) '_FOOOF_ECOG_v1.jpg'],'Resolution',600);

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
        
        subj_counts2{iband}(isub,:) = bar_data(:);
    end
    set(h(1),'yTickLabel',ctx_bins_labels);
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
save('./data/subj_counts_fooof_ecog.mat','subj_counts2');

