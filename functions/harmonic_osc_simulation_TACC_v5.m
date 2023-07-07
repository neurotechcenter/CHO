function [ACC sensitivity specificity avg_SNR] = harmonic_osc_simulation_TACC_v5(sim_param)

warning('off');
ACC = [];
sensitivity = [];
specificity = [];

%% simulation details for this video

srate  = sim_param.srate; % sampling rate in Hz
time   = -1:1/srate:4;
pnts   = length(time);
osc_amp = sim_param.osc_amp*2;

freq_range = sim_param.freq_range;

noise_level = sim_param.noise_amp;

fprintf('SNR = %d\n',osc_amp*0.5/noise_level);

n_iteration = sim_param.n_iteration;

final_center_freq1 = [];
final_center_freq2 = [];
all_outputs = [];
frequency_range = [1:40];


for sim_iter = 1:n_iteration
    lineLength = fprintf('%d/%d',sim_iter,n_iteration);
    osc_freq = randi(freq_range,1,1);
    sim_freqs(sim_iter) = osc_freq;
    w_params = [ 2.5/osc_freq 1 3];
%     w_params = [1/osc_freq 2.5/osc_freq 1 2];
    
    
    Lw = round(srate*1/osc_freq); % add less than 2 cycle oscillation
    H = tukeywin(Lw,0.40);
    H = H(1:Lw);
    H = H/max(H);
    
    masking_signal = zeros(size(time));
    rand_tp = randi(pnts-Lw);
    masking_signal(rand_tp:rand_tp+Lw-1) = H;
    
    beta_scaler = sim_param.asym_ratio; %    0< beta_scaler < 1
    skew = 0; % -pi/4 < skew < pi/4
    rand_phase = rand*2*pi;
    y_x=noise_level*sin(2*pi*osc_freq*time);
    y=noise_level*sin(2*pi*osc_freq*time-skew*y_x + rand_phase);
    y(y > 0) =y(y > 0)*beta_scaler;
    y(y< 0) =y(y < 0)*(1-beta_scaler);
    short_spike = y.*masking_signal;
    
    noise = noise_level*pinknoise(size(time,1),size(time,2)) + short_spike;
    
    for w_iter = 1:length(w_params)
        % frequencies to simulate
        Lw = round(srate*w_params(w_iter));
        %     H = hann(Lw);
        H = tukeywin(Lw,0.40);
        H = H(1:Lw);
        H = H/max(H);
        
        masking_signal = zeros(size(time));
        rand_tp = randi(pnts-Lw);
        masking_signal(rand_tp:rand_tp+Lw-1) = H;
        
        
        %%
        beta_scaler = sim_param.asym_ratio; %    0< beta_scaler < 1
        skew = 0; % -pi/4 < skew < pi/4
        rand_phase = rand*2*pi;
        y_x=osc_amp*sin(2*pi*osc_freq*time);
        y=osc_amp*sin(2*pi*osc_freq*time-skew*y_x + rand_phase);
        y(y > 0) =y(y > 0)*beta_scaler;
        y(y< 0) =y(y < 0)*(1-beta_scaler);
        oscillations{w_iter} = y.*masking_signal;
        SNRs(sim_iter,w_iter) = snr(oscillations{w_iter},noise);
        sim_timewindow{w_iter,sim_iter} = masking_signal;
        
        %     signal{w_iter} = noise+ oscillations{w_iter} + oscillations1{w_iter};
        signal{w_iter} = noise + oscillations{w_iter};

        
    end
    
    %% CHO
    max_idxs = [];
    for w_iter = 1:length(w_params)
        % frequencies to simulation
        param.plot = 0;
        param.minimum_cycles = 2;
        param.ovlp_threshold = 0.50;
        param.frequency_vector = 1:40;
        
        [cho_outputs(w_iter)] = CHO_v21(signal{w_iter}', srate, param);
        
        peak_vals = [];
        center_freqs = [];
        box_signals = [];
        for iboxes = 1:length(cho_outputs(w_iter).bounding_boxes)
            peak_vals(iboxes) = cho_outputs(w_iter).bounding_boxes(iboxes).peak_val;
            center_freqs(iboxes) = cho_outputs(w_iter).bounding_boxes(iboxes).center_fp;
            box_starts(iboxes) = cho_outputs(w_iter).bounding_boxes(iboxes).start;
            box_ends(iboxes) = cho_outputs(w_iter).bounding_boxes(iboxes).stop;
            box_signals(iboxes,:) = zeros(size(time));
            box_signals(iboxes,box_starts(iboxes):box_ends(iboxes)) = 1;
        end
        if isempty(peak_vals)
            max_idxs(w_iter) = -1;
            final_center_freq{1}{w_iter,sim_iter} = [];
            final_center_wind{1}{w_iter,sim_iter} = [];
        else
            final_center_freq{1}{w_iter,sim_iter} = center_freqs;
            final_center_wind{1}{w_iter,sim_iter} = box_signals;
        end
    end
    
    %% OEvent
    max_idxs = [];
    for w_iter = 1:length(w_params)
        % frequencies to simulation
        param.plot = 0;
        [OEvent_outputs(w_iter)] = OEvent_Neymotin_v1(signal{w_iter}', srate, param.plot);
        peak_vals = [];
        center_freqs = [];
        box_signals = [];
        for iboxes = 1:length(OEvent_outputs(w_iter).bounding_boxes)
            peak_vals(iboxes) = OEvent_outputs(w_iter).bounding_boxes(iboxes).peak_val;
            center_freqs(iboxes) = OEvent_outputs(w_iter).bounding_boxes(iboxes).center_fp;
            box_starts(iboxes) = OEvent_outputs(w_iter).bounding_boxes(iboxes).start;
            box_ends(iboxes) = OEvent_outputs(w_iter).bounding_boxes(iboxes).stop;
            box_signals(iboxes,:) = zeros(size(time));
            box_signals(iboxes,box_starts(iboxes):box_ends(iboxes)) = 1;
        end
        if isempty(peak_vals)
            max_idxs(w_iter) = -1;
            final_center_freq{2}{w_iter,sim_iter} = [];
            final_center_wind{2}{w_iter,sim_iter} = [];
        else
            final_center_freq{2}{w_iter,sim_iter} = center_freqs;
            final_center_wind{2}{w_iter,sim_iter} = box_signals;
        end
    end
    

    
    %%
    fprintf(repmat('\b',1,lineLength));
end
fprintf('\n');

%% spectral sensitivity and specificity
clc;
n_methods = 2;
for imethod = 1:n_methods
    TP{imethod} = zeros(length(w_params),1);
    TN{imethod} = zeros(length(w_params),1);
    FP{imethod} = zeros(length(w_params),1);
    FN{imethod} = zeros(length(w_params),1);
end

concede_threshold = 1.5;

for w_iter = 1:length(w_params)
    for sim_iter = 1:n_iteration
        
        tmp1 = [];
        for imethod = 1:n_methods
            tmp1{w_iter,sim_iter} = abs(final_center_freq{imethod}{w_iter,sim_iter} - sim_freqs(sim_iter));
            tp_idxs = find(tmp1{w_iter,sim_iter} <= concede_threshold);
            n_boxes = size(final_center_wind{imethod}{w_iter,sim_iter},1);            
            fp_idxs = find(tmp1{w_iter,sim_iter} > concede_threshold);
            n_correct_boxes = 0;
            if ~isempty(tp_idxs)
                for iboxes = 1:n_boxes
                    x = sim_timewindow{w_iter,sim_iter};
                    y = final_center_wind{imethod}{w_iter,sim_iter}(iboxes,:);
                    [cval pval] = corr(x',y');                    
                    if cval > 0.5 && pval < 0.05 && tmp1{w_iter,sim_iter}(iboxes) <= 1
                        n_correct_boxes = n_correct_boxes + 1;  
                    else
                        FP{imethod}(w_iter) = FP{imethod}(w_iter) + 1;
                    end
                end
                
                if n_correct_boxes > 0
                    TP{imethod}(w_iter) = TP{imethod}(w_iter) + 1;
                end

            else
                FN{imethod}(w_iter) = FN{imethod}(w_iter) + 1;
            end
            
            FP{imethod}(w_iter) = FP{imethod}(w_iter) + length(fp_idxs);
            
            if isempty(fp_idxs) && length(tp_idxs) >= 1
                TN{imethod}(w_iter) = TN{imethod}(w_iter) + 1;
            end           
            
            
%             if isempty(tp_idxs)
%                 FN{imethod}(w_iter) = FN{imethod}(w_iter) + 1;
%             else
%                 for iboxes = 1:n_boxes
%                     x = sim_timewindow{w_iter,sim_iter};
%                     y = final_center_wind{imethod}{w_iter,sim_iter}(iboxes,:);
%                     [cval pval] = corr(x',y');
%                     
%                     if cval > 0.5 && pval < 0.05 && tmp1{w_iter,sim_iter}(iboxes) <= 1
%                         n_correct_boxes = n_correct_boxes + 1;
%                     else
%                         FP{imethod}(w_iter) = FP{imethod}(w_iter) + 1;
%                     end
%                 end
%                 
%                 if n_correct_boxes > 0
%                     TP{imethod}(w_iter) = TP{imethod}(w_iter) + 1;
%                 end
%                 if n_correct_boxes == 1 && n_boxes == 1
%                     TN{imethod}(w_iter) = TN{imethod}(w_iter) + 1;
%                 end
%             end
        end   

    end
    for imethod = 1:n_methods
        sensitivity(imethod,w_iter) = TP{imethod}(w_iter)/(TP{imethod}(w_iter)+FN{imethod}(w_iter));
        specificity(imethod,w_iter) = TN{imethod}(w_iter)/(TN{imethod}(w_iter)+FP{imethod}(w_iter));
        ACC(imethod,w_iter) = (TP{imethod}(w_iter)+TN{imethod}(w_iter))/(TP{imethod}(w_iter)+FN{imethod}(w_iter)+TN{imethod}(w_iter)+FP{imethod}(w_iter));
    end
end


%%
disp('Spetro-temporal Acc');
ACC
sensitivity
specificity
avg_SNR = mean(SNRs)

%%
% figure();
% scatter(1-sensitivity(1),specificity(1));
% xlim([0 1]);
% ylim([0 1]);