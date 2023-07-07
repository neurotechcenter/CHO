function [ACC sensitivity specificity avg_SNR] = harmonic_osc_simulation_for_different_method_v7(sim_param)

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

% w_params = [2.5/osc_freq 1 3];
% w_params = sim_param.w_params;

all_outputs = [];
frequency_range = [1:40];


for sim_iter = 1:n_iteration
    lineLength = fprintf('%d/%d',sim_iter,n_iteration);
    osc_freq = randi(freq_range,1,1);
    sim_freqs(sim_iter) = osc_freq;
    w_params = [2.5/osc_freq 1 3];
%     w_params = [0 1 3];
           
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
    
    pink_noise = pinknoise(size(time,1),size(time,2));
%     pink_noise = pink_noise/max(abs(pink_noise));
    noise = noise_level*pink_noise + short_spike;
    
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
        
        r = snr(oscillations{w_iter},noise);
%         r1 = mag2db(rssq(oscillations{w_iter})/rssq(noise));
        
        %     signal{w_iter} = noise+ oscillations{w_iter} + oscillations1{w_iter};
        signal{w_iter} = noise + oscillations{w_iter};
        SNR(sim_iter,w_iter) = r;

        %% run fooof
        fparam = [];
        fparam.srate = srate;
        fparam.slide_window_length_sec = 1;
        fparam.window_overlap = 0.5;
        fparam.freq_range = [1 40];

        [fractal,original] = ft_fooof_matlab_v1(signal{w_iter},fparam);
        
        all_outputs(sim_iter,w_iter).fractal = fractal;
        all_outputs(sim_iter,w_iter).original = original;
        
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
        

        %%
        peak_vals = [];
        center_freqs = [];
        for iboxes = 1:length(cho_outputs(w_iter).bounding_boxes)
            peak_vals(iboxes) = cho_outputs(w_iter).bounding_boxes(iboxes).peak_val;
            center_freqs(iboxes) = cho_outputs(w_iter).bounding_boxes(iboxes).center_fp;
        end
        if isempty(peak_vals)
            max_idxs(w_iter) = -1;
            final_center_freq{1}{w_iter,sim_iter} = [];
        else
            final_center_freq{1}{w_iter,sim_iter} = center_freqs;
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
        for iboxes = 1:length(OEvent_outputs(w_iter).bounding_boxes)
            peak_vals(iboxes) = OEvent_outputs(w_iter).bounding_boxes(iboxes).peak_val;
            center_freqs(iboxes) = OEvent_outputs(w_iter).bounding_boxes(iboxes).center_fp;
        end
        if isempty(peak_vals)
            max_idxs(w_iter) = -1;
            final_center_freq{2}{w_iter,sim_iter} = [];
        else
            final_center_freq{2}{w_iter,sim_iter} = center_freqs;
        end
    end   
    
    
    %% fooof
    max_idxs = [];
    for w_iter = 1:length(w_params)
        % frequencies to simulation
        peak_vals = [];
        center_freqs = [];
        for igaussian = 1:size(all_outputs(sim_iter,w_iter).fractal.fooofparams.peak_params,1)
            peak_vals(igaussian) = all_outputs(sim_iter,w_iter).fractal.fooofparams.peak_params(igaussian,2);
            center_freqs(igaussian) = all_outputs(sim_iter,w_iter).fractal.fooofparams.peak_params(igaussian,1);
        end
        if isempty(peak_vals)
            max_idxs(w_iter) = -1;
            final_center_freq{3}{w_iter,sim_iter} = [];
        else
            final_center_freq{3}{w_iter,sim_iter} = center_freqs;
        end
        %%
%         osc_freq
%         center_freqs
    end
    
    %% SPRiNT
    max_idxs = [];
    for w_iter = 1:length(w_params)
        % Inputs 
        F = signal{w_iter};                             % Input time series
        % STFT opts
        opt.sfreq = srate;                    % Input sampling rate
        opt.WinLength = 1;                  % STFT window length
        opt.WinOverlap = 50;                % Overlap between sliding windows (in %)
        opt.WinAverage = 5;                 % Number of sliding windows averaged by time point
        % specparam opts
        opt.freq_range          = [1 40];
        opt.peak_width_limits   = [1.5 6];
        opt.max_peaks           = 3;
        opt.min_peak_height     = 6 / 10; % convert from dB to B
        opt.aperiodic_mode      = 'fixed'; % alternative: knee
        opt.peak_threshold      = 2.0;   % 2 std dev: parameter for interface simplification
        % Matlab-only options
        opt.peak_type           = 'gaussian'; % alternative: cauchy
        opt.proximity_threshold = 2;
        opt.guess_weight        = 'none';
        opt.thresh_after        = true;   % Threshold after fitting, always selected for Matlab 
                                          % (mirrors the Python FOOOF closest by removing peaks
                                          % that do not satisfy a user's predetermined conditions)
                                          % only used in the absence of the
        if license('test','optimization_toolbox') % check for optimization toolbox
            opt.hOT = 1;
            disp('Using constrained optimization, Guess Weight ignored.')
        else
            opt.hOT = 0;
            disp('Using unconstrained optimization, with Guess Weights.')
        end
        opt.rmoutliers          = 'yes';
        opt.maxfreq             = 2.5;
        opt.maxtime             = 6;
        opt.minnear             = 3;  

        Freqs = 0:1/opt.WinLength:opt.sfreq/2;
        channel = struct('data',[],'peaks',[],'aperiodics',[],'stats',[]);
        % Compute short-time Fourier transform
        [TF, ts] = SPRiNT_stft(F,opt);
        outputStruct = struct('opts',opt,'freqs',Freqs,'channel',channel);

        % Parameterize STFTs
        s_data = SPRiNT_specparam_matlab(TF,outputStruct.freqs,outputStruct.opts,ts);
        STF = squeeze(TF);
        SPRiNT_output(w_iter) = s_data;
        
        %
        %%
%         n_plots = 4;
%         figure('position',[10 10 400 1000]);
        
%         subplot(n_plots,1,1);
%         plot(F);
%         hold on;
%         plot(oscillations{w_iter});
%         xlabel('time');
%         title(num2str(osc_freq));
%         
%         subplot(n_plots,1,2);
%         imagesc(ts,Freqs,10*log10(STF'));
%         set(gca, 'YDir','normal');
%         ylim([0 40]);
%         xlabel('time');
%         
%         subplot(n_plots,1,3);
%         plot_data = squeeze(s_data.SPRiNT.peak_models)';
%         imagesc(ts,s_data.Freqs,plot_data);
%         % ylim([0 40]);
%         set(gca, 'YDir','normal');
%         xlabel('time');
%         
%         subplot(n_plots,1,4);
%         % plot_data = squeeze(s_data.SPRiNT.SPRiNT_models)';
%         plot_data = squeeze(s_data.SPRiNT.peak_models)';
%         plot(10*log10(plot_data));
%         xlim([0 40]);
%         xlabel('Hz');

        
        %%
        peak_vals = [];
        center_freqs = [];
        for ipeaks = 1:length(s_data.SPRiNT.channel.peaks)
            peak_vals(ipeaks) = SPRiNT_output(w_iter).SPRiNT.channel.peaks(ipeaks).amplitude;
            center_freqs(ipeaks) = SPRiNT_output(w_iter).SPRiNT.channel.peaks(ipeaks).center_frequency;
        end
        if isempty(peak_vals)
            max_idxs(w_iter) = -1;
            final_center_freq{4}{w_iter,sim_iter} = [];
        else
            final_center_freq{4}{w_iter,sim_iter} = center_freqs;
        end
    end
    
    %%
    fprintf(repmat('\b',1,lineLength));
end
fprintf('\n');

%% spectral sensitivity and specificity
clc;
n_methods = 4;
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
            fp_idxs = find(tmp1{w_iter,sim_iter} > concede_threshold);
            if length(tp_idxs) > 0
                TP{imethod}(w_iter) = TP{imethod}(w_iter) + 1;
            else
                FN{imethod}(w_iter) = FN{imethod}(w_iter) + 1;
            end
            FP{imethod}(w_iter) = FP{imethod}(w_iter) + length(fp_idxs);
            if isempty(fp_idxs) && length(tp_idxs) >= 1
                TN{imethod}(w_iter) = TN{imethod}(w_iter) + 1;
            end
        end             
    end
    
    for imethod = 1:n_methods
        sensitivity(imethod,w_iter) = TP{imethod}(w_iter)/(TP{imethod}(w_iter)+FN{imethod}(w_iter));
        specificity(imethod,w_iter) = TN{imethod}(w_iter)/(TN{imethod}(w_iter)+FP{imethod}(w_iter));
        ACC(imethod,w_iter) = (TP{imethod}(w_iter)+TN{imethod}(w_iter))/(TP{imethod}(w_iter)+FN{imethod}(w_iter)+TN{imethod}(w_iter)+FP{imethod}(w_iter));
    end
end

%%
disp('Spectral Acc');
ACC
sensitivity
specificity
avg_SNR = mean(SNR)

end
%%  Step 1: Produce STFTs
function [TF, ts] = SPRiNT_stft(F,opts)
% SPRiNT_stft: Compute a locally averaged short-time Fourier transform (for
% use in SPRiNT)
% 
% Segments of this function were adapted from the Brainstorm software package:
% https://neuroimage.usc.edu/brainstorm
% Tadel et al. (2011) 
% 
% Copyright (c)2000-2020 University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
% 
% Author: Luc Wilson (2022)

    sfreq = opts.sfreq;                 % sample rate, in Hertz
    WinLength = opts.WinLength;         % window length, in seconds
    WinOverlap = opts.WinOverlap;       % window overlap, in percent
    avgWin = opts.WinAverage;           % number of windows being averaged per PSD
    nTime = size(F,2);
    % ===== WINDOWING =====
    Lwin  = round(WinLength * sfreq); % number of data points in windows
    Loverlap = round(Lwin * WinOverlap / 100); % number of data points in overlap
    % If window is too small
    if (Lwin < 50)
        return;
    % If window is bigger than the data
    elseif (Lwin > nTime)
        Lwin = size(F,2);
        Lwin = Lwin - mod(Lwin,2); % Make sure the number of samples is even
        Loverlap = 0;
        Nwin = 1;
    % Else: there is at least one full time window
    else
        Lwin = Lwin - mod(Lwin,2);    % Make sure the number of samples is even
        Nwin = floor((nTime - Loverlap) ./ (Lwin - Loverlap));
    end
    % Next power of 2 from length of signal
    NFFT = Lwin;                    % No zero-padding: Nfft = Ntime 
    % Positive frequency bins spanned by FFT
    FreqVector = sfreq / 2 * linspace(0,1,NFFT/2+1);
    % Determine hann window shape/power
    Win = hann(Lwin)';
    WinNoisePowerGain = sum(Win.^2);
    % Initialize STFT,time matrices
    ts = nan(Nwin-(avgWin-1),1);
    TF = nan(size(F,1), Nwin-(avgWin-1), size(FreqVector,2));
    TFtmp = nan(size(F,1), avgWin, size(FreqVector,2));
    % ===== CALCULATE FFT FOR EACH WINDOW =====
    TFfull = zeros(size(F,1),Nwin,size(FreqVector,2));
    for iWin = 1:Nwin
        % Build indices
        iTimes = (1:Lwin) + (iWin-1)*(Lwin - Loverlap);
        center_time = floor(median((iTimes-(avgWin-1)./2*(Lwin - Loverlap))))./sfreq;
        % Select indices
        Fwin = F(:,iTimes);
        % No need to enforce removing DC component (0 frequency).
        Fwin = Fwin - mean(Fwin,2);
        % Apply a Hann window to signal
        Fwin = Fwin .* Win;
        % Compute FFT
        Ffft = fft(Fwin, NFFT, 2);
        % One-sided spectrum (keep only first half)
        % (x2 to recover full power from negative frequencies)
        TFwin = Ffft(:,1:NFFT/2+1) * sqrt(2 ./ (sfreq * WinNoisePowerGain));
        % x2 doesn't apply to DC and Nyquist.
        TFwin(:, [1,end]) = TFwin(:, [1,end]) ./ sqrt(2);
        % Permute dimensions: time and frequency
        TFwin = permute(TFwin, [1 3 2]);
        % Convert to power
        TFwin = abs(TFwin).^2;
        TFfull(:,iWin,:) = TFwin;
        TFtmp(:,mod(iWin,avgWin)+1,:) = TFwin;
        if isnan(TFtmp(1,1,1)) 
            continue % Do not record anything until transient is gone
        else
    %     Save STFTs for window
        TF(:,iWin-(avgWin-1),:) = mean(TFtmp,2);
        ts(iWin-(avgWin-1)) = center_time;
        end
    end
end

%%  Step 2: parameterize spectrograms spectra
function [s_data] = SPRiNT_specparam_matlab(TF, fs, opt, ts)
% SPRiNT_specparam_matlab: Compute time-resolved specparam models for
% short-time Fourier transformed signals.
% 
% The spectral parameterization algorithm used herein (specparam) can be 
% cited as:
%   Donoghue, T., Haller, M., Peterson, E.J. et al., Parameterizing neural 
%   power spectra into periodic and aperiodic components. Nat Neurosci 23,
%   1655–1665 (2020). https://doi.org/10.1038/s41593-020-00744-x
%
% Author: Luc Wilson (2022)

    fMask = logical(round(fs.*10)./10 >= round(opt.freq_range(1).*10)./10 & (round(fs.*10)./10 <= round(opt.freq_range(2).*10)./10));
    fs = fs(fMask);
    s_data.Freqs = fs;
    nChan = size(TF,1);
    nTimes = size(TF,2);
    % Adjust TF plots to only include modelled frequencies
    TF = TF(:,:,fMask);
    % Initalize FOOOF structs
    channel(nChan) = struct();
    SPRiNT = struct('options',opt,'freqs',fs,'channel',channel,'SPRiNT_models',nan(size(TF)),'peak_models',nan(size(TF)),'aperiodic_models',nan(size(TF)));
    % Iterate across channels
    for chan = 1:nChan
        channel(chan).data(nTimes) = struct(...
            'time',             [],...
            'aperiodic_params', [],...
            'peak_params',      [],...
            'peak_types',       '',...
            'ap_fit',           [],...
            'fooofed_spectrum', [],...
            'power_spectrum',   [],...
            'peak_fit',         [],...
            'error',            [],...
            'r_squared',        []);
        channel(chan).peaks(nTimes*opt.max_peaks) = struct(...
            'time',             [],...
            'center_frequency', [],...
            'amplitude',        [],...
            'st_dev',           []);
        channel(chan).aperiodics(nTimes) = struct(...
            'time',             [],...
            'offset',           [],...
            'exponent',         []);
        channel(chan).stats(nTimes) = struct(...
            'MSE',              [],...
            'r_squared',        [],...
            'frequency_wise_error', []);
        spec = log10(squeeze(TF(chan,:,:))); % extract log spectra for a given channel
        if nTimes == 1
            spec =reshape(spec,1,length(spec));
        end
        % Iterate across time
        i = 1; % For peak extraction
        ag = -(spec(1,end)-spec(1,1))./log10(fs(end)./fs(1)); % aperiodic guess initialization
        for time = 1:nTimes
            % Fit aperiodic 
            aperiodic_pars = robust_ap_fit(fs, spec(time,:), opt.aperiodic_mode, ag);
            % Remove aperiodic
            flat_spec = flatten_spectrum(fs, spec(time,:), aperiodic_pars, opt.aperiodic_mode);
            % Fit peaks
            [peak_pars, peak_function] = fit_peaks(fs, flat_spec, opt.max_peaks, opt.peak_threshold, opt.min_peak_height, ...
                opt.peak_width_limits/2, opt.proximity_threshold, opt.peak_type, opt.guess_weight,opt.hOT);
            if opt.thresh_after && ~opt.hOT  % Check thresholding requirements are met for unbounded optimization
                peak_pars(peak_pars(:,2) < opt.min_peak_height,:)     = []; % remove peaks shorter than limit
                peak_pars(peak_pars(:,3) < opt.peak_width_limits(1)/2,:)  = []; % remove peaks narrower than limit
                peak_pars(peak_pars(:,3) > opt.peak_width_limits(2)/2,:)  = []; % remove peaks broader than limit
                peak_pars = drop_peak_cf(peak_pars, opt.proximity_threshold, opt.freq_range); % remove peaks outside frequency limits
                peak_pars(peak_pars(:,1) < 0,:) = []; % remove peaks with a centre frequency less than zero (bypass drop_peak_cf)
                peak_pars = drop_peak_overlap(peak_pars, opt.proximity_threshold); % remove smallest of two peaks fit too closely
            end
            % Refit aperiodic
            aperiodic = spec(time,:);
            for peak = 1:size(peak_pars,1)
                aperiodic = aperiodic - peak_function(fs,peak_pars(peak,1), peak_pars(peak,2), peak_pars(peak,3));
            end
            aperiodic_pars = simple_ap_fit(fs, aperiodic, opt.aperiodic_mode, aperiodic_pars(end));
            ag = aperiodic_pars(end); % save exponent estimate for next iteration
            % Generate model fit
            ap_fit = gen_aperiodic(fs, aperiodic_pars, opt.aperiodic_mode);
            model_fit = ap_fit;
            for peak = 1:size(peak_pars,1)
                model_fit = model_fit + peak_function(fs,peak_pars(peak,1),...
                    peak_pars(peak,2),peak_pars(peak,3));
            end
            % Calculate model error
            MSE = sum((spec(time,:) - model_fit).^2)/length(model_fit);
            rsq_tmp = corrcoef(spec(time,:),model_fit).^2;
            % Return FOOOF results
            aperiodic_pars(2) = abs(aperiodic_pars(2));
            channel(chan).data(time).time                = ts(time);
            channel(chan).data(time).aperiodic_params    = aperiodic_pars;
            channel(chan).data(time).peak_params         = peak_pars;
            channel(chan).data(time).peak_types          = func2str(peak_function);
            channel(chan).data(time).ap_fit              = 10.^ap_fit;
            aperiodic_models(chan,time,:)                = 10.^ap_fit;
            channel(chan).data(time).fooofed_spectrum    = 10.^model_fit;
            SPRiNT_models(chan,time,:)                   = 10.^model_fit;
            channel(chan).data(time).power_spectrum   	 = 10.^spec(time,:);
            channel(chan).data(time).peak_fit            = 10.^(model_fit-ap_fit); 
            peak_models(chan,time,:)                     = 10.^(model_fit-ap_fit); 
            channel(chan).data(time).error               = MSE;
            channel(chan).data(time).r_squared           = rsq_tmp(2);
            % Extract peaks
            if ~isempty(peak_pars) & any(peak_pars)
                for p = 1:size(peak_pars,1)
                    channel(chan).peaks(i).time = ts(time);
                    channel(chan).peaks(i).center_frequency = peak_pars(p,1);
                    channel(chan).peaks(i).amplitude = peak_pars(p,2);
                    channel(chan).peaks(i).st_dev = peak_pars(p,3);
                    i = i +1;
                end
            end
            % Extract aperiodic
            channel(chan).aperiodics(time).time = ts(time);
            channel(chan).aperiodics(time).offset = aperiodic_pars(1);
            if length(aperiodic_pars)>2 % Legacy specparam alters order of parameters
                channel(chan).aperiodics(time).exponent = aperiodic_pars(3);
                channel(chan).aperiodics(time).knee_frequency = aperiodic_pars(2);
            else
                channel(chan).aperiodics(time).exponent = aperiodic_pars(2);
            end
            channel(chan).stats(time).MSE = MSE;
            channel(chan).stats(time).r_squared = rsq_tmp(2);
            channel(chan).stats(time).frequency_wise_error = abs(spec(time,:)-model_fit);
        end
        channel(chan).peaks(i:end) = [];
    end
    SPRiNT.channel = channel;
    SPRiNT.aperiodic_models = aperiodic_models;
    SPRiNT.SPRiNT_models = SPRiNT_models;
    SPRiNT.peak_models = peak_models;
    if strcmp(opt.rmoutliers,'yes')
        SPRiNT = remove_outliers(SPRiNT,peak_function,opt);
    end
    SPRiNT = cluster_peaks_dynamic2(SPRiNT); % Cluster peaks
    s_data.SPRiNT = SPRiNT;
end

function SPRiNT = remove_outliers(SPRiNT,peak_function,opt)
% Helper function to remove outlier peaks according to user-defined specifications 
% Author: Luc Wilson

    timeRange = opt.maxtime.*opt.WinLength.*(1-opt.WinOverlap./100);
    nC = length(SPRiNT.channel);
    for c = 1:nC
        ts = [SPRiNT.channel(c).data.time];
        remove = 1;
        while any(remove) 
            remove = zeros(length([SPRiNT.channel(c).peaks]),1);
            for p = 1:length([SPRiNT.channel(c).peaks])
                if sum((abs([SPRiNT.channel(c).peaks.time] - SPRiNT.channel(c).peaks(p).time) <= timeRange) &...
                        (abs([SPRiNT.channel(c).peaks.center_frequency] - SPRiNT.channel(c).peaks(p).center_frequency) <= opt.maxfreq)) < opt.minnear+1 % includes current peak
                    remove(p) = 1;
                end
            end
            SPRiNT.channel(c).peaks(logical(remove)) = [];
        end
        
        for t = 1:length(ts)
            
            if SPRiNT.channel(c).data(t).peak_params(1) == 0
                continue % never any peaks to begin with
            end
            p = [SPRiNT.channel(c).peaks.time] == ts(t);
            if sum(p) == size(SPRiNT.channel(c).data(t).peak_params,1)
                continue % number of peaks has not changed
            end
            peak_fit = zeros(size(SPRiNT.freqs));
            if any(p)
                SPRiNT.channel(c).data(t).peak_params = [[SPRiNT.channel(c).peaks(p).center_frequency]' [SPRiNT.channel(c).peaks(p).amplitude]' [SPRiNT.channel(c).peaks(p).st_dev]'];
                peak_pars = SPRiNT.channel(c).data(t).peak_params;
                for peak = 1:size(peak_pars,1)
                    peak_fit = peak_fit + peak_function(SPRiNT.freqs,peak_pars(peak,1),...
                        peak_pars(peak,2),peak_pars(peak,3));
                end
                ap_spec = log10(SPRiNT.channel(c).data(t).power_spectrum) - peak_fit;
                ap_pars = simple_ap_fit(SPRiNT.freqs, ap_spec, opt.aperiodic_mode, SPRiNT.channel(c).data(t).aperiodic_params(end));
                ap_fit = gen_aperiodic(SPRiNT.freqs, ap_pars, opt.aperiodic_mode);
                MSE = sum((ap_spec - ap_fit).^2)/length(SPRiNT.freqs);
                rsq_tmp = corrcoef(ap_spec+peak_fit,ap_fit+peak_fit).^2;
                % Return FOOOF results
                ap_pars(2) = abs(ap_pars(2));
                SPRiNT.channel(c).data(t).ap_fit = 10.^(ap_fit);
                SPRiNT.channel(c).data(t).fooofed_spectrum = 10.^(ap_fit+peak_fit);
                SPRiNT.channel(c).data(t).peak_fit = 10.^(peak_fit);
                SPRiNT.aperiodic_models(c,t,:) = SPRiNT.channel(c).data(t).ap_fit;
                SPRiNT.SPRiNT_models(c,t,:) = SPRiNT.channel(c).data(t).fooofed_spectrum;
                SPRiNT.peak_models(c,t,:) = SPRiNT.channel(c).data(t).peak_fit;
                SPRiNT.channel(c).aperiodics(t).offset = ap_pars(1);
                if length(ap_pars)>2 % Legacy FOOOF alters order of parameters
                    SPRiNT.channel(c).aperiodics(t).exponent = ap_pars(3);
                    SPRiNT.channel(c).aperiodics(t).knee_frequency = ap_pars(2);
                else
                    SPRiNT.channel(c).aperiodics(t).exponent = ap_pars(2);
                end
                SPRiNT.channel(c).stats(t).MSE = MSE;
                SPRiNT.channel(c).stats(t).r_squared = rsq_tmp(2);
                SPRiNT.channel(c).stats(t).frequency_wise_error = abs(ap_spec-ap_fit);
                
            else
                SPRiNT.channel(c).data(t).peak_params = [0 0 0];
                ap_spec = log10(SPRiNT.channel(c).data(t).power_spectrum) - peak_fit;
                ap_pars = simple_ap_fit(SPRiNT.freqs, ap_spec, opt.aperiodic_mode, SPRiNT.channel(c).data(t).aperiodic_params(end));
                ap_fit = gen_aperiodic(SPRiNT.freqs, ap_pars, opt.aperiodic_mode);
                MSE = sum((ap_spec - ap_fit).^2)/length(SPRiNT.freqs);
                rsq_tmp = corrcoef(ap_spec+peak_fit,ap_fit+peak_fit).^2;
                % Return FOOOF results
                ap_pars(2) = abs(ap_pars(2));
                SPRiNT.channel(c).data(t).ap_fit = 10.^(ap_fit);
                SPRiNT.channel(c).data(t).fooofed_spectrum = 10.^(ap_fit+peak_fit);
                SPRiNT.channel(c).data(t).peak_fit = 10.^(peak_fit);
                SPRiNT.aperiodic_models(c,t,:) = SPRiNT.channel(c).data(t).ap_fit;
                SPRiNT.SPRiNT_models(c,t,:) = SPRiNT.channel(c).data(t).fooofed_spectrum;
                SPRiNT.peak_models(c,t,:) = SPRiNT.channel(c).data(t).peak_fit;
                SPRiNT.channel(c).aperiodics(t).offset = ap_pars(1);
                if length(ap_pars)>2 % Legacy FOOOF alters order of parameters
                    SPRiNT.channel(c).aperiodics(t).exponent = ap_pars(3);
                    SPRiNT.channel(c).aperiodics(t).knee_frequency = ap_pars(2);
                else
                    SPRiNT.channel(c).aperiodics(t).exponent = ap_pars(2);
                end
                SPRiNT.channel(c).stats(t).MSE = MSE;
                SPRiNT.channel(c).stats(t).r_squared = rsq_tmp(2);
                SPRiNT.channel(c).stats(t).frequency_wise_error = abs(ap_spec-ap_fit);
            end
        end
    end
end

function oS = cluster_peaks_dynamic2(oS)
% Second Helper function to cluster peaks across time
% Author: Luc Wilson (2022)
    pthr = oS.options.proximity_threshold;
    for chan = 1:length(oS.channel)
        clustLead = [];
        nCl = 0;
        oS.channel(chan).clustered_peaks = struct();
        times = unique([oS.channel(chan).peaks.time]);
        all_peaks = oS.channel(chan).peaks;
        for time = 1:length(times)
            time_peaks = all_peaks([all_peaks.time] == times(time));
            % Initialize first clusters
            if time == 1
                nCl = length(time_peaks);
                for Cl = 1:nCl
                    oS.channel(chan).clustered_peaks(Cl).cluster = Cl;
                    oS.channel(chan).clustered_peaks(Cl).peaks(Cl) = time_peaks(Cl);
                    clustLead(Cl,1) = time_peaks(Cl).time;
                    clustLead(Cl,2) = time_peaks(Cl).center_frequency;
                    clustLead(Cl,3) = time_peaks(Cl).amplitude;
                    clustLead(Cl,4) = time_peaks(Cl).st_dev;
                    clustLead(Cl,5) = Cl;
                end
                continue
            end
            
            % Cluster "drafting stage": find points that make good matches to each cluster.  
            for Cl = 1:nCl
                match = abs([time_peaks.center_frequency]-clustLead(Cl,2))./clustLead(Cl,4) < pthr;
                idx_tmp = find(match);
                if any(match)
                    % Add the best candidate peak 
                    % Note: Auto-adds only peaks, but adds best candidate
                    % for multiple options
                    [tmp,idx] = min(([time_peaks(match).center_frequency] - clustLead(Cl,2)).^2 +...
                            ([time_peaks(match).amplitude] - clustLead(Cl,3)).^2 +...
                            ([time_peaks(match).st_dev] - clustLead(Cl,4)).^2);
                    oS.channel(chan).clustered_peaks(clustLead(Cl,5)).peaks(length(oS.channel(chan).clustered_peaks(clustLead(Cl,5)).peaks)+1) = time_peaks(idx_tmp(idx)); 
                    clustLead(Cl,1) = time_peaks(idx_tmp(idx)).time;
                    clustLead(Cl,2) = time_peaks(idx_tmp(idx)).center_frequency;
                    clustLead(Cl,3) = time_peaks(idx_tmp(idx)).amplitude;
                    clustLead(Cl,4) = time_peaks(idx_tmp(idx)).st_dev;
                    % Don't forget to remove the candidate from the pool
                    time_peaks(idx_tmp(idx)) = [];
                end
            end
            % Remaining peaks get sorted into their own clusters
            if ~isempty(time_peaks)
                for peak = 1:length(time_peaks)
                    nCl = nCl + 1;
                    Cl = nCl;
                    clustLead(Cl,1) = time_peaks(peak).time;
                    clustLead(Cl,2) = time_peaks(peak).center_frequency;
                    clustLead(Cl,3) = time_peaks(peak).amplitude;
                    clustLead(Cl,4) = time_peaks(peak).st_dev;
                    clustLead(Cl,5) = Cl;
                    oS.channel(chan).clustered_peaks(Cl).cluster = Cl;
                    oS.channel(chan).clustered_peaks(Cl).peaks(length(oS.channel(chan).clustered_peaks(clustLead(Cl,5)).peaks)+1) = time_peaks(peak); 
                end
            end     
            % Sort clusters based on most recent
            clustLead = sortrows(clustLead,1,'descend');
        end
    end
end

%% ===== GENERATE APERIODIC =====
function ap_vals = gen_aperiodic(freqs,aperiodic_params,aperiodic_mode)
%       Generate aperiodic values, from parameter definition.
%
%       Parameters
%       ----------
%       freqs : 1xn array
%       	Frequency vector to create aperiodic component for.
%       aperiodic_params : 1x3 array
%           Parameters that define the aperiodic component.
%       aperiodic_mode : {'fixed', 'knee'}
%           Defines absence or presence of knee in aperiodic component.
%
%       Returns
%       -------
%       ap_vals : 1d array
%           Generated aperiodic values.

    switch aperiodic_mode
        case 'fixed'  % no knee
            ap_vals = expo_nk_function(freqs,aperiodic_params);
        case 'knee'
            ap_vals = expo_function(freqs,aperiodic_params);
    end
end


%% ===== CORE MODELS =====
function ys = gaussian(freqs, mu, hgt, sigma)
%       Gaussian function to use for fitting.
%
%       Parameters
%       ----------
%       freqs : 1xn array
%           Frequency vector to create gaussian fit for.
%       mu, hgt, sigma : doubles
%           Parameters that define gaussian function (centre frequency,
%           height, and standard deviation).
%
%       Returns
%       -------
%       ys :    1xn array
%       Output values for gaussian function.

    ys = hgt*exp(-(((freqs-mu)./sigma).^2) /2);

end

function ys = cauchy(freqs, ctr, hgt, gam)
%       Cauchy function to use for fitting.
% 
%       Parameters
%       ----------
%       freqs : 1xn array
%           Frequency vector to create cauchy fit for.
%       ctr, hgt, gam : doubles
%           Parameters that define cauchy function (centre frequency,
%           height, and "standard deviation" [gamma]).
%
%       Returns
%       -------
%       ys :    1xn array
%       Output values for cauchy function.

    ys = hgt./(1+((freqs-ctr)/gam).^2);

end

function ys = expo_function(freqs,params)
%       Exponential function to use for fitting 1/f, with a 'knee' (maximum at low frequencies).
%
%       Parameters
%       ----------
%       freqs : 1xn array
%           Input x-axis values.
%       params : 1x3 array (offset, knee, exp)
%           Parameters (offset, knee, exp) that define Lorentzian function:
%           y = 10^offset * (1/(knee + x^exp))
%
%       Returns
%       -------
%       ys :    1xn array
%           Output values for exponential function.

    ys = params(1) - log10(abs(params(2)) +freqs.^params(3));

end

function ys = expo_nk_function(freqs, params)
%       Exponential function to use for fitting 1/f, without a 'knee'.
%
%       Parameters
%       ----------
%       freqs : 1xn array
%           Input x-axis values.
%       params : 1x2 array (offset, exp)
%           Parameters (offset, exp) that define Lorentzian function:
%           y = 10^offset * (1/(x^exp))
%
%       Returns
%       -------
%       ys :    1xn array
%           Output values for exponential (no-knee) function.

    ys = params(1) - log10(freqs.^params(2));

end

%% ===== FITTING ALGORITHM =====
function aperiodic_params = simple_ap_fit(freqs, power_spectrum, aperiodic_mode, aperiodic_guess)
%       Fit the aperiodic component of the power spectrum.
%
%       Parameters
%       ----------
%       freqs : 1xn array
%           Frequency values for the power spectrum, in linear scale.
%       power_spectrum : 1xn array
%           Power values, in log10 scale.
%       aperiodic_mode : {'fixed','knee'}
%           Defines absence or presence of knee in aperiodic component.
%       aperiodic_guess: double
%           SPRiNT specific - feeds previous timepoint aperiodic slope as
%           guess
%
%       Returns
%       -------
%       aperiodic_params : 1xn array
%           Parameter estimates for aperiodic fit.

%       Set guess params for lorentzian aperiodic fit, guess params set at init
    options = optimset('Display', 'off', 'TolX', 1e-4, 'TolFun', 1e-6, ...
        'MaxFunEvals', 5000, 'MaxIter', 5000);

    switch (aperiodic_mode)
        case 'fixed'  % no knee
            guess_vec = [power_spectrum(1), aperiodic_guess];
            aperiodic_params = fminsearch(@error_expo_nk_function, guess_vec, options, freqs, power_spectrum);
        case 'knee'
            guess_vec = [power_spectrum(1),0, aperiodic_guess];
            aperiodic_params = fminsearch(@error_expo_function, guess_vec, options, freqs, power_spectrum);
    end

end

function aperiodic_params = robust_ap_fit(freqs, power_spectrum, aperiodic_mode, aperiodic_guess)
%       Fit the aperiodic component of the power spectrum robustly, ignoring outliers.
%
%       Parameters
%       ----------
%       freqs : 1xn array
%           Frequency values for the power spectrum, in linear scale.
%       power_spectrum : 1xn array
%           Power values, in log10 scale.
%       aperiodic_mode : {'fixed','knee'}
%           Defines absence or presence of knee in aperiodic component.
%       aperiodic_guess: double
%           SPRiNT specific - feeds previous timepoint aperiodic slope as
%           guess
%
%       Returns
%       -------
%       aperiodic_params : 1xn array
%           Parameter estimates for aperiodic fit.

    % Do a quick, initial aperiodic fit
    popt = simple_ap_fit(freqs, power_spectrum, aperiodic_mode, aperiodic_guess);
    initial_fit = gen_aperiodic(freqs, popt, aperiodic_mode);

    % Flatten power_spectrum based on initial aperiodic fit
    flatspec = power_spectrum - initial_fit;

    % Flatten outliers - any points that drop below 0
    flatspec(flatspec(:) < 0) = 0;

    % Use percential threshold, in terms of # of points, to extract and re-fit
    perc_thresh = prctile(flatspec, 0.025);
    perc_mask = flatspec <= perc_thresh;
    freqs_ignore = freqs(perc_mask);
    spectrum_ignore = power_spectrum(perc_mask);

    % Second aperiodic fit - using results of first fit as guess parameters

    options = optimset('Display', 'off', 'TolX', 1e-4, 'TolFun', 1e-6, ...
        'MaxFunEvals', 5000, 'MaxIter', 5000);
    guess_vec = popt;

    switch (aperiodic_mode)
        case 'fixed'  % no knee
            aperiodic_params = fminsearch(@error_expo_nk_function, guess_vec, options, freqs_ignore, spectrum_ignore);
        case 'knee'
            aperiodic_params = fminsearch(@error_expo_function, guess_vec, options, freqs_ignore, spectrum_ignore);
    end
end

function spectrum_flat = flatten_spectrum(freqs, power_spectrum, robust_aperiodic_params, aperiodic_mode)
%       Flatten the power spectrum by removing the aperiodic component.
%
%       Parameters
%       ----------
%       freqs : 1xn array
%           Frequency values for the power spectrum, in linear scale.
%       power_spectrum : 1xn array
%           Power values, in log10 scale.
%       robust_aperiodic_params : 1x2 or 1x3 array (see aperiodic_mode)
%           Parameter estimates for aperiodic fit.
%       aperiodic_mode : 1 or 2
%           Defines absence or presence of knee in aperiodic component.
%
%       Returns
%       -------
%       spectrum_flat : 1xn array
%           Flattened (aperiodic removed) power spectrum.


spectrum_flat = power_spectrum - gen_aperiodic(freqs,robust_aperiodic_params,aperiodic_mode);

end

function [model_params,peak_function] = fit_peaks(freqs, flat_iter, max_n_peaks, peak_threshold, min_peak_height, gauss_std_limits, proxThresh, peakType, guess_weight, hOT)
%       Iteratively fit peaks to flattened spectrum.
%
%       Parameters
%       ----------
%       freqs : 1xn array
%           Frequency values for the power spectrum, in linear scale.
%       flat_iter : 1xn array
%           Flattened (aperiodic removed) power spectrum.
%       max_n_peaks : double
%           Maximum number of gaussians to fit within the spectrum.
%       peak_threshold : double
%           Threshold (in standard deviations of noise floor) to detect a peak.
%       min_peak_height : double
%           Minimum height of a peak (in log10).
%       gauss_std_limits : 1x2 double
%           Limits to gaussian (cauchy) standard deviation (gamma) when detecting a peak.
%       proxThresh : double
%           Minimum distance between two peaks, in st. dev. (gamma) of peaks.
%       peakType : {'gaussian', 'cauchy'}
%           Which types of peaks are being fitted
%       guess_weight : {'none', 'weak', 'strong'}
%           Parameter to weigh initial estimates during optimization (None, Weak, or Strong)
%       hOT : 0 or 1
%           Defines whether to use constrained optimization, fmincon, or
%           basic simplex, fminsearch.
%
%       Returns
%       -------
%       gaussian_params : mx3 array, where m = No. of peaks.
%           Parameters that define the peak fit(s). Each row is a peak, as [mean, height, st. dev. (gamma)].
    switch peakType 
        case 'gaussian' % gaussian only
            peak_function = @gaussian; % Identify peaks as gaussian
            % Initialize matrix of guess parameters for gaussian fitting.
            guess_params = zeros(max_n_peaks, 3);
            % Save intact flat_spectrum
            flat_spec = flat_iter;
            % Find peak: Loop through, finding a candidate peak, and fitting with a guess gaussian.
            % Stopping procedure based on either the limit on # of peaks,
            % or the relative or absolute height thresholds.
            for guess = 1:max_n_peaks
                % Find candidate peak - the maximum point of the flattened spectrum.
                max_ind = find(flat_iter == max(flat_iter));
                max_height = flat_iter(max_ind);

                % Stop searching for peaks once max_height drops below height threshold.
                if max_height <= peak_threshold * std(flat_iter)
                    break
                end

                % Set the guess parameters for gaussian fitting - mean and height.
                guess_freq = freqs(max_ind);
                guess_height = max_height;

                % Halt fitting process if candidate peak drops below minimum height.
                if guess_height <= min_peak_height
                    break
                end

                % Data-driven first guess at standard deviation
                % Find half height index on each side of the center frequency.
                half_height = 0.5 * max_height;

                le_ind = sum(flat_iter(1:max_ind) <= half_height);
                ri_ind = length(flat_iter) - sum(flat_iter(max_ind:end) <= half_height)+1;

                % Keep bandwidth estimation from the shortest side.
                % We grab shortest to avoid estimating very large std from overalapping peaks.
                % Grab the shortest side, ignoring a side if the half max was not found.
                % Note: will fail if both le & ri ind's end up as None (probably shouldn't happen).
                short_side = min(abs([le_ind,ri_ind]-max_ind));

                % Estimate std from FWHM. Calculate FWHM, converting to Hz, get guess std from FWHM
                fwhm = short_side * 2 * (freqs(2)-freqs(1));
                guess_std = fwhm / (2 * sqrt(2 * log(2)));

                % Check that guess std isn't outside preset std limits; restrict if so.
                % Note: without this, curve_fitting fails if given guess > or < bounds.
                if guess_std < gauss_std_limits(1)
                    guess_std = gauss_std_limits(1);
                end
                if guess_std > gauss_std_limits(2)
                    guess_std = gauss_std_limits(2);
                end

                % Collect guess parameters.
                guess_params(guess,:) = [guess_freq, guess_height, guess_std];

                % Subtract best-guess gaussian.
                peak_gauss = gaussian(freqs, guess_freq, guess_height, guess_std);
                flat_iter = flat_iter - peak_gauss;

            end
            % Remove unused guesses
            guess_params(guess_params(:,1) == 0,:) = [];

            % Check peaks based on edges, and on overlap
            % Drop any that violate requirements.
            guess_params = drop_peak_cf(guess_params, proxThresh, [min(freqs) max(freqs)]);
            guess_params = drop_peak_overlap(guess_params, proxThresh);

            % If there are peak guesses, fit the peaks, and sort results.
            if ~isempty(guess_params)
                model_params = fit_peak_guess(guess_params, freqs, flat_spec, 1, guess_weight, gauss_std_limits,hOT);
            else
                model_params = zeros(1, 3);
            end
            
        case 'cauchy' % cauchy only
            peak_function = @cauchy; % Identify peaks as cauchy
            guess_params = zeros(max_n_peaks, 3);
            flat_spec = flat_iter;
            for guess = 1:max_n_peaks
                max_ind = find(flat_iter == max(flat_iter));
                max_height = flat_iter(max_ind);
                if max_height <= peak_threshold * std(flat_iter)
                    break
                end
                guess_freq = freqs(max_ind);
                guess_height = max_height;
                if guess_height <= min_peak_height
                    break
                end
                half_height = 0.5 * max_height;
                le_ind = sum(flat_iter(1:max_ind) <= half_height);
                ri_ind = length(flat_iter) - sum(flat_iter(max_ind:end) <= half_height);
                short_side = min(abs([le_ind,ri_ind]-max_ind));

                % Estimate gamma from FWHM. Calculate FWHM, converting to Hz, get guess gamma from FWHM
                fwhm = short_side * 2 * (freqs(2)-freqs(1));
                guess_gamma = fwhm/2;
                % Check that guess gamma isn't outside preset limits; restrict if so.
                % Note: without this, curve_fitting fails if given guess > or < bounds.
                if guess_gamma < gauss_std_limits(1)
                    guess_gamma = gauss_std_limits(1);
                end
                if guess_gamma > gauss_std_limits(2)
                    guess_gamma = gauss_std_limits(2);
                end

                % Collect guess parameters.
                guess_params(guess,:) = [guess_freq(1), guess_height, guess_gamma];

                % Subtract best-guess cauchy.
                peak_cauchy = cauchy(freqs, guess_freq(1), guess_height, guess_gamma);
                flat_iter = flat_iter - peak_cauchy;

            end
            guess_params(guess_params(:,1) == 0,:) = [];
            guess_params = drop_peak_cf(guess_params, proxThresh, [min(freqs) max(freqs)]);
            guess_params = drop_peak_overlap(guess_params, proxThresh);

            % If there are peak guesses, fit the peaks, and sort results.
            if ~isempty(guess_params)
                model_params = fit_peak_guess(guess_params, freqs, flat_spec, 2, guess_weight, gauss_std_limits,hOT);
            else
                model_params = zeros(1, 3);
            end
    end
            
end

function guess = drop_peak_cf(guess, bw_std_edge, freq_range)
%       Check whether to drop peaks based on center's proximity to the edge of the spectrum.
%
%       Parameters
%       ----------
%       guess : mx3 array, where m = No. of peaks.
%           Guess parameters for peak fits.
%
%       Returns
%       -------
%       guess : qx3 where q <= m No. of peaks.
%           Guess parameters for peak fits.

    cf_params = guess(:,1)';
    bw_params = guess(:,3)' * bw_std_edge;

    % Check if peaks within drop threshold from the edge of the frequency range.

    keep_peak = abs(cf_params-freq_range(1)) > bw_params & ...
        abs(cf_params-freq_range(2)) > bw_params;

    % Drop peaks that fail the center frequency edge criterion
    guess = guess(keep_peak,:);

end

function guess = drop_peak_overlap(guess, proxThresh)
%       Checks whether to drop gaussians based on amount of overlap.
%
%       Parameters
%       ----------
%       guess : mx3 array, where m = No. of peaks.
%           Guess parameters for peak fits.
%       proxThresh: double
%           Proximity threshold (in st. dev. or gamma) between two peaks.
%
%       Returns
%       -------
%       guess : qx3 where q <= m No. of peaks.
%           Guess parameters for peak fits.
%
%       Note
%       -----
%       For any gaussians with an overlap that crosses the threshold,
%       the lowest height guess guassian is dropped.

    % Sort the peak guesses, so can check overlap of adjacent peaks
    guess = sortrows(guess);

    % Calculate standard deviation bounds for checking amount of overlap

    bounds = [guess(:,1) - guess(:,3) * proxThresh, ...
        guess(:,1), guess(:,1) + guess(:,3) * proxThresh];

    % Loop through peak bounds, comparing current bound to that of next peak
    drop_inds =  [];

    for ind = 1:size(bounds,1)-1

        b_0 = bounds(ind,:);
        b_1 = bounds(ind + 1,:);

        % Check if bound of current peak extends into next peak
        if b_0(2) > b_1(1)
            % If so, get the index of the gaussian with the lowest height (to drop)
            drop_inds = [drop_inds (ind - 1 + find(guess(ind:ind+1,2) == ...
                min(guess(ind,2),guess(ind+1,2))))];
        end
    end
    % Drop any peaks guesses that overlap too much, based on threshold.
    guess(drop_inds,:) = [];
end

function peak_params = fit_peak_guess(guess, freqs, flat_spec, peak_type, guess_weight, std_limits, hOT)
%     Fits a group of peak guesses with a fit function.
%
%     Parameters
%     ----------
%       guess : mx3 array, where m = No. of peaks.
%           Guess parameters for peak fits.
%       freqs : 1xn array
%           Frequency values for the power spectrum, in linear scale.
%       flat_iter : 1xn array
%           Flattened (aperiodic removed) power spectrum.
%       peakType : {'gaussian', 'cauchy', 'best'}
%           Which types of peaks are being fitted.
%       guess_weight : 'none', 'weak', 'strong'
%           Parameter to weigh initial estimates during optimization.
%       std_limits: 1x2 array
%           Minimum and maximum standard deviations for distribution.
%       hOT : 0 or 1
%           Defines whether to use constrained optimization, fmincon, or
%           basic simplex, fminsearch.
%
%       Returns
%       -------
%       peak_params : mx3, where m =  No. of peaks.
%           Peak parameters post-optimization.

    
    if hOT % Use OptimToolbox for fmincon
        options = optimset('Display', 'off', 'TolX', 1e-3, 'TolFun', 1e-5, ...
        'MaxFunEvals', 3000, 'MaxIter', 3000); % Tuned options
        lb = [guess(:,1)-guess(:,3)*2,zeros(size(guess(:,2))),ones(size(guess(:,3)))*std_limits(1)];
        ub = [guess(:,1)+guess(:,3)*2,inf(size(guess(:,2))),ones(size(guess(:,3)))*std_limits(2)];
        peak_params = fmincon(@error_model_constr,guess,[],[],[],[], ...
            lb,ub,[],options,freqs,flat_spec, peak_type);
    else % Use basic simplex approach, fminsearch, with guess_weight
        options = optimset('Display', 'off', 'TolX', 1e-4, 'TolFun', 1e-5, ...
        'MaxFunEvals', 5000, 'MaxIter', 5000);
        peak_params = fminsearch(@error_model,...
            guess, options, freqs, flat_spec, peak_type, guess, guess_weight);
    end
end


%% ===== ERROR FUNCTIONS =====
function err = error_expo_nk_function(params,xs,ys)
    ym = -log10(xs.^params(2)) + params(1);
    err = sum((ys - ym).^2);
end

function err = error_expo_function(params,xs,ys)
    ym = expo_function(xs,params);
    err = sum((ys - ym).^2);
end

function err = error_model(params, xVals, yVals, peak_type, guess, guess_weight)
    fitted_vals = 0;
    weak = 1E2;
    strong = 1E7;
    for set = 1:size(params,1)
        switch (peak_type)
            case 1 % Gaussian
                fitted_vals = fitted_vals + gaussian(xVals, params(set,1), params(set,2), params(set,3));
            case 2 % Cauchy
                fitted_vals = fitted_vals + cauchy(xVals, params(set,1), params(set,2), params(set,3));
        end
    end
    switch guess_weight
        case 'none'
            err = sum((yVals - fitted_vals).^2);
        case 'weak' % Add small weight to deviations from guess m and amp
            err = sum((yVals - fitted_vals).^2) + ...
                 weak*sum((params(:,1)-guess(:,1)).^2) + ...
                 weak*sum((params(:,2)-guess(:,2)).^2);
        case 'strong' % Add large weight to deviations from guess m and amp
            err = sum((yVals - fitted_vals).^2) + ...
                 strong*sum((params(:,1)-guess(:,1)).^2) + ...
                 strong*sum((params(:,2)-guess(:,2)).^2);
    end
end

function err = error_model_constr(params, xVals, yVals, peak_type)
    fitted_vals = 0;
    for set = 1:size(params,1)
        switch (peak_type)
            case 1 % Gaussian
                fitted_vals = fitted_vals + gaussian(xVals, params(set,1), params(set,2), params(set,3));
            case 2 % Cauchy
                fitted_vals = fitted_vals + cauchy(xVals, params(set,1), params(set,2), params(set,3));
        end
    end
    err = sum((yVals - fitted_vals).^2);
end


