function [fractal,original] = ft_fooof_matlab_v1(signal_data,param)
%FT_FOOOF_MATLAB_V1 Summary of this function goes here
%   Detailed explanation goes here
% signal_data should be 1 by time matrix
srate = param.srate;
slide_window_length_sec = param.slide_window_length_sec;
window_overlap = param.window_overlap;
freq_range = param.freq_range;

%% fooof in fieldtrip
% raw_data'
data = [];
n_tp = length(signal_data);
data.trial{1} = signal_data;
data.time{1,1}  = (1:n_tp)/srate;
data.label{1} = 'chan';
data.trialinfo(1,1) = 1;

%%
% chunk into 2-second segments
cfg               = [];
cfg.length        = slide_window_length_sec;
cfg.overlap       = window_overlap;
data              = ft_redefinetrial(cfg, data);

% compute the fractal and original spectra
cfg               = [];
cfg.foilim        = freq_range;
cfg.pad           = 4;
cfg.tapsmofrq     = 2;
cfg.method        = 'mtmfft';
cfg.output        = 'fooof_aperiodic';
fractal = ft_freqanalysis(cfg, data);
cfg.output        = 'pow';
original = ft_freqanalysis(cfg, data);

end

