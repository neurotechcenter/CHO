clear;
clc;


load data_for_figure8.mat
srate = 500;


%%
param.plot = 1;
param.minimum_cycles = 2;
param.ovlp_threshold = 0.50;
param.frequency_axis = 1:40;

[cho_output] = CHO_v22(raw_data, srate, param);
% [cho_output] = CHO_v25(raw_data, srate, param);


% %%
% target_freq = 50;
% filtered_data1 = wavelet_single_freq_filtering_v4(raw_data, srate, target_freq, 8);
% filtered_data2 = wavelet_single_freq_filtering_v3(raw_data, srate, target_freq, 8);
% 
% %%
% figure();
% plot(raw_data);
% hold on;
% plot(filtered_data1);
% plot(filtered_data2);
% 
% 
% 
% %%
% Fs = 500;            % Sampling frequency                    
% T = 1/Fs;             % Sampling period       
% L = 1500;             % Length of signal
% t = (0:L-1)*T;        % Time vector
% 
% S = 0.7*sin(2*pi*11*t) + 0.9*sin(2*pi*27*t);
% 
% raw_data = S + 0.1*randn(size(t));
% 
% figure();
% plot(1000*t,raw_data)
% title("Signal Corrupted with Zero-Mean Random Noise")
% xlabel("t (milliseconds)")
% ylabel("X(t)")
% 
% 
% %%
% target_freq = 11
% filtered_data0 = butter_bandpass_filtering(raw_data,srate,[10 12],3);
% filtered_data1 = wavelet_single_freq_filtering_v4(raw_data, srate, target_freq, 8);
% filtered_data2 = wavelet_single_freq_filtering_v3(raw_data, srate, target_freq, 8);
% 
% 
% %%
% figure();
% plot(raw_data);
% hold on;
% plot(filtered_data1);
% plot(filtered_data2);