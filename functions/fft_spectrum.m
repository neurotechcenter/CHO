function [power_spectrum, f] = fft_spectrum(data, srate)
% data: time by channel
Fs = srate;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = size(data,1);             % Length of signal
t = (0:L-1)*T;        % Time vector
X = data;
Y = fft(X);
f = Fs*(0:(L/2))/L;
power_spectrum = abs(Y(1:length(f),:));
end