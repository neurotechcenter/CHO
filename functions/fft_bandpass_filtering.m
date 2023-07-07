function [ filtered_data ] = fft_bandpass_filtering(data, srate, fband )
%FFT_BANDPASS_FILTERING 이 함수의 요약 설명 위치
%   자세한 설명 위치

% data : time by channel or trials matrix
% srate : sampling rate
% fband : [low_bound high_bound]

n_tp = size(data,1);
n_ch = size(data,2);
Fs = srate;            % Sampling frequency
T = 1/Fs;             % Sampling period
L = n_tp;             % Length of signal
t = (0:L-1)*T;        % Time vector
f = Fs*(0:(L/2))/L;

X = data;
Y = fft(X);

f1 = fband(1);
f2 = fband(2);

target_fidx = find(f1 <= f & f<= f2);
nontarget_fidx = setdiff(1:length(f),target_fidx);

Y(nontarget_fidx,:) = 0;
Y(end-nontarget_fidx+1,:) = 0;

% Y(nontarget_fidx,:) = Y(nontarget_fidx,:)*(1/(10^100));
% Y(end-nontarget_fidx+1,:) = Y(end-nontarget_fidx+1,:)*(1/(10^100));

filtered_data = real(ifft(Y));

end

