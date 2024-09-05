clear;
clc;
% close all;


load example_sinusoidal_data_v1.mat
%%


srate = 400;
timex = linspace(-1500,0,srate*1.5);
param.plot = 1;
param.minimum_cycles = 2;
param.ovlp_threshold = 0.50;
param.frequency_axis = 1:40;

[cho_outputs] = CHO_v22(raw_data, srate, param);


%%
[pxx, f] = pwelch(raw_data',[],srate/2,1:srate/2,srate);
frequency_range = [1:40];
X = f(frequency_range)';
Y = log10(pxx(frequency_range))';
[gaussians,outputs] = fooof_matlab_v9(Y,X);

%% example and power spectrum 
figure('position',[10 10 1000 300]);
subaxis(1,2,1);
plot(timex,raw_data);
ylim([-250 250]);

subaxis(1,2,2);
plot(f,log10(pxx));
xlim([1 40]);

hold on;
plot(frequency_range,outputs.fooofed_ap_fit);
plot(frequency_range,outputs.final_fit);

legend('Power spectrum','Multi-gaussian fitting','1/f noise');

% plot(original.freq,10*log10(original.powspctrm));
% hold on;
% plot(fractal.freq,10*log10(fractal.powspctrm));

%% 
bandpassed_data{1} = fft_bandpass_filtering(raw_data,srate,[4 8]);
bandpassed_data{2} = fft_bandpass_filtering(raw_data,srate,[10 14]);
bandpassed_data{3} = fft_bandpass_filtering(raw_data,srate,[16 20]);

for i = 1:3
    hilbert_data{i} = hilbert(bandpassed_data{i});
    ag_data{i} = angle(hilbert_data{i});
    
    cos_data_1vs{i} = cos(i*ag_data{1});
    hilbert_1vs{i} = hilbert(cos_data_1vs{i});
    ag_data_1vs{i} = angle(hilbert_1vs{i});
end


n = 4;
figure('position',[10 10 400 1000]);
subplot(n,1,1);
plot(timex,raw_data);
ylim([-250 200]);

subplot(n,1,2);
plot(timex, bandpassed_data{1});
ylim([-250 200]);

subplot(n,1,3);
plot(timex, bandpassed_data{2});
ylim([-250 200]);
hold on;
plot(timex, cos_data_1vs{2}*20);

subplot(n,1,4);
plot(timex, bandpassed_data{3});
ylim([-250 200]);

hold on;
plot(timex, cos_data_1vs{3}*20);

%%
n = 4;
figure('position',[10 10 400 1000]);
subplot(n,1,1);
plot(timex,raw_data);
ylim([-250 200]);

subplot(n,1,2);
plot(timex, ag_data{1});

subplot(n,1,3);
plot(timex, ag_data{2});
hold on;
plot(timex, ag_data_1vs{2});

subplot(n,1,4);
plot(timex, ag_data{3});
hold on;
plot(timex, ag_data_1vs{3});

%% Sawtooth example to real data

dt = 1/srate;
T = 1.5; % total simulation time
t = dt:dt:T; % time vector

fundamental_band = eegfilt(raw_data',srate,4,8,0,150);
% harmonic_band = eegfilt(raw_data',srate,10,14,0,150); % theta and alpha
harmonic_band = eegfilt(raw_data',srate,16,20,0,150); % theta and beta

fundamental_phase = angle(hilbert(fundamental_band));
harmonic_phase = angle(hilbert(harmonic_band));

%% Surrogate Analysis

EpochLength = 0.5;

N_sample = 300;
nm = 1:10;

Roriginal = zeros(N_sample,length(nm));
Rsurr = zeros(N_sample,length(nm));

intervals = linspace(1,EpochLength/2,EpochLength/4);

for n_sample = 1:N_sample
    
    display(n_sample)
    
%     temp = randi(EpochLength/2);
%     I=temp:temp+round(EpochLength/2);
%     
%     temp2 = randi(EpochLength/2);
%     I2= temp2:temp2+round(EpochLength/2);

%     I = randi(EpochLength,1,EpochLength/2);
%     I2 = randi(EpochLength,1,EpochLength/2);

    temp = randi(round((T-EpochLength)*srate));
    I=temp:temp+round(EpochLength*srate);
    
    temp2 = randi(round((T-EpochLength)*srate));
    I2= temp2:temp2+round(EpochLength*srate);

%     I = round(intervals + randi(EpochLength/2));
%     I2 = round(intervals + randi(EpochLength/2));
%     
    for m = nm
        Roriginal(n_sample,m) = abs(mean(exp(1i*(unwrap(harmonic_phase(I))-m*unwrap(fundamental_phase(I))))));
        Rsurr(n_sample,m) = abs(mean(exp(1i*(harmonic_phase(I)-m*fundamental_phase(I2)))));
    end
end


%%
figure();
subplot(221)
plot(t,raw_data/max(raw_data),'k')
hold on
plot(t,fundamental_band/max(fundamental_band)-2,'b')
plot(t,harmonic_band/max(harmonic_band)-4,'r')
% plot(t,0.25*theta_phase/pi-5.5,'b.')
% plot(t,0.25*gamma_phase/pi-7,'r.')
hold off
% xlim([10 10.5])
axis off

subplot(222)
plot(t,fundamental_phase/max(fundamental_phase),'b.')
hold on
plot(t,harmonic_phase/max(harmonic_phase)-3,'r.')
hold off
% xlim([10 10.5])
axis off

subplot(2,3,[4 5])
plot(mean(Roriginal),'g-o','markerfacecolor','g')%,'markersize',10)
hold on
plot(mean(Rsurr),'r-o','markerfacecolor','r') %,'markersize',10)
hold off
ylim([0 1])
xlim([0 max(nm)])
box off

subplot(2,3,6)
m = 3;
boxplot([Roriginal(:,m),Rsurr(:,m)],'symbol','')

% Unpaired t-test to check if Original values are higher than Single Run Random
% Permutation Surrogates
[h p]=ttest2(Roriginal(:,m),Rsurr(:,m),'Tail','right','Alpha',0.01);
% title(p)
ylim([0 1])
set(gca,'ytick',[0:0.35:0.7])
box off


%%
[acf,lags,bounds] = autocorr(raw_data,'numlags',srate);

n = 2;
figure('position',[10 10 400 600]);
subplot(n,1,1);
plot(timex,raw_data);
ylim([-250 200]);

subplot(n,1,2);
plot(lags, acf);
ylim([-1 1]);