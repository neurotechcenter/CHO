clear;
clc;
close all;



%%
warning('off');

load example_sinusoidal_data_v1.mat
srate = 400;
timex = linspace(-1500,0,srate*1.5);


%%
[acf,lags,bounds] = autocorr(raw_data,'numlags',srate);
timex2 = linspace(0,1500,srate*1.5);
n = 2;
figure('position',[10 10 400 600]);
subplot(n,1,1);
plot(timex2,raw_data);
ylim([-250 200]);

subplot(n,1,2);
plot(lags, acf);
ylim([-1 1]);

