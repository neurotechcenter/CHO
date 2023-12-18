clear;
clc;


load data_for_figure8.mat
srate = 500;


%%
param.plot = 1;
param.minimum_cycles = 2;
param.ovlp_threshold = 0.50;
param.frequency_vector = 1:40;

[cho_output] = CHO_v22(raw_data, srate, param);


