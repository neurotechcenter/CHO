clear;
clc;

osc_amp_list = 5:20;
n_iteration = 10;

for snr_iter = 1:length(osc_amp_list)
    
    for iter = 1:n_iteration
        fprintf('Iteration #%d\n',iter);
        sim_param.srate = 500;
        sim_param.osc_amp = osc_amp_list(snr_iter);
        sim_param.freq_range = [2 30];
        sim_param.noise_amp = 10;
        sim_param.n_iteration = 10;
        sim_param.asym_ratio = 0.1; % 0< beta_scaler < 1
        [ACC{snr_iter,iter} SENS{snr_iter,iter} SPEC{snr_iter,iter} SNRs{snr_iter,iter}] = harmonic_osc_simulation_TACC_v5(sim_param);
    end    
end
% 
%%
save('./data/HS_TACC_simulation_MSNR_v5.mat', 'ACC', 'SENS','SPEC','SNRs');
%%
load HS_TACC_simulation_MSNR_v5.mat
%%
SNR_data = [];
n_methods = 2;
for snr_iter = 1:length(osc_amp_list)
    for iter = 1:n_iteration
        SNR_data = [SNR_data SNRs{snr_iter,iter}]; 
    end
end


for imethod = 1:n_methods
    ACC_data{imethod} = []; 
    SPEC_data{imethod} = [];
    SENS_data{imethod} = [];
    
    for snr_iter = 1:length(osc_amp_list)
        for iter = 1:n_iteration
            ACC_data{imethod} = [ACC_data{imethod} ACC{snr_iter,iter}(imethod,:)];
            SPEC_data{imethod} = [SPEC_data{imethod} SPEC{snr_iter,iter}(imethod,:)];
            SENS_data{imethod} = [SENS_data{imethod} SENS{snr_iter,iter}(imethod,:)];
        end
    end
end

%%
% figure();
% plot(SNR_data',SPEC_data{1}', '.');

%%
y_axis = 0:0.1:1;
snr_axis = -26:4:0;
% figure();
% hist3([SNR_data; SPEC_data{1}]','Ctrs',{snr_axis y_axis},'CDataMode','auto');

% prepare bins    
Method_SPEC = [];
Method_SENS = [];
Method_ACC = [];
for imethod = 1:n_methods
    for j = 1:length(snr_axis)
        SNR_bins_SPEC{j} = [];
        SNR_bins_SENS{j} = [];
        SNR_bins_ACC{j} = [];
    end
    for i = 1:length(SNR_data)        
        if SNR_data(i) <= snr_axis(1)
            SNR_bins_SPEC{1} = [SNR_bins_SPEC{j}; SPEC_data{imethod}(i)];
            SNR_bins_SENS{1} = [SNR_bins_SENS{j}; SENS_data{imethod}(i)];
            SNR_bins_ACC{1} = [SNR_bins_ACC{j}; ACC_data{imethod}(i)];
        elseif SNR_data(i) > snr_axis(end)
            SNR_bins_SPEC{length(snr_axis)} = [SNR_bins_SPEC{j}; SPEC_data{imethod}(i)];
            SNR_bins_SENS{length(snr_axis)} = [SNR_bins_SENS{j}; SENS_data{imethod}(i)];
            SNR_bins_ACC{length(snr_axis)} = [SNR_bins_ACC{j}; ACC_data{imethod}(i)];
        else
            for j = 1:length(snr_axis)-1
                if SNR_data(i) > snr_axis(j) && SNR_data(i) <= snr_axis(j+1)
                    SNR_bins_SPEC{j} = [SNR_bins_SPEC{j}; SPEC_data{imethod}(i)];
                    SNR_bins_SENS{j} = [SNR_bins_SENS{j}; SENS_data{imethod}(i)];
                    SNR_bins_ACC{j} = [SNR_bins_ACC{j}; ACC_data{imethod}(i)];
                end
            end
        end        
    end    
    clear bin_length_SPEC bin_length_SENS bin_length_ACC
    for j = 1:length(snr_axis)
        bin_length_SPEC(j) = length(SNR_bins_SPEC{j});
        bin_length_SENS(j) = length(SNR_bins_SPEC{j});
        bin_length_ACC(j) = length(SNR_bins_ACC{j});
    end
    max_bin_length_SPEC = max(bin_length_SPEC);  
    max_bin_length_SENS = max(bin_length_SENS);
    max_bin_length_ACC = max(bin_length_ACC);
    Method_SPEC{imethod} = nan(max_bin_length_SPEC,length(snr_axis));
    Method_SENS{imethod} = nan(max_bin_length_SENS,length(snr_axis));
    Method_ACC{imethod} = nan(max_bin_length_ACC,length(snr_axis));
    for j = 1:length(snr_axis)
        for i = 1:length(SNR_bins_SPEC{j})
            Method_SPEC{imethod}(i,j) = SNR_bins_SPEC{j}(i);
        end
        for i = 1:length(SNR_bins_SENS{j})
            Method_SENS{imethod}(i,j) = SNR_bins_SENS{j}(i);
        end
        for i = 1:length(SNR_bins_ACC{j})
            Method_ACC{imethod}(i,j) = SNR_bins_ACC{j}(i);
        end
    end
end

%%
clear aline
title_methods = {'CHO', 'OEvent', 'Fooof', 'SPRiNT'};
n_plot = 3;
nrSamples = 2;
cMap = lines(nrSamples);
figure('position',[10 10 600 900]);
subplot(n_plot,1,1);
for method_iter = 1 : nrSamples
    aline(method_iter) = stdshade(Method_ACC{method_iter},0.1,cMap(method_iter,:));
%     set(gca,'xtick',1:length(snr_axis));
    set(gca,'xticklabel',snr_axis+2);
    hold on;
end
legend(aline,title_methods,'Location','northeastoutside');
xlabel('SNR (dB)');
ylabel('ACCURACY');
ylim([0 1]);

subplot(n_plot,1,2);
for method_iter = 1 : nrSamples
    aline(method_iter) = stdshade(Method_SENS{method_iter},0.1,cMap(method_iter,:));
    set(gca,'xticklabel',snr_axis+2);
    hold on;
end
legend(aline,title_methods,'Location','northeastoutside');
xlabel('SNR (dB)');
ylabel('SENSITIVITY');
ylim([0 1]);


subplot(n_plot,1,3);
for method_iter = 1 : nrSamples
    aline(method_iter) = stdshade(Method_SPEC{method_iter},0.1,cMap(method_iter,:));
    set(gca,'xticklabel',snr_axis+2);
    hold on;
end
legend(aline,title_methods,'Location','northeastoutside'); 
xlabel('SNR (dB)');
ylabel('SPECIFICITY');
ylim([0 1]);
