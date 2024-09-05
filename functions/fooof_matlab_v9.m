function [gaussians,outputs] = fooof_matlab_v9(psd_in_db,frequency_domain)
%FOOOF_MATLAB_V1 Summary of this function goes here
%   Detailed explanation goes here
% psd_in_db and frequency_domain should be column vector
% frequency_domain should avoid 0

%% functions for fitting
modelfun = @(b,x) b(1) - log10(b(2) + x(:,1).^b(3));
% modelfun_l = @(b,x) b(1) * x(:,1) + b(2);
% modelfun = @(b,x) b(1) - b(2)*log10(x(:,1).^b(3));
modelfun_g = @(b,x) b(1) * exp(-(x(:, 1) - b(2)).^2/b(3)) + b(4);
% modelfun_mg = @(b,x) b(1) * exp(-(x(:, 1) - b(2)).^2/b(3)) + b(4) * exp(-(x(:, 1) - b(5)).^2/b(6)) + b(7) * exp(-(x(:, 1) - b(8)).^2/b(9));  
modelfun_mg = @(b,x) b(1) * exp(-(x(:, 1) - b(2)).^2/b(3)) + b(4) * exp(-(x(:, 1) - b(5)).^2/b(6)) + b(7) * exp(-(x(:, 1) - b(8)).^2/b(9)) + b(10) * exp(-(x(:, 1) - b(11)).^2/b(12));  
const_for_std = (2*(2*log(2))^0.5);

%% settings
gaussians = [];
outputs = [];
options = optimoptions('lsqcurvefit','Display','off');
MAX_GAUSSIANS = 4;
X = frequency_domain;
Y = psd_in_db;
ap_percentile_thresh = 10;

%% Fit aperiodic signal - crucial and not trivial

beta0 = [min(Y) 0 1];
lb = [min(Y)-2 0 0];
ub = [max(Y)+2 0 1000 ];

mdl_coef_init = lsqcurvefit(modelfun,beta0,X,Y,lb,ub,options);
initial_fit = modelfun(mdl_coef_init, X);

flatspec = Y - initial_fit;

flatspec(flatspec<0) = 0;
perc_thresh = prctile(flatspec,ap_percentile_thresh);
perc_mask = flatspec <= perc_thresh;
freqs_ignore = X(perc_mask);
spectrum_ignore = Y(perc_mask);

mdl_coef = lsqcurvefit(modelfun,beta0,freqs_ignore,spectrum_ignore,lb,ub,options);

% robust_fit = modelfun(mdl_coef, X);
% beta0 = [-1 0 ];
% lb = [-1000 -1000];
% ub = [1000 1000];
% 
% mdl_coef_l = lsqcurvefit(modelfun_l,beta0,log10(X([1 length(X)])),Y([1 length(X)]),lb,ub,options);
% log_ap = modelfun_l(mdl_coef_l,log10(X));
% 
% beta0 = [min(Y) 0 1];
% lb = [min(Y)-2 0 0];
% ub = [max(Y)+2 0 1000 ];
% 
% mdl_coef = lsqcurvefit(modelfun,beta0,X,log_ap,lb,ub,options);

%% Fit and remove Gaussains
n_gaussians = 0;
offsets(1) = 0;
while(1)
    if n_gaussians == 0
        ap_removed_signal = Y - modelfun(mdl_coef, X);
        error_rate_init_ap = std(ap_removed_signal)^2/std(psd_in_db)^2*100;
    else
        ap_removed_signal = Y - modelfun_g(mdl_coef_g{n_gaussians}, X);
        error_rate_init_ap = std(ap_removed_signal)^2/std(psd_in_db)^2*100;
    end
    
    threshold(n_gaussians+1) = std(abs(ap_removed_signal));
    [pks,locs]= findpeaks(ap_removed_signal,'MinPeakHeight',threshold(n_gaussians+1));
    
    if(length(locs) == 0 | n_gaussians == MAX_GAUSSIANS | error_rate_init_ap < .1)
        break;
    end
    
    [maxpks maxloc] = max(pks);
    [sorted_pks sorted_idx] = sort(pks,'descend');
    maxpks = sorted_pks(1);
    maxloc = sorted_idx(1);
    
    half_max_power = maxpks/2;
    [minval minidxs] = min(abs(ap_removed_signal - half_max_power));
    approx_std = 2*(X(locs(maxloc)) - X(minidxs))/const_for_std;
    
    %% select the best peak
    n_pks = length(pks);
    candidates = [];
    for ipks = 1:n_pks
        current_m = X(locs(maxloc));
        half_max_power = pks(ipks)/2;
        right_flank = locs(maxloc);
        left_flank = locs(maxloc);
        for fp = locs(maxloc):length(X)
            if ap_removed_signal(fp) < half_max_power
                right_flank = fp;
                break;
            end
        end
        for fp = locs(maxloc):-1:1
            if ap_removed_signal(fp) < half_max_power
                left_flank = fp;
                break;
            end
        end
        min_val = min([X(left_flank) X(right_flank)]);
        current_std = 2*abs(current_m-min_val)/const_for_std;
        
        if n_gaussians == 0
            cond2 = current_m - current_std < X(1) && X(1) < current_m + current_std;
            cond3 = current_m - current_std < X(end) && X(end) < current_m + current_std;
            if cond2 | cond3
                candidates(ipks) = 0;
            else
                candidates(ipks) = 1;
            end
        else
            previous_m = mdl_coef_g{n_gaussians}(2);
            previous_std = (mdl_coef_g{n_gaussians}(3)/2)^0.5;
            
            cond1 = previous_m - 0.75*previous_std < current_m && current_m < previous_m + 0.75*previous_std;
            cond2 = current_m - current_std < X(1) && X(1) < current_m + current_std;
            cond3 = current_m - current_std < X(end) && X(end) < current_m + current_std;
            if cond1 | cond2 | cond3
                candidates(ipks) = 0;
            else
                candidates(ipks) = 1;
            end
        end
    end
    % select max power peak
    valid_pk_idxs = find(candidates>0);
    final_pks_list = pks(valid_pk_idxs);
    final_pk_locs_list = locs(valid_pk_idxs);
    [final_maxpks final_maxloc] = max(final_pks_list);
    
    final_pk = X(final_pk_locs_list(final_maxloc));
    final_pk_power = final_maxpks;
    
    if length(final_pk) == 0
        break;
    end

%%
%     final_pk = X(locs(maxloc));
    
    %%
    Y = ap_removed_signal;
    beta0 = [final_pk_power final_pk 1 0];
    lb = [0 final_pk-0.1 0 min(ap_removed_signal)-1];
    ub = [final_pk_power+5 final_pk+0.1 50 max(ap_removed_signal)+1];
%     lb = [0 locs(maxloc)-0.1 0 -10];
%     ub = [100 locs(maxloc)+0.1 10000 10];
    
    mdl_coef_g{n_gaussians+1} = lsqcurvefit(modelfun_g,beta0,X,Y,lb,ub,options);
%     yFitted_g(n_gaussians+1,:) = modelfun_g(mdl_coef_g{n_gaussians+1}, D');
    
    offsets(n_gaussians+1) = mdl_coef_g{n_gaussians+1}(end);
    n_gaussians = n_gaussians + 1;
    %%
%     figure();
%     plot(X,ap_removed_signal);
%     hold on;
%     vline(X(locs(maxloc)));
%     vline(X(minidxs));
end

%% Multi-gaussian fit using iteration parameters
Y = psd_in_db;
ap_removed_signal_org = Y - modelfun(mdl_coef, X) - sum(offsets);
% ap_removed_signal_org = Y - modelfun(mdl_coef, X);
if n_gaussians == 0
    mdl_coef_mg = zeros(1,MAX_GAUSSIANS*3);
else
    beta0 = zeros(1,MAX_GAUSSIANS*3);
    lb = zeros(1,MAX_GAUSSIANS*3);
    ub = zeros(1,MAX_GAUSSIANS*3);
    for iter = 1:n_gaussians
        beta0((iter-1)*3 + 1) = mdl_coef_g{iter}(1);
        beta0((iter-1)*3 + 2) = mdl_coef_g{iter}(2);
        beta0((iter-1)*3 + 3) = mdl_coef_g{iter}(3);
        
        lb((iter-1)*3 + 1) = 0;
        lb((iter-1)*3 + 2) = mdl_coef_g{iter}(2)-.1;
        lb((iter-1)*3 + 3) = -100;
        
        ub((iter-1)*3 + 1) = 100;
        ub((iter-1)*3 + 2) = mdl_coef_g{iter}(2)+.1;
        ub((iter-1)*3 + 3) = 100;
    end
    
    mdl_coef_mg = lsqcurvefit(modelfun_mg,beta0,X,ap_removed_signal_org,lb,ub,options);
end

%% finalize number of gaussians
% correction_g = zeros(1,MAX_GAUSSIANS);
% for giter = 1:MAX_GAUSSIANS
%     gpower(giter) = mdl_coef_mg((giter-1)*3 + 1);
%     if gpower(giter) > 0
%         correction_g(giter) = 1;
% %         fprintf('corrected!\n');
%     end
% end
% n_gaussians = sum(correction_g);
%% remove gaussians from original PSD and Re-fit aperiodic
pure_ap_signal = Y - modelfun_mg(mdl_coef_mg, X);

Y = pure_ap_signal;
beta0 = mdl_coef;
lb = [-100 0 0 ];
ub = [100 0 100];
mdl_coef_new_ap = lsqcurvefit(modelfun,beta0,X,Y,lb,ub,options);

%% combine and calculate final residuals
ap_fit = modelfun(mdl_coef_new_ap, X);
mg_fit = modelfun_mg(mdl_coef_mg, X);

final_fit = ap_fit+mg_fit;
residuals = psd_in_db - final_fit;
error_rate = std(residuals)^2/std(psd_in_db)^2*100;

%% results
if n_gaussians == 0
    outputs.initial_ap_fit_coef = mdl_coef;
    outputs.indiv_gaussian_coefs{1} = [0 0 0];
    outputs.multi_gaussian_coef = zeros(1,MAX_GAUSSIANS*3);
    outputs.final_ap_fit_coef = mdl_coef_new_ap;
    outputs.n_gaussians = 0;
    outputs.error_rate = error_rate;
    outputs.fooofed_ap_fit = ap_fit;
    outputs.final_fit = final_fit;
else
    for iter = 1:n_gaussians
        gaussians{iter}.freq = mdl_coef_mg((iter-1)*3+2);
        gaussians{iter}.power = mdl_coef_mg((iter-1)*3+1);
        gaussians{iter}.bandwidth = ((mdl_coef_mg((iter-1)*3+3)/2).^0.5)*2;
    end

    outputs.initial_ap_fit_coef = mdl_coef;
    outputs.indiv_gaussian_coefs = mdl_coef_g;
    outputs.multi_gaussian_coef = mdl_coef_mg;
    outputs.final_ap_fit_coef = mdl_coef_new_ap;
    outputs.n_gaussians = n_gaussians;
    outputs.error_rate = error_rate;
    outputs.fooofed_ap_fit = ap_fit;
    outputs.final_fit = final_fit;
end

%%
% figure();
% plot(X,psd_in_db);
% hold on;
% plot(X,ap_fit+mg_fit);
% plot(X,ap_fit,'--');
% fprintf('error rate = %d\n',error_rate);

end

