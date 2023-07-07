function [outputs] = oof_fitting_v2(psd_in_db,frequency_domain)
%FOOOF_MATLAB_V1 Summary of this function goes here
%   Detailed explanation goes here
% psd_in_db and frequency_domain should be column vector
% frequency_domain should avoid 0

%% functions for fitting
modelfun = @(b,x) b(1) - log10(b(2) + x(:,1).^b(3));


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

if length(spectrum_ignore) > 3

    mdl_coef = lsqcurvefit(modelfun,beta0,freqs_ignore,spectrum_ignore,lb,ub,options);
    initial_fit2 = modelfun(mdl_coef, X);
    
    %% results
    outputs.initial_ap_fit_coef = mdl_coef;
    outputs.init_ap_fit = initial_fit2;
    % outputs.final_fit = final_fit;
else
    outputs.initial_ap_fit_coef = nan;
    outputs.init_ap_fit = NaN(size(psd_in_db));
end

%%
% figure();
% plot(X,psd_in_db);
% hold on;
% plot(X,initial_fit2);
% fprintf('error rate = %d\n',error_rate);

end

