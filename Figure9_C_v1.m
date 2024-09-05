clear;
clc;
close all;
%%
hist_data = [];
n_subjects = 6;
srate = 400;

    
for isub = 1:6

    load(['ze' num2str(isub) '_cho_method_ecog_rest_v1.mat']);

    load([ 'ze' num2str(isub) '_brodmann.mat']);

    
    %     ch_pk_freq_prob = [];
    start_time_ch = [];
    start_time_tp = [];
        
    n_trials = size(method_outputs,2);
    n_ch = size(method_outputs,1);
    
    for ch = 1:n_ch
        
        if isBRODMANN(ch,SecondaryLabel,[44 45])
            for tr = 1:n_trials
                n_boxes = length(method_outputs(ch,tr).bounding_boxes);
                
                
                for ibox = 1:n_boxes
                    center_fp = method_outputs(ch,tr).bounding_boxes(ibox).center_fp;
                    start_time = method_outputs(ch,tr).bounding_boxes(ibox).start;
                    stop_time = method_outputs(ch,tr).bounding_boxes(ibox).stop;
                    power = method_outputs(ch,tr).bounding_boxes(ibox).peak_val;
                    
                    hist_data = [hist_data; center_fp (stop_time-start_time)/srate*1000];
                end
                
            end
%         else
            
        end
    end
    
    
    
end



%%

y_axis = 1:40;
x_axis = 0:50:5000;

% figure();
% hist(hist_data(:,1),y_axis);
% 
% %%
% figure();
% hist(hist_data(:,2),x_axis);



% h = histogram2(hist_data(:,2),hist_data(:,1),'DisplayStyle','tile','ShowEmptyBins','off');
hh3 = hist3(hist_data,'Ctrs',{y_axis x_axis},'CDataMode','auto');

figure();
subplot(3,3, [2 3 5 6]);
imagesc(x_axis,y_axis,hh3);
set(gca,'YDir','normal');
xlim([0 3000]);

xlabel('Duration (ms)');
ylabel('Frequency (Hz)');


subplot(3,3, [1 4]);
hdata=hist(hist_data(:,1),y_axis);
hist(hist_data(:,1),y_axis);
[mv maxi_hist] = max(hdata);
vline(maxi_hist,'r-');
h = findobj(gca,'Type','patch');
h.FaceColor = [0 0.5 0.5];
xlim([1 40]);
camroll(90);
set(gca,'XDir','reverse');
set(gca,'visible','off');
text(maxi_hist+2,mv,sprintf('%d Hz',maxi_hist));

subplot(3,3, [8 9]);
hdata = hist(hist_data(:,2),x_axis);
hist(hist_data(:,2),x_axis);
[mv maxi_hist] = max(hdata);
vline(x_axis(maxi_hist),'r-');
h = findobj(gca,'Type','patch');
h.FaceColor = [0 0.5 0.5];
xlim([0 3000]);
set(gca,'visible','off');
text(x_axis(maxi_hist)+100,mv,sprintf('%d ms',x_axis(maxi_hist)));



%%
frequency = 1:40;
two_cycle_duration_ms = 2./frequency *1000;


cycle_data = [two_cycle_duration_ms,frequency];

figure();
subplot(3,3, [2 3 5 6]);
plot(two_cycle_duration_ms,frequency);
xlim([0 3000]);

