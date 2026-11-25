function [cho_output, srate, param, figureHandle] = example_CHO_v27()
%EXAMPLE_CHO_V27 Run the wavelet detector on the public hippocampal fixture.
% From the repository root: add_path; cho_output = example_CHO_v27;
root = fileparts(mfilename('fullpath'));
fixture = load(fullfile(root, 'data', 'data_for_figure8.mat'), 'raw_data');
srate = 500;
param = cho_v27_parameters();
param.return_maps = 1;
cho_output = CHO_v27(fixture.raw_data, srate, param);
figureHandle = plotExample(fixture.raw_data, srate, cho_output);
fprintf('CHO v27: %d episodes, %g Hz sampling, %g-%g Hz analysis.\n', ...
    cho_output.n_final_boxes, srate, min(cho_output.freqs), max(cho_output.freqs));
end

function fig = plotExample(signal, fs, output)
% Display the low-frequency detections clearly; analysis still spans 1-150 Hz.
times = (0:numel(signal)-1) / fs;
fig = figure('Color','w','Position',[80 80 1150 700]);
layout = tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');
ax = nexttile(layout);
plot(ax,times,signal,'Color',[0.12 0.25 0.40],'LineWidth',0.8);
ylabel(ax,'Amplitude (data units)');
title(ax,sprintf('CHO v27: fixed 8-cycle preset, %d detected episodes',output.n_final_boxes));
xlim(ax,[times(1) times(end)]); box(ax,'off');
ax = nexttile(layout);
shown = output.freqs <= 40;
imagesc(ax,times,output.freqs(shown),output.denoised_env_data(:,shown)');
set(ax,'YDir','normal');
clim(ax,[0 3*output.std_residue]);
colormap(ax,parula);
hold(ax,'on');
for k=1:numel(output.bounding_boxes)
    b=output.bounding_boxes(k);
    x1=(b.start-1)/fs; x2=(b.stop-1)/fs;
    rectangle(ax,'Position',[x1 b.minF x2-x1 b.maxF-b.minF], ...
        'EdgeColor','w','LineWidth',1.4);
    text(ax,(x1+x2)/2,b.maxF+0.6,sprintf('%.1f Hz',b.center_fp), ...
        'HorizontalAlignment','center','VerticalAlignment','bottom', ...
        'Color','w','BackgroundColor',[0.12 0.16 0.22],'Margin',1,'FontSize',9);
end
xlim(ax,[times(1) times(end)]); ylim(ax,[1 40]);
xlabel(ax,'Time (s)'); ylabel(ax,'Frequency (Hz)');
title(ax,'Analysis: 1-150 Hz; display: 1-40 Hz');
bar=colorbar(ax); bar.Label.String='Residual log_{10} wavelet magnitude';
end
