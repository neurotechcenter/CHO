function tests = testWaveletV27
tests = functiontests(localfunctions);
end

function setup(testCase)
root=fileparts(fileparts(mfilename('fullpath')));
previous=struct('folder',pwd,'path',path,'warnings',warning,'random',rng, ...
    'visibility',get(groot,'DefaultFigureVisible'), ...
    'figures',findall(groot,'Type','figure'));
testCase.addTeardown(@() restoreSession(previous));
cd(root); add_path
set(groot,'DefaultFigureVisible','off');
testCase.TestData.root=root;
end

function testPublicResearchPreset(testCase)
previousWarnings=warning;
[output,fs,param,fig]=example_CHO_v27();
testCase.verifyEqual(warning,previousWarnings);
testCase.verifyEqual(which('CHO_v27'), ...
    fullfile(testCase.TestData.root,'functions','v27','CHO_v27.m'));
testCase.verifyEqual(fs,500);
testCase.verifyEqual(param.wavelet_cycles,8);
testCase.verifyEqual(param.merge_mode,'fast_tf_overlap');
testCase.verifyEqual(param.minimum_cycles,3);
testCase.verifyEqual(param.osc_threshold_std,2);
testCase.verifyEqual(output.freqs,1:150);
testCase.verifyEqual(output.param.wavelet_cycles_by_freq,8*ones(1,150));
testCase.verifyEqual(output.status,"ok");
testCase.verifyEqual(output.n_final_boxes,7);
testCase.verifySize(output.denoised_env_data,[5000 150]);
testCase.verifyTrue(all(isfinite(output.ap_fit)));
validateBoxes(testCase,output.bounding_boxes,5000,[1 150]);
testCase.assertNumElements(findall(fig,'Type','axes'),2);
testCase.verifyNumElements(findall(fig,'Type','rectangle'),7);
destination=fullfile(testCase.TestData.root,'artifacts');
if ~isfolder(destination),mkdir(destination);end
exportgraphics(fig,fullfile(destination,'cho-v27-hippocampal.png'),'Resolution',150);
save(fullfile(destination,'cho-v27-output.mat'),'output','fs','param');
fprintf('V27 fixed-8 preset: %d episodes, %g Hz sampling, 1-150 Hz.\n',output.n_final_boxes,fs);
end

function testDetectorDefaultsRemainDistinct(testCase)
fixture=load(fullfile(testCase.TestData.root,'data','data_for_figure8.mat'));
previousWarnings=warning;
output=CHO_v27(fixture.raw_data,500);
testCase.verifyEqual(warning,previousWarnings);
testCase.verifyEqual(output.param.wavelet_cycles,'linear_3_10');
testCase.verifyEqual(output.param.merge_mode,'correlation');
testCase.verifyEqual(output.freqs,1:40);
testCase.verifyEqual(output.param.wavelet_cycles_by_freq,linspace(3,10,40));
testCase.verifyEqual(output.n_final_boxes,5);
validateBoxes(testCase,output.bounding_boxes,5000,[1 40]);
end

function testFlatInputAndNyquistGuard(testCase)
param=cho_v27_parameters();
empty=CHO_v27(zeros(1000,1),100,param);
testCase.verifyEqual(empty.status,"invalid_or_flat_signal");
testCase.verifyEmpty(empty.bounding_boxes);
fixture=load(fullfile(testCase.TestData.root,'data','data_for_figure8.mat'));
output=CHO_v27(fixture.raw_data,100,param);
testCase.verifyEqual(output.freqs,1:49);
testCase.verifyEqual(output.n_final_boxes,numel(output.bounding_boxes));
validateBoxes(testCase,output.bounding_boxes,5000,[1 49]);
end

function validateBoxes(testCase,boxes,nSamples,limits)
testCase.assertNotEmpty(boxes);
testCase.verifyTrue(all(isfinite([boxes.center_fp])));
testCase.verifyTrue(all(isfinite([boxes.peak_val])));
testCase.verifyTrue(all([boxes.start]>=1 & [boxes.stop]<=nSamples));
testCase.verifyTrue(all([boxes.start]<[boxes.stop]));
testCase.verifyTrue(all([boxes.minF]>=limits(1) & [boxes.maxF]<=limits(2)));
end

function restoreSession(previous)
delete(setdiff(findall(groot,'Type','figure'),previous.figures));
set(groot,'DefaultFigureVisible',previous.visibility);
warning(previous.warnings); rng(previous.random);
path(previous.path); cd(previous.folder);
end
