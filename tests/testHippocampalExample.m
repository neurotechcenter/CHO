function tests = testHippocampalExample
% Regression check for the public README example and its bundled fixture.
tests = functiontests(localfunctions);
end

function testReadmeExample(testCase)
root = fileparts(fileparts(mfilename('fullpath')));
originalFolder = pwd;
originalPath = path;
originalVisibility = get(groot, 'DefaultFigureVisible');
figuresBefore = findall(groot, 'Type', 'figure');
cleanup = onCleanup(@() restoreEnvironment(originalFolder, originalPath, ...
    originalVisibility, figuresBefore)); %#ok<NASGU>
cd(root);
set(groot, 'DefaultFigureVisible', 'off');

requiredFunctions = {'autocorr', 'hilbert', 'bwlabel', ...
    'optimoptions', 'lsqcurvefit'};
for k = 1:numel(requiredFunctions)
    testCase.assertNotEmpty(which(requiredFunctions{k}), ...
        ['Required toolbox function is unavailable: ' requiredFunctions{k}]);
end

% The script's clear command is isolated inside this helper function.
[output, signal, samplingRate, parameters] = executeReadmeCommands();
testCase.verifyEqual(which('CHO_v22'), fullfile(root, 'functions', 'CHO_v22.m'));
fixturePath = fullfile(root, 'data', 'data_for_figure8.mat');
testCase.verifyEqual(which('data_for_figure8.mat'), fixturePath);
fixture = load(fixturePath, 'raw_data');
testCase.verifyEqual(signal, fixture.raw_data);
testCase.verifySize(signal, [5000 1]);
testCase.verifyEqual(samplingRate, 500);
testCase.verifyEqual(parameters.frequency_axis, 1:40);
testCase.assertClass(output, 'struct');
testCase.assertSize(output, [1 1]);
testCase.assertTrue(all(isfield(output, {'bounding_boxes', ...
    'initial_ap_fit_coef', 'init_ap_fit', 'signal_std'})));
testCase.verifySize(output.initial_ap_fit_coef, [1 3]);
testCase.verifySize(output.init_ap_fit, [40 1]);
testCase.verifyTrue(all(isfinite(output.initial_ap_fit_coef(:))));
testCase.verifyTrue(all(isfinite(output.init_ap_fit(:))));
testCase.verifyEqual(output.signal_std, std(signal), 'RelTol', 1e-12);

boxes = output.bounding_boxes;
testCase.assertSize(boxes, [1 6]);
testCase.assertTrue(all(isfield(boxes, {'center_fp', 'center_tp', ...
    'peak_val', 'minF', 'maxF', 'start', 'stop'})));
testCase.verifyEqual(sort([boxes.center_fp]), [7 7 8 8 9 10]);
testCase.verifyTrue(all([boxes.start] >= 1));
testCase.verifyTrue(all([boxes.stop] <= numel(signal)));
testCase.verifyTrue(all([boxes.start] < [boxes.stop]));
testCase.verifyTrue(all(isfinite([boxes.peak_val])));

newFigures = setdiff(findall(groot, 'Type', 'figure'), figuresBefore);
testCase.assertNumElements(newFigures, 1);
testCase.verifyNumElements(findall(newFigures, 'Type', 'axes'), 4);
artifactFolder = fullfile(root, 'artifacts');
if ~isfolder(artifactFolder)
    mkdir(artifactFolder);
end
pngPath = fullfile(artifactFolder, 'hippocampal-example.png');
exportgraphics(newFigures, pngPath, 'Resolution', 150);
pngInfo = dir(pngPath);
testCase.assertNumElements(pngInfo, 1);
testCase.verifyGreaterThan(pngInfo.bytes, 0);
save(fullfile(artifactFolder, 'hippocampal-output.mat'), ...
    'output', 'samplingRate', 'parameters');
fprintf('CHO smoke test: %d episodes, %d figure, 500 Hz, 1:40 Hz.\n', ...
    numel(boxes), numel(newFigures));
end

function [output, signal, samplingRate, parameters] = executeReadmeCommands
add_path
figure8BCDE
output = cho_output;
signal = raw_data;
samplingRate = srate;
parameters = param;
end

function restoreEnvironment(folder, searchPath, visibility, previousFigures)
newFigures = setdiff(findall(groot, 'Type', 'figure'), previousFigures);
delete(newFigures);
set(groot, 'DefaultFigureVisible', visibility);
path(searchPath);
cd(folder);
end
