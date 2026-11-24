cho_root = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(cho_root, 'functions')));
addpath(genpath(fullfile(cho_root, 'data')));

if exist('ft_defaults', 'file') == 2
    ft_defaults
end
if exist('eeglab', 'file') == 2
    eeglab
end
clear cho_root
