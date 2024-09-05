function [ ] = ecog_topo_plot(topo_param,topo_vector)

cmapstruct = topo_param.cmapstruct;
cortex = topo_param.cortex;
tala = topo_param.tala;
viewstruct = topo_param.viewstruct;
ix = topo_param.ix;
n_ch = length(tala.electrodes);
right_hemi = topo_param.right_hemi;
topo_param.tala.activations = topo_vector;
view = topo_param.view;
alpha_level = topo_param.alpha;

%% Use activate brain package
if strcmpi('lateral',view)
    if right_hemi > 0
        viewstruct.viewvect = [-90 0];
        viewstruct.lightpos = [-200,0,0];
    else
        viewstruct.viewvect = [90 0];
        viewstruct.lightpos = [200,0,0];
    end
elseif strcmpi('superior',view)
    viewstruct.viewvect = [90 90];
    viewstruct.lightpos = [0,0,200];    
elseif strcmpi('inferior',view)
    viewstruct.viewvect = [-90 -90];
    viewstruct.lightpos = [0,0,-200];
else
    
end

viewstruct.what2view    = {'brain'};

viewstruct.material     = 'dull';
viewstruct.enablelight  = 1;
viewstruct.enableaxis   = 0;
viewstruct.lightingtype = 'gouraud';


% cmapstruct.basecol = [0.7, 0.7, 0.7];

activateBrain(cortex, topo_param.vcontribs, tala, ix, cmapstruct, viewstruct);
alpha(alpha_level);
% plotBalls(topo_param.tala.trielectrodes, 'k', 1);
% NeuralAct(cortex, [], tala, 1, cmapstruct, viewstruct);



% function [ ] = ecog_topo_plot(topo_info,right_hemi)
% 
% cmapstruct = topo_info.cmapstruct;
% cortex = topo_info.cortex;
% tala = topo_info.tala;
% vcontribs = topo_info.vcontribs;
% viewstruct = topo_info.viewstruct;
% ix = topo_info.ix;
% 
% %% Use activate brain package
% if right_hemi > 0
%     viewstruct.viewvect = [-90 0];
%     viewstruct.lightpos = [-200,0,0];
% else
%     viewstruct.viewvect = [90 0];
%     viewstruct.lightpos = [200,0,0];
% end
% 
% viewstruct.what2view = {'brain' 'activations'};
% 
% activateBrain(cortex, vcontribs, tala, ix, cmapstruct, viewstruct);
% % NeuralAct(cortex, vcontribs, tala, 1, cmapstruct, viewstruct);

