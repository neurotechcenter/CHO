function [yes_or_no] = isBRODMANN(elec_number, elec_label_list, bnumber_list)
%createColormapFromAnnotations - Creates a colormap for Annotation of the
%Surface 
% surface - Surface Data object
% returns:
% annotation_remap - remapped annotations from 1 to max
% cmap - colormap associated with annotations
% See also Model3DView, plot3DModel


% generate blabel list
for i = 1:length(bnumber_list)
    brod_labels{i} = sprintf('Brodmann area %d',bnumber_list(i));
end


elec_label = elec_label_list{elec_number};
targets= find(cellfun(@(x)any(strcmpi(x,elec_label)),brod_labels));

if isempty(targets)
    yes_or_no = 0;
else
    yes_or_no = 1;
end
    

end