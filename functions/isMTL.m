function [yes_or_no] = isMTL(elec_number, elec_label_list)
%createColormapFromAnnotations - Creates a colormap for Annotation of the
%Surface 
% surface - Surface Data object
% returns:
% annotation_remap - remapped annotations from 1 to max
% cmap - colormap associated with annotations
% See also Model3DView, plot3DModel




limbic_labels = {'Corticoamygdaloid-transitio','Anterior-amygdaloid-area-AAA',...
    'hippocampal_fissure',...
    'CA3-head','CA3-body',...
    'CA1-head','CA1-body','CA4-body', 'CA4-head',...
    'GC-ML-DG-head', 'GC-ML-DG-body', 'GC-ML-DG-body',...
    'molecular_layer_HP-body','molecular_layer_HP-head',...    
    'Accessory-Basal-nucleus','Basal-nucleus','Lateral-nucleus','Right-Inf-Lat-Vent',...
    'fimbria','HATA',...
    'ctx-rh-entorhinal','presubiculum-head','subiculum-head','subiculum-body',...
    'ctx-rh-parahippocampal','ctx-lh-parahippocampal',...
    };

elec_label = elec_label_list{elec_number};
targets= find(cellfun(@(x)any(strcmpi(x,elec_label)),limbic_labels));

if isempty(targets)
    yes_or_no = 0;
else
    yes_or_no = 1;
end
    

end