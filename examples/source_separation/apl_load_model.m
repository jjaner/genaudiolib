function [W , param] = apl_load_model(nmf_model, param)
% function [W , param] = apl_load_model(nmf_model, param)
%
% Loads an existing supervised NMF model (stored in HDF5 format)

svnmf.sv = nmf_model; %targetModel
hdf5info(svnmf.sv);

model_params.F0Table = hdf5read(svnmf.sv, 'F0Table');
[~, param.firstNote] = midiFreqConversion([],model_params.F0Table(1));
model_params.stepNotes = 1;
model_params.basisPerF0 = int32(hdf5read(svnmf.sv, '/parameters/nmfBasisPerF0Count'));
model_params.basisWidth = int32(hdf5read(svnmf.sv, '/parameters/nmfBasisBlockSize'));
param.F0Table = model_params.F0Table;
param.basisPerF0 = model_params.basisPerF0;
param.stepNotes = model_params.stepNotes;

W = hdf5read(svnmf.sv, 'Wtarget'); % model

param.norm_w = getNormalizationMode(svnmf.sv);
param.width = size(W,3); % frames per basis

end




%-------------------------------------------------------------
%-------------------------------------------------------------
% Additional methods
%-------------------------------------------------------------


% get the normalization mode of the CNMF model in the HDF5 format
function norm_w = getNormalizationMode(modelFile)

norm_w = 3; % Type of CNMF basis normalization.
mode = hdf5read(modelFile, '/parameters/nmfNormalizationMode');
if isnumeric (mode)
norm_w = mode;
else
if ischar(mode.Data)
switch mode.Data,
case '1-norm'
norm_w = 1;
case '2-norm'
norm_w = 2;
case '2-norm-patch',
norm_w = 3;
case 'energy-patch'
norm_w = 4;
end

end
end
end