%% make table from ROI output
ExportPATH = '/Volumes/Research/CLPS_Shenhav_Lab/ScanningData/BASB/Export/';



%% FOR BAS ROIS
load(sprintf('%sBAS_S1_S2_ROIS.mat', ExportPATH))

%%

allMs = permute(allMs, [3  2  1]);

%%
dataTable = reshape(allMs, [size(allMs,1)*size(allMs,2), size(allMs,3)]);

%%

dataTable= asinh(dataTable);
%%
%[allMs(:,2,1), dataTable(145:288,1)] credibility check: passed!

ROIdataTable = array2table(dataTable);
%% 
for ncnames= 1: length(Rs)
orgname=Rs{ncnames}(16: end-4);
renamed= strrep(orgname, '-','_');
hdr{ncnames}=renamed;
end

%% 
ROIdataTable.Properties.VariableNames=hdr;

%% export table for analysis in R:
save(sprintf('%sROItable.mat', ExportPATH), 'ROIdataTable')

writetable(ROIdataTable,'~/Dropbox (Brown)/ShenhavLab/experiments/bas/Analysis/BASBWfMRI/BASB_508ROIData_053118.xls');

%% Bartra ROIS
load(sprintf('%sBARTRA_FU_ROIS.mat', ExportPATH))

%%

allMs = permute(allMs, [3  2  1]);

%%
dataTable = reshape(allMs, [size(allMs,1)*size(allMs,2), size(allMs,3)]);

%%

dataTable= asinh(dataTable);


ROIdataTable = array2table(dataTable);
%% 
for ncnames= 1: length(Rs)
orgname=Rs{ncnames}(13: end-4);
renamed= strrep(orgname, '-','_');
hdr{ncnames}=renamed;
end

%% 
ROIdataTable.Properties.VariableNames=hdr;

%% export table for analysis in R:
save(sprintf('%sBARTRAROItable.mat', ExportPATH), 'ROIdataTable')

writetable(ROIdataTable,'~/Dropbox (Brown)/ShenhavLab/experiments/bas/Analysis/BASBWfMRI/BASB_508BARTRAROIData_053118.xls');

%% Export BARTRA resp locked


load(sprintf('%sBARTRA_FULL_ROIS_resp.mat', ExportPATH))

%%

allMs = permute(allMs, [3  2  1]);

%%
dataTable = reshape(allMs, [size(allMs,1)*size(allMs,2), size(allMs,3)]);

%%

dataTable= asinh(dataTable);
%%

ROIdataTable = array2table(dataTable);
%% 
for ncnames= 1: length(Rs)
orgname=Rs{ncnames}(13: end-4);
renamed= strrep(orgname, '-','_');
hdr{ncnames}=renamed;
end

%% 
ROIdataTable.Properties.VariableNames=hdr;

%% export table for analysis in R:
save(sprintf('%sBARTRAROIRESPtable.mat', ExportPATH), 'ROIdataTable')

writetable(ROIdataTable,'~/Dropbox (Brown)/ShenhavLab/experiments/bas/Analysis/BASBWfMRI/BASB_508BARTRARESPROIData.xls');


%% The other Bartras
load(sprintf('%sBARTRA_*n_ROIS.mat', ExportPATH))

%%
allMs = permute(allMs, [3  2  1]);

%%
dataTable = reshape(allMs, [size(allMs,1)*size(allMs,2), size(allMs,3)]);

%%

dataTable= asinh(dataTable);
%%


ROIdataTable = array2table(dataTable);
%% 
for ncnames= 1: length(Rs)
orgname=Rs{ncnames}(1: 11);
renamed= strrep(orgname, '-','_');
hdr{ncnames}=renamed;
end

%% 
ROIdataTable.Properties.VariableNames=hdr;

%% export table for analysis in R:
save(sprintf('%sBARTRAROIvmPFCvstrtable.mat', ExportPATH), 'ROIdataTable')

writetable(ROIdataTable,'~/Dropbox (Brown)/ShenhavLab/experiments/bas/Analysis/BASBWfMRI/BASB_508BARTRAROIsubRegData_112118.xls');

%% export the sphere stuff

load(sprintf('%sleft_righ_ROIS.mat', ExportPATH))

%%
allMs = permute(allMs, [3  2  1]);

%%
dataTable = reshape(allMs, [size(allMs,1)*size(allMs,2), size(allMs,3)]);

%%

dataTable= asinh(dataTable);
%%

ROIdataTable = array2table(dataTable);
%% 
for ncnames= 1: length(Rs)
orgname=Rs{ncnames}(24: end-4);
renamed= strrep(orgname, '-','_');
hdr{ncnames}=renamed;
end

%% 
ROIdataTable.Properties.VariableNames=hdr;

%% export table for analysis in R:
save(sprintf('%sBASB_508_left_righ_Str_ROIS.mat', ExportPATH), 'ROIdataTable')

writetable(ROIdataTable,'~/Dropbox (Brown)/ShenhavLab/experiments/bas/Analysis/BASBWfMRI/BASB_508_left_righ_Str_ROIS.xls');


%% OVr effect ROI OVr_L_R_ROIS.mat

load(sprintf('%sOVr_L_R_ROIS.mat', ExportPATH))

%%
allMs = permute(allMs, [3  2  1]);

%%
dataTable = reshape(allMs, [size(allMs,1)*size(allMs,2), size(allMs,3)]);

%%

dataTable= asinh(dataTable);
%%

ROIdataTable = array2table(dataTable);
%% 
for ncnames= 1: length(Rs)
orgname=Rs{ncnames}(1: end-4);
renamed= strrep(orgname, '-','_');
hdr{ncnames}=renamed;
end

%% 
ROIdataTable.Properties.VariableNames=hdr;

%% export table for analysis in R:
save(sprintf('%sBASB_508_OVr_ROI.mat', ExportPATH), 'ROIdataTable')

writetable(ROIdataTable,'~/Dropbox (Brown)/ShenhavLab/experiments/bas/Analysis/BASBWfMRI/BASB_508_OVr_ROI.xls');


%% pgACC/mOFC  ROI Anat_ProbAt_ROIS.mat

load(sprintf('%sAnat_ProbAt_ROIS.mat', ExportPATH))

%%
allMs = permute(allMs, [3  2  1]);

%%
dataTable = reshape(allMs, [size(allMs,1)*size(allMs,2), size(allMs,3)]);

%%

dataTable= asinh(dataTable);
%%
%[allMs(:,2,1), dataTable(145:288,1)] credibility check: passed!

ROIdataTable = array2table(dataTable);
%% 
for ncnames= 1: length(Rs)
orgname=Rs{ncnames}(1: end-4);
renamed= strrep(orgname, '-','_');
hdr{ncnames}=renamed;
end

%% 
ROIdataTable.Properties.VariableNames=hdr;

%% export table for analysis in R:
save(sprintf('%sVMPFC_OFC_ROI.mat', ExportPATH), 'ROIdataTable')

writetable(ROIdataTable,'~/Dropbox (Brown)/ShenhavLab/experiments/bas/Analysis/BASBWfMRI/VMPFC_OFC_ROI.xls');
