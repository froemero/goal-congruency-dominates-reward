%+---------------------------------------
%|
%| Robert C. Welsh
%| 2002.10.20
%| Copyright 2002,2003,2004
%|
%|
%| Version 1.00
%|
%|
%| University of Michigan
%| Department of Radiology
%|
%| A little toolbox to apply a mask to
%| a t-map or beta-map or con-map and 
%| report the mean and variance.
%| 
%|
%| To use this you must supply a mask image
%| that is defined with the exact same 
%| dimensions as the files to extract values from.
%| 
%| You can do up to "N" group extractions.
%| For each extraction you supply a mask image
%| and a list of images to extract values from. 
%|
%| These can be Beta, T, or Con images. Actually
%| they may be any form of image, it is up to you
%| to interpret.
%|
%|
%| For this code to appear in the toolbox area it must
%| live in a subdirectory called "ExtractVals" in the 
%| spm99/toolbox directory.
%| 
%| I have no idea if it works with spm2, but I would hazard
%| that it does.
%| 
%| Please acknowledge in publications, as that will urge
%| me to write more toolboxes and release them to the public.
%|
%| Kindly,
%|
%| Robert C. Welsh
%|
%+---------------------------------------

function [AllmeanVal,AllvariVal,AllnVoxels,spmMaskFiles] = AS_ExtractVals_auto(spmMask,spmVals,printText,contBasedMasks,contBasedNumVox)
%% SHOULD either input equal ROIs/file-sets, 1 ROI & multiple file-sets, or 1 file-set & multiple ROIs


if ~exist('printText','var')
    printText = 1;
end
if ~exist('contBasedMasks','var')
    for rixi = 1:length(spmMask)
        contBasedMasks{rixi} = [];
    end
end
if ~exist('contBasedNumVox','var')
    contBasedNumVox = zeros(1,length(contBasedMasks));
end

% here's probably the problem...
if (length(spmMask)~=length(spmVals)) % so spmVals is the number of subjects while spmMasks is the nb of masks
    if (length(spmMask)>length(spmVals)) % When I have more masks than subjects I loop through masks
        maskOpt = 2;
        nExtractions = length(spmMask);
    else
        maskOpt = 1;
        nExtractions = length(spmVals); % else I loop through subjects, which I'm doing now
    end
else
    nExtractions = length(spmMask); % of they're equal, I loop through masks
    maskOpt = 0;
end
spmMaskFiles = {};
spmValsFiles = {};


for iExtraction = 1:nExtractions
    % compiles the different mask files
    if ((maskOpt~=1) || (iExtraction == 1))
        spmMaskFiles{iExtraction} = spmMask{iExtraction};
    elseif (iExtraction > 1)
        spmMaskFiles{iExtraction} = spmMaskFiles{1}; 
    end
    % compiles the subject files
    if ((maskOpt~=2) || (iExtraction == 1))
        spmValsFiles{iExtraction}  = spmVals{iExtraction};
    elseif (iExtraction > 1)
        spmValsFiles{iExtraction}  = spmValsFiles{1}; % ATTENTION, changed/ back
    end
end


% Now extract the files.
if printText
    fprintf('Subject NVoxels Mean Variance\n');
end

%% iExtraction IS subject number when using same masks across subs
for iExtraction = 1:nExtractions
    
    % this gets me the mask file name for a given extraction (which is
    % strange bc these are supposed to be subjects)
    if strcmp(spmMaskFiles{iExtraction}(end-2:end),'img')
        maskHdr = spm_vol(spmMaskFiles{iExtraction});
    elseif strcmp(spmMaskFiles{iExtraction}(end-2:end),'nii')
        maskHdr = spm_vol_nifti(spmMaskFiles{iExtraction});
    else
        disp('Unrecognizable mask format! (not img or nii)');
    end
    
    % reads in the the mask
    maskVol = spm_read_vols(maskHdr);
        
    if printText
        fprintf('Extraction #%d\n',iExtraction);
    end
    
    %% iSubject IS NOT subject no. when using same masks across subs, instead fixed at 1
    for iSubject = 1:size(spmValsFiles{iExtraction},1);
           
        %%% Not running contrast based masks

            maskIDX = find(maskVol);

        try
            if strcmp(spmValsFiles{iExtraction}(end-2:end),'img')
                valsHdr = spm_vol(spmValsFiles{iExtraction}(iSubject,:));
            elseif strcmp(spmValsFiles{iExtraction}{1}(end-3:end),'.nii')
                valsHdr = spm_vol_nifti(spmValsFiles{iExtraction}{iSubject});
            else
                disp('Unrecognizable mask format! (not img or nii)');
            end            
            valsVol = spm_read_vols(valsHdr);
            meanVal = nanmean(valsVol(maskIDX));
            variVal = nanvar(valsVol(maskIDX));
            nVoxels = length(maskIDX);
            
            AllmeanVal(iExtraction,iSubject) = meanVal;
            AllvariVal(iExtraction,iSubject) = variVal;
            AllnVoxels(iExtraction,iSubject) = nVoxels;
        catch me
            keyboard
        end
        
        [d1 fnametoprint d2] = fileparts(valsHdr.fname);
        if printText
            fprintf('%d %03d %03d %+6.4f %+6.4f %s\n',iSubject, iExtraction,nVoxels,meanVal,variVal,fnametoprint);
            iSubject
        end
    end
end

fprintf('\nFinished extracting values.\n');

