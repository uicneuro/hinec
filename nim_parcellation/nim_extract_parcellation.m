function nim_extract_parcellation()

fslmaths $FSLDIR/data/atlases/HarvardOxford/HarvardOxford-cort-maxprob-thr0-1mm.nii.gz -thr 0 -bin parcellation_mask.nii.gz
end