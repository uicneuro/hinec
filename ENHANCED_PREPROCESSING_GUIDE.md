# Enhanced Preprocessing with Field Map Correction

This guide explains how to use the enhanced preprocessing pipeline to fix tractography edge artifacts through proper distortion correction.

## 🎯 **Problem Solved**

The original tractography showed:
- ❌ **Edge artifacts**: Random short tracks at brain boundaries
- ❌ **Missing connections**: Gaps in major white matter tracts (corpus callosum)
- ❌ **Poor data quality**: Corrupted tensor estimation from preprocessing issues

## 🚀 **Enhanced Solution**

### **Root Cause Analysis**
1. **Missing eddy current correction** → Volume misalignment → Corrupted tensors
2. **No susceptibility correction** → Spatial distortions → Edge artifacts
3. **Poor seeding strategy** → Boundary contamination → Spurious tracks

### **Enhanced Pipeline Features**
- ✅ **Field map distortion correction** (susceptibility artifacts)
- ✅ **Alternative eddy current correction** (volume alignment)
- ✅ **White matter seeding masks** (boundary protection)
- ✅ **Improved tissue segmentation** (quality control)

## 📋 **Quick Start**

### **Step 1: Enhanced Preprocessing**

```matlab
% Configure enhanced preprocessing
options = struct();
options.run_fieldmap_correction = true;
options.fieldmap_file = 'ISMRM/ISMRM_fmap_Hz.nii.gz';  % Your field map
options.phase_encoding_dir = 'y';  % Adjust for your acquisition
options.dwell_time = 0.00058;      % Adjust for your scanner
options.eddy_method = 'eddy_correct';  % Fallback method
options.create_wm_mask = true;
% Run enhanced preprocessing
nim_preprocessing_enhanced('ISMRM/ISMRM', options);
```

### **Step 2: DTI Analysis**
```matlab
% Run DTI calculation with enhanced data
main('ISMRM/ISMRM.nii.gz', 'ISMRM_enhanced.mat');
```

### **Step 3: Enhanced Tractography**
```matlab
% Run tractography with improved seeding
runTractography('ISMRM_enhanced.mat');
```

### **Step 4: Visualize Results**
```matlab
% Command-line slice viewer
visualizeTractographySlices('tractography_results/tracks_standard.mat', 'ISMRM_enhanced.mat', 64, 64, 32);

% Or use Python GUI
% python tractography_slice_gui.py
```

## 🔧 **Configuration Options**

### **Field Map Correction Parameters**

```matlab
options.run_fieldmap_correction = true;           % Enable correction
options.fieldmap_file = 'path/to/fieldmap.nii.gz'; % Field map file
options.fieldmap_units = 'Hz';                    % 'Hz' or 'rad/s'
options.phase_encoding_dir = 'y';                 % 'x', 'y', 'z', 'x-', 'y-', 'z-'
options.dwell_time = 0.00058;                     % Echo spacing (seconds)
```

**Common dwell time values by scanner:**
- **Siemens**: ~0.00058s (typical)
- **GE**: ~0.000476s (typical)
- **Philips**: ~0.000694s (typical)

### **Eddy Correction Options**

```matlab
options.eddy_method = 'auto';        % Auto-detect best method
options.eddy_method = 'eddy_correct'; % Simple fallback (no acqp/index needed)
options.eddy_method = 'eddy';        % Advanced method (requires acqp/index)
```

### **Enhanced Seeding Options**

```matlab
options.create_wm_mask = true;```

## 📁 **File Requirements**

### **Required Files**
- `{prefix}_raw.nii.gz` - Raw DWI data
- `{prefix}.bvec` - B-vectors
- `{prefix}.bval` - B-values
- `{prefix}_fmap_Hz.nii.gz` or `{prefix}_fmap_RadPerSec.nii.gz` - Field map

### **Generated Files**
- `{prefix}.nii.gz` - Corrected DWI data
- `{prefix}_WM_mask.nii.gz` - White matter seeding mask
- `{prefix}_enhanced_preprocessing_report.mat` - Processing report

## 🎯 **Expected Improvements**

### **Before Enhancement**
- Edge artifacts from uncorrected distortions
- Missing tract connections from poor data quality
- Random short tracks at boundaries

### **After Enhancement**
- ✅ **Reduced edge artifacts** from field map correction
- ✅ **Restored major tract connectivity** from better alignment
- ✅ **Cleaner boundaries** from white matter seeding
- ✅ **Improved fiber coherence** from distortion correction

## 🔍 **Troubleshooting**

### **Common Issues**

**Field map not found:**
```
⚠ No field map found. Field map correction will be skipped.
```
**Solution:** Provide explicit field map path in `options.fieldmap_file`

**Eddy correction fails:**
```
⚠ Eddy correction failed: ...
```
**Solution:** Use fallback method: `options.eddy_method = 'eddy_correct'`

**FSL not found:**
```
❌ FSLDIR environment variable is not set
```
**Solution:** Install FSL and set environment: `export FSLDIR=/usr/local/fsl`

### **Parameter Tuning**

**If phase encoding is wrong:**
- Check acquisition protocol
- Try opposite direction (e.g., `'y'` → `'y-'`)

**If dwell time is wrong:**
- Check scanner parameters
- Look in DICOM headers or acquisition logs

**If correction overcorrects:**
- Reduce field map smoothing
- Check field map units (Hz vs rad/s)

## 📊 **Quality Assessment**

### **Check Preprocessing Report**
```matlab
load('ISMRM_enhanced_preprocessing_report.mat');
disp(preprocessing_report.enhanced_features);
disp(preprocessing_report.warnings);
```

### **Verify Output Quality**
1. **FA maps**: Should show clear white matter structure
2. **Track count**: Should be reasonable (~10K-100K tracks)
3. **Edge artifacts**: Should be significantly reduced
4. **Major tracts**: Corpus callosum should be continuous

### **Compare Results**
Use the slice viewer to compare before/after:
```matlab
% Original results
visualizeTractographySlices('tracks_original.mat', 'original.mat', 64, 64, 32);

% Enhanced results
visualizeTractographySlices('tracks_enhanced.mat', 'enhanced.mat', 64, 64, 32);
```

## 🚀 **Advanced Usage**

### **Batch Processing**
```matlab
subjects = {'subj01', 'subj02', 'subj03'};
for i = 1:length(subjects)
    fprintf('Processing %s...\n', subjects{i});
    nim_preprocessing_enhanced(subjects{i}, options);
    main([subjects{i} '.nii.gz'], [subjects{i} '_enhanced.mat']);
    runTractography([subjects{i} '_enhanced.mat']);
end
```

### **Custom Field Map Processing**
```matlab
% For multi-echo field maps
options.fieldmap_file = 'custom_fieldmap.nii.gz';
options.TE = 0.00492;  % Echo time difference
options.smooth_sigma = 1.5;  % Reduced smoothing
```

### **Research Quality Settings**
```matlab
% Maximum quality preprocessing
options.run_fieldmap_correction = true;
options.run_denoising = true;
options.denoise_method = 'dwidenoise';
options.run_motion_correction = true;
options.eddy_method = 'eddy';  % If you have acqp/index files
```

## 📚 **References**

- **FSL Documentation**: https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/
- **FUGUE Manual**: FSL's susceptibility distortion correction
- **Eddy Manual**: FSL's eddy current correction
- **ISMRM Guidelines**: Best practices for diffusion MRI preprocessing

## 🆘 **Support**

If you encounter issues:
1. Check FSL installation: `echo $FSLDIR`
2. Verify file formats: `fslinfo your_file.nii.gz`
3. Test with simple parameters first
4. Check preprocessing report for specific errors
5. Compare with working example (ISMRM dataset)

The enhanced preprocessing pipeline addresses the fundamental data quality issues that cause tractography artifacts, providing a robust foundation for high-quality fiber tracking analysis.