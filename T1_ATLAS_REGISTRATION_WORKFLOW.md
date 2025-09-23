# T1-Based Atlas Registration Implementation Workflow

## Overview

Implement proper T1-based atlas registration to fix parcellation mask remapping issues. This workflow will edit existing files to integrate the robust MNI→T1→DWI registration chain using FSL tools.

## Core Implementation Strategy

Replace the problematic direct atlas resampling with a proper registration chain that preserves label integrity and provides anatomically accurate atlas mapping.

## Phase 1: T1 Data Detection and Validation

### Step 1: Update File Detection Logic
- Modify main.m to detect T1 files using pattern {prefix}_T1.nii.gz
- Add T1 file validation alongside existing DWI detection
- Integrate T1 availability check into preprocessing decision logic

### Step 2: Enhance Configuration Options
- Add T1-based registration options to nim_preprocessing.m default configuration
- Include enable/disable flags for T1-based atlas registration
- Add fallback options when T1 is not available

## Phase 2: T1 Preprocessing Functions

### Step 3: Create T1 Brain Extraction Function
- Implement preproc_t1_brain_extraction.m function
- Use FSL BET with robust skull stripping parameters
- Generate T1 brain mask and brain-extracted T1 image

### Step 4: Create T1-DWI Registration Function
- Implement preproc_t1_dwi_registration.m function
- Use FSL epi_reg for boundary-based registration between T1 and DWI
- Handle EPI distortions and generate T1_to_dwi.mat transformation

### Step 5: Create T1-MNI Registration Function
- Implement preproc_t1_mni_registration.m function
- Perform linear pre-alignment using FLIRT
- Execute nonlinear registration using FNIRT
- Generate and invert warp fields for MNI→T1 transformation

## Phase 3: DWI Reference Creation

### Step 6: Enhance DWI Reference Generation
- Modify existing B0 extraction to create proper DWI reference
- Extract B0 volume after eddy correction
- Apply brain extraction and masking to create clean reference
- Ensure reference is distortion-corrected and brain-masked

## Phase 4: Atlas Registration Pipeline

### Step 7: Replace Atlas Resampling Function
- Completely rewrite preproc_atlas_resampling.m
- Remove problematic FLIRT with -usesqform approach
- Implement proper registration chain using applywarp

### Step 8: Implement Composite Transform Application
- Use applywarp to combine MNI→T1 warp with T1→DWI affine transformation
- Apply nearest neighbor interpolation to preserve label values
- Generate atlas in native DWI space with proper anatomical alignment

### Step 9: Add Label Integrity Validation
- Implement atlas label validation checks
- Verify integer label values are preserved
- Check atlas coverage and anatomical reasonableness
- Generate quality metrics for registration assessment

## Phase 5: Integration and Workflow Updates

### Step 10: Update Main Preprocessing Pipeline
- Modify nim_preprocessing.m to integrate T1-based registration
- Add conditional logic to use T1 registration when available
- Maintain backward compatibility with non-T1 datasets

### Step 11: Update Configuration and Options
- Add T1 registration parameters to preprocessing options
- Include quality thresholds and validation criteria
- Provide user control over registration method selection

### Step 12: Enhance Error Handling and Fallbacks
- Add comprehensive error handling for registration failures
- Implement fallback to original method when T1 registration fails
- Provide clear diagnostic messages for troubleshooting

## Phase 6: Quality Assurance and Validation

### Step 13: Implement Registration Quality Metrics
- Add registration quality assessment functions
- Generate transformation quality reports
- Validate atlas label preservation and anatomical accuracy

### Step 14: Update Preprocessing Report
- Enhance preprocessing report to include T1 registration details
- Document transformation files and quality metrics
- Provide registration method used and success indicators

### Step 15: Update Documentation
- Modify PREPROCESSING.md to reflect new T1-based registration
- Document T1 file requirements and naming conventions
- Update configuration examples and troubleshooting guides

## Implementation Priorities

### Critical Path Dependencies
- T1 detection must be implemented before registration functions
- DWI reference creation must be completed before atlas registration
- Registration functions must be tested individually before integration

### Parallel Development Opportunities
- T1 preprocessing functions can be developed independently
- Documentation updates can proceed alongside implementation
- Error handling can be implemented in parallel with core functions

### Quality Gates
- Each registration function must validate output file existence
- Label integrity checks must pass before accepting registration results
- Transformation quality metrics must meet defined thresholds

## Risk Mitigation

### Technical Risks
- Registration failures due to poor image quality or preprocessing
- Transform composition errors leading to incorrect atlas mapping
- Performance impact from additional registration steps

### Mitigation Strategies
- Implement robust fallback to original registration method
- Add comprehensive validation at each registration step
- Provide user control over registration method selection
- Include detailed error reporting for troubleshooting

## Success Criteria

### Functional Requirements
- T1 files are automatically detected and used when available
- Atlas registration produces anatomically accurate parcellation masks
- Label values are preserved throughout registration process
- Pipeline maintains backward compatibility with existing datasets

### Quality Requirements
- Registration quality metrics indicate successful anatomical alignment
- Atlas coverage in DWI space matches expected anatomical boundaries
- Parcellation regions maintain proper size and location relationships
- Processing time increase is acceptable for improved accuracy

## File Modification Targets

### Primary Files to Edit
- main.m: Add T1 detection and conditional registration logic
- nim_preprocessing.m: Integrate T1 registration workflow
- preproc_atlas_resampling.m: Complete rewrite for proper registration

### New Files to Create
- preproc_t1_brain_extraction.m: T1 skull stripping
- preproc_t1_dwi_registration.m: EPI registration with boundary-based alignment
- preproc_t1_mni_registration.m: Nonlinear T1 to MNI registration

### Documentation to Update
- PREPROCESSING.md: Document new T1-based registration approach
- CLAUDE.md: Update with T1 file requirements and conventions