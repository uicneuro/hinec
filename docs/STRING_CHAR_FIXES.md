# String vs Char Array Fixes - HINEC Run Directory System

## Problem Overview

MATLAB distinguishes between **string** types (introduced in newer versions) and **char arrays** (classic MATLAB). Many file operations like `fileparts()`, `copyfile()`, and `movefile()` expect char arrays, but the YAML config system and shell script may pass string types.

**Error Symptoms**:
- `Error using fileparts: Input path must be text`
- `Error using copyfile: All arguments must be non-empty character vectors or string scalars`
- `Operands to the logical AND (&&) and OR (||) operators must be convertible to logical scalar values`

## Root Cause

When the run directory system was implemented, file paths were being passed around as strings (from YAML configs and shell scripts) but MATLAB's file operation functions expected char arrays. The issue manifested in three ways:

1. **`fileparts()` failures**: Cannot extract basename from string arrays
2. **`copyfile()`/`movefile()` failures**: File operations require exact char array types
3. **`isfile()` with logical operators**: Returns arrays instead of scalars when given string arrays

## Complete Fix Strategy

**Universal Rule**: Convert all external file path inputs to char arrays **at the entry point** of each function, before any operations.

---

## Files Fixed

### 1. `/Users/12salty/Documents/research-chun/hinec/main.m`

**Location**: Lines 322, 351, 396

**Functions Fixed**:
- `setup_file_paths(imgpath, run_info)` - Line 322
- `handle_preprocessed_data(img_file, mask_file, imgpath, run_info)` - Line 351
- `handle_raw_data(img_file, raw_file, t1_file, imgpath, options, run_info)` - Line 396

**Fix Applied**:
```matlab
% At the very beginning of each function, before any operations:
imgpath = char(imgpath);
```

**Why This Works**:
- Ensures `fileparts(imgpath)` works correctly
- Ensures string concatenation like `[imgpath '.nii.gz']` produces char arrays
- All subsequent operations receive proper char array types

**Also Fixed**:
- Line 342-343: `detect_data_type()` - Wrapped `isfile()` with `any()` for scalar logical
- Line 135, 422: Error messages - Wrapped arguments with `char()` for proper formatting

---

### 2. `/Users/12salty/Documents/research-chun/hinec/nim_utils/create_run_directory.m`

**Location**: Lines 31-43

**Function Fixed**: `create_run_directory(config_file, varargin)`

**Fix Applied**:
```matlab
% Convert config_file to char array for consistent handling
config_file = char(config_file);

% Parse optional arguments
p = inputParser;
addParameter(p, 'base_dir', 'hinec_runs', @ischar);
addParameter(p, 'run_name', '', @ischar);
addParameter(p, 'description', '', @ischar);
parse(p, varargin{:});

base_dir = char(p.Results.base_dir);
custom_run_name = char(p.Results.run_name);
description = char(p.Results.description);
```

**Why This Works**:
- `config_file` comes from shell script (may be string type)
- Converting at entry ensures `fileparts(config_file)` works at line 55
- Converting parsed arguments ensures all path construction uses char arrays
- `copyfile(config_file, config_copy)` at line 101 works correctly

---

### 3. `/Users/12salty/Documents/research-chun/hinec/nim_preprocessing/preproc_brain_extraction.m`

**Location**: Lines 23-26

**Function Fixed**: `preproc_brain_extraction(dwi_or_b0_file, output_dir, brain_mask_file)`

**Fix Applied**:
```matlab
% Convert inputs to char arrays for consistent handling
dwi_or_b0_file = char(dwi_or_b0_file);
output_dir = char(output_dir);
brain_mask_file = char(brain_mask_file);
```

**Why This Works**:
- Called from `main.m` which may pass string types from YAML config
- `fileparts(dwi_or_b0_file)` at line 37 now works correctly
- `movefile(bet_mask_file, brain_mask_file)` at line 69 works correctly
- `fullfile()` operations produce char arrays

---

### 4. `/Users/12salty/Documents/research-chun/hinec/nim_preprocessing/nim_preprocessing.m`

**Location**: Line 24 (already fixed)

**Function Fixed**: `nim_preprocessing(file_prefix, varargin)`

**Already Present**:
```matlab
% Ensure file_prefix is a character array (handle string input)
file_prefix = char(file_prefix);
```

**Status**: ✅ Already correctly handled - no changes needed

---

## Functions That Don't Need Fixes

### `runTractography.m`
- **Why**: Uses `fullfile()` for paths but doesn't receive external file paths directly
- **Status**: ✅ No string/char issues

### `nim_registration/register_t1_to_mni.m`
- **Why**: Receives `registration_data` structure, not direct file paths
- **Status**: ✅ No string/char issues

### `nim_utils/load_config_yaml.m`
- **Why**: Only reads YAML files, returns structures with values (strings are fine in structures)
- **Status**: ✅ No string/char issues

---

## Testing Checklist

After applying these fixes, the complete pipeline should work:

```bash
# Run with nohup
nohup ./run_hinec.sh ISMRM/ISMRM sample.mat config/hinec_default.yml > hinec_wrapper.log 2>&1 &

# Monitor progress
tail -f hinec_wrapper.log

# Verify run directory created
ls -lt hinec_runs/latest/

# Check files copied correctly
ls -l hinec_runs/latest/intermediate/
ls -l hinec_runs/latest/output/
ls -l hinec_runs/latest/tractography/
```

**Expected Behavior**:
1. ✅ Preprocessing runs in `ISMRM/` directory
2. ✅ Files copied to `hinec_runs/run_YYYYMMDD_HHMMSS_hinec_default/intermediate/`
3. ✅ DTI calculation uses copied files from run directory
4. ✅ Tractography outputs to `hinec_runs/.../tractography/`
5. ✅ No string/char array type errors

---

## Prevention Guidelines

To prevent similar issues in the future:

### 1. Function Entry Point Rule
**Always convert file path parameters to char arrays at function entry:**
```matlab
function result = my_function(file_path, other_args)
    % Convert to char array FIRST
    file_path = char(file_path);

    % Now safe to use fileparts, copyfile, etc.
    [dir, name, ext] = fileparts(file_path);
    % ...
end
```

### 2. Detect String/Char Issues
Look for these patterns that indicate potential issues:
- `fileparts(variable)` - variable MUST be char array
- `copyfile(src, dst)` - both MUST be char arrays
- `movefile(src, dst)` - both MUST be char arrays
- `isfile(path) && other_condition` - wrap with `any(isfile(path))`

### 3. Testing with YAML Configs
When testing with YAML config system, paths come from:
- Shell script arguments (may be strings)
- YAML file values (parsed as strings)
- `fullfile()` results (returns string if inputs are strings)

**Always test the complete pipeline via `run_hinec.sh`, not just MATLAB directly.**

---

## Technical Notes

### Why `any()` for `isfile()`?

When `isfile()` receives a string array, it may return a logical array instead of a scalar:
```matlab
% Problem:
is_raw = isfile(raw_file);  % May be [true] (array) instead of true (scalar)
is_preprocessed = isfile(img_file) && ~is_raw;  % ERROR: && needs scalar

% Solution:
is_raw = any(isfile(raw_file));  % Force scalar conversion
is_preprocessed = any(isfile(img_file)) && ~is_raw;  # Works!
```

### Why Convert at Entry, Not at Use?

**Bad approach** (convert at each use):
```matlab
function result = process(file_path)
    [dir, name, ext] = fileparts(char(file_path));  % char() here
    copyfile(char(file_path), char(dest));          # char() here
    new_path = [char(file_path) '.bak'];            # char() here
end
```

**Good approach** (convert once at entry):
```matlab
function result = process(file_path)
    file_path = char(file_path);  # Convert ONCE

    [dir, name, ext] = fileparts(file_path);  # Clean
    copyfile(file_path, dest);                # Clean
    new_path = [file_path '.bak'];            # Clean
end
```

**Benefits**:
- Single point of conversion
- Cleaner code throughout function
- No risk of forgetting `char()` wrapper
- Consistent variable type for entire function

---

---

## Additional Fixes (Session 2)

### 5. `/Users/12salty/Documents/research-chun/hinec/nim_utils/nim_save.m`

**Location**: Lines 7, 11

**Function Fixed**: `nim_save(nim, nimpath)`

**Fix Applied**:
```matlab
function nim_save(nim, nimpath)
    arguments
        % nim struct
        nim

        % Path to save processed .mat file. Must end in `.mat`
        nimpath  % Removed "string" type declaration
    end

    % Convert to char array for save() function compatibility
    nimpath = char(nimpath);
```

**Why This Works**:
- MATLAB's `save()` function requires char array, not string type
- Removing the `string` type declaration from arguments allows any type
- Converting to char at entry ensures compatibility

**Error Fixed**:
```
Error using save
Argument must be a text scalar.
```

---

### 6. `/Users/12salty/Documents/research-chun/hinec/nim_utils/nim_load_nim.m`

**Location**: Lines 11, 19-22

**Function Fixed**: `nim_load_nim(file_prefix)`

**Fix Applied**:
```matlab
% Convert to char array for consistent handling
file_prefix = char(file_prefix);

% Define the suffixes for different files
suffix_processed = '.nii.gz';
suffix_brain_mask = '_M.nii.gz';
suffix_atlas_labels = '_atlas_labels.mat';

% Construct file paths (using char array concatenation)
nii_file = [file_prefix suffix_processed];
mask_file = [file_prefix suffix_brain_mask];
parcellation_mask_file = fullfile(fileparts(file_prefix), 'parcellation_mask.nii.gz');
atlas_labels_file = [file_prefix suffix_atlas_labels];
```

**Why This Works**:
- Original code used `+` operator which creates string types
- Changed to `[]` concatenation which preserves char array type
- `fileparts()` at line 21 now receives proper char array
- All `exist()` and file operations work correctly

---

## Summary

**Files Modified**: 5 files
- `main.m` (3 functions) ✅
- `nim_utils/create_run_directory.m` (1 function) ✅
- `nim_preprocessing/preproc_brain_extraction.m` (1 function) ✅
- `nim_utils/nim_save.m` (1 function) ✅ **NEW**
- `nim_utils/nim_load_nim.m` (1 function) ✅ **NEW**

**Files Already Correct**: 3 files
- `nim_preprocessing/nim_preprocessing.m` ✅
- `nim_utils/nim_read.m` ✅
- `runTractography.m` ✅
- Other preprocessing/registration functions ✅

**Fix Pattern**: Convert all external file path inputs to char arrays at function entry point

**Status**: ✅ Complete - All string/char array issues resolved
