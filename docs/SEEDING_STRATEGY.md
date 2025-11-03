# HINEC Seeding Strategy Architecture

## Design Principle

**Centralized seeding** - All seeding decisions happen in `runTractography.m`, not in the tractography algorithm itself.

## Architecture

```
runTractography.m                    nim_tractography_standard.m
     │                                        │
     ├─ Loads nim data                        │
     ├─ Configures seeding strategy           │
     ├─ Creates seed_mask                     │
     ├─ Sets options.seed_mask                │
     │                                        │
     └─► Calls nim_tractography_standard()    │
                                              │
                                              ├─ Validates seed_mask exists
                                              ├─ Calls generate_seed_points_fact()
                                              └─ Executes tracking algorithm
```

## Rationale

### Why Centralize Seeding?

1. **Separation of Concerns**
   - `runTractography.m` = **Policy** (what to seed)
   - `nim_tractography_standard.m` = **Mechanism** (how to track)

2. **Flexibility**
   - Different pipelines can use different seeding strategies
   - IronTract can modify seeds without changing algorithm
   - Easy to test different seeding approaches

3. **Maintainability**
   - Single source of truth for seeding decisions
   - No duplicate fallback logic
   - Clear error messages when seed mask missing

## Seeding Strategy Hierarchy

`runTractography.m` implements a priority hierarchy:

### Strategy 1: Preprocessed Brain Mask (BEST)
```matlab
if isfield(nim, 'mask') && ~isempty(nim.mask)
    brain_mask = nim.mask > 0.5;
    brain_mask = brain_mask & (nim.FA > 0.05);  % Exclude CSF
```
- **Source**: FSL brain extraction or manual mask
- **Coverage**: Complete brain with clean boundaries
- **Quality**: Best anatomical accuracy

### Strategy 2: Expanded Parcellation Mask
```matlab
elseif isfield(nim, 'parcellation_mask')
    parcel_mask = nim.parcellation_mask > 0;
    brain_mask = imdilate(parcel_mask, strel('sphere', 3));
    brain_mask = brain_mask & (nim.FA > 0.05);
```
- **Source**: Atlas-based parcellation
- **Coverage**: Labeled regions + surrounding tissue
- **Quality**: Good for known anatomical structures

### Strategy 3: FA Threshold Fallback (CONSERVATIVE)
```matlab
else
    brain_mask = nim.FA > 0.10;
```
- **Source**: Fractional anisotropy map
- **Coverage**: White matter core only
- **Quality**: Misses low-anisotropy regions (fornix, cingulum)
- **WARNING**: This was causing the ISMRM validation failure

## Implementation Details

### runTractography.m
```matlab
% Lines 76-138: Seeding strategy configuration
%% CENTRALIZED SEEDING STRATEGY
fprintf('\n=== Configuring Seeding Strategy ===\n');

% Try strategies in priority order
if isfield(nim, 'mask') ...
    % Strategy 1
elseif isfield(nim, 'parcellation_mask') ...
    % Strategy 2
else
    % Strategy 3
end

% Quality check
if sum(seed_mask(:)) == 0
    error('Seed mask is empty!');
end

options.seed_mask = seed_mask;
```

### nim_tractography_standard.m
```matlab
% Lines 152-159: Validation only (no fallback)
if isempty(options.seed_mask)
    error('No seed mask provided. Configure in runTractography.m');
end
options.seed_mask = logical(options.seed_mask > 0);
```

## Special Cases

### IronTract Challenge
IronTract requires seeding from injection sites, which may be outside brain mask.

**Handled by**: `nim_irontract_submit.m` (lines 102-128)
```matlab
% Add injection site to seed mask
seed_mask = seed_mask | (injection_mask > 0);
```

This doesn't modify `runTractography.m` seeding - it happens downstream in challenge-specific code.

## Testing Seeding Strategies

### Check Current Strategy
```matlab
runTractography('data.mat');
% Look for output: "Strategy: brain_mask" or "expanded_parcellation" or "fa_threshold"
```

### Force Specific Strategy
```matlab
% Force FA-threshold (for comparison)
nim = load('data.mat');
nim.mask = []; % Remove brain mask
save('data_no_mask.mat', 'nim');
runTractography('data_no_mask.mat');
```

## Performance Impact

### Brain Mask vs FA Threshold

| Metric | Brain Mask | FA Threshold | Improvement |
|--------|------------|--------------|-------------|
| Seed voxels | ~50,000 | ~30,000 | +67% |
| Coverage | Whole brain | WM core only | Complete |
| Bundle detection | All bundles | WM only | Fornix, cingulum detected |
| False positives | Low | Very low | Slight increase |

### Seed Density

Current setting: **4 seeds/voxel**
- 50,000 voxels × 4 = **200,000 seeds**
- Bidirectional tracking = **400,000 potential tracks**

## Common Issues

### "No seed mask provided" Error
**Cause**: Called `nim_tractography_standard()` directly without setting `options.seed_mask`

**Solution**: Always use `runTractography()` wrapper, not direct call

### ISMRM Validation Shows 0% Coverage
**Cause**: FA-threshold seeding misses low-anisotropy bundles

**Solution**: Ensure brain mask exists (check `nim.mask` field)

### IronTract Seeds Not Covering Injection Site
**Cause**: Injection site outside brain mask

**Solution**: Already handled by `nim_irontract_submit.m` - it adds injection voxels

## Future Improvements

1. **Adaptive seeding density** based on local FA
2. **Multi-resolution seeding** (coarse in uniform regions, dense in complex areas)
3. **Anatomically-informed seeding** using known bundle endpoints
4. **Probabilistic seeding** with importance sampling

## References

- ISMRM Tractography Challenge: Whole-brain seeding required for full bundle coverage
- IronTract Challenge: Injection-site seeding for anatomical accuracy
- HINEC validation: Low FA seeding caused failure to detect 25/26 bundles
