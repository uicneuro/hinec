# HINEC High-Order Tractography Implementation Workflow

## Overview

This workflow describes how to implement high-order tractography enhancements by copying and modifying `nim_tractography_standard.m` to create `nim_tractography_hinec.m`.

**Three Core Enhancements**:
1. Trilinear interpolation of diffusion direction fields
2. Runge-Kutta 4th order (RK4) numerical integration
3. Anatomically Constrained Tractography (ACT)

**Approach**: Copy existing FACT implementation and systematically replace discrete voxel tracking with high-order methods.

---

## Phase 1: File Setup and Parameter Configuration

### Step 1.1: Create New File
**Action**: Copy `nim_tractography_standard.m` to `nim_tractography_hinec.m`
**Modification**: Update function name from `nim_tractography_standard` to `nim_tractography_hinec`

### Step 1.2: Update File Header Documentation
**Location**: Lines 1-40 (header comments)
**Modifications**:
- Change description from "FACT" to "High-Order Deterministic Tractography with ACT"
- Update algorithm description to mention interpolation, RK4, and ACT
- Remove FACT-specific boundary intersection comments
- Add references to high-order methods and ACT papers

### Step 1.3: Add New Parameters to Options Structure
**Location**: Lines 42-91 (default parameter section)
**Additions**:
- `interp_method` - Interpolation method ('trilinear' as default, 'none')
- `integration_order` - Integration method (1=Euler, 2=RK2, 4=RK4 (default))
- `act_enabled` - Enable anatomically constrained tracking (boolean)
- `wm_mask` - White matter mask for ACT seeding
- `gm_mask` - Gray matter mask for valid termination
- `csf_mask` - CSF mask for invalid termination

### Step 1.4: Update Default Values
**Modifications**:
- Change default `step_size` from 0.5 to 0.2 (smaller steps for RK4 accuracy)
- Change default `interp_method` from "none" to "trilinear"
- Set default `integration_order` to 4 (RK4)
- Set default `act_enabled` to true

---

## Phase 2: Implement Trilinear Interpolation

### Step 2.1: Create Direction Interpolation Function
**Action**: Replace `get_voxel_direction_fact()` function (lines 619-670)
**New Function**: `interpolate_direction_trilinear()`
**Purpose**: Interpolate eigenvector and FA at sub-voxel positions

**Key Changes**:
- Remove discrete voxel rounding
- Extract primary eigenvector components (x, y, z) into separate 3D arrays
- Use MATLAB `interp3()` for trilinear interpolation of each component
- Interpolate FA value at same position
- Normalize resulting direction vector
- Add bounds checking for positions near volume edges

### Step 2.2: Pre-Extract Eigenvector Components
**Location**: Lines 126-150 (data verification section)
**Addition**: Extract and store eigenvector components for efficient interpolation access
**Components**:
- Extract v1_x, v1_y, v1_z from nim.evec(:,:,:,:,1)
- Store in nim structure for repeated access during tracking
- Reduces memory access overhead during interpolation

### Step 2.3: Update Direction Retrieval Calls
**Locations**: Throughout tracking functions
**Modifications**:
- Replace all calls to `get_voxel_direction_fact()`
- Use `interpolate_direction_trilinear()` instead
- Update function signatures to pass interpolation parameters

---

## Phase 3: Implement RK4 Integration

### Step 3.1: Create RK4 Step Function
**Action**: Add new function `rk4_integration_step()`
**Location**: After helper functions (around line 550)
**Purpose**: Compute single RK4 integration step

**Algorithm Structure**:
- Calculate k1: direction at current position
- Calculate k2: direction at position + 0.5*step_size*k1
- Calculate k3: direction at position + 0.5*step_size*k2
- Calculate k4: direction at position + step_size*k3
- Combine with weighting: new_position = position + step_size*(k1 + 2*k2 + 2*k3 + k4)/6

**Requirements**:
- Each k evaluation uses trilinear interpolation
- Enforce direction consistency (flip if dot product negative)
- Handle interpolation failures at intermediate points
- Fallback to previous valid direction if needed

### Step 3.2: Add Integration Method Selection
**Action**: Create `advance_position()` function
**Purpose**: Route to appropriate integration method based on options
**Methods**:
- Order 1: Euler integration (position + step_size * direction)
- Order 2: RK2/midpoint method
- Order 4: Full RK4 integration

### Step 3.3: Replace Position Updates in Tracking Loop
**Location**: Lines 473-543 (`track_fiber_fact()` function)
**Modifications**:
- Remove FACT boundary intersection logic (`find_voxel_boundary_exit()`)
- Replace with call to `advance_position()` using selected integration method
- Update position by fixed step size instead of voxel boundary jumps
- Maintain direction consistency checking

---

## Phase 4: Implement Anatomically Constrained Tractography (ACT)

### Step 4.1: Add Tissue Mask Loading
**Location**: Lines 102-116 (data loading section)
**Additions**:
- Check for tissue masks in nim structure
- Load WM, GM, CSF masks if available
- Validate mask dimensions match diffusion data
- Create flags indicating which masks are available

### Step 4.2: Constrain Seed Generation to White Matter
**Location**: `generate_seed_points_fact()` function (lines 351-422)
**Modifications**:
- If ACT enabled and WM mask available, intersect seed_mask with WM mask
- Only generate seeds within white matter regions
- Update seed statistics reporting

### Step 4.3: Add Tissue Type Checking Function
**Action**: Create new function `check_tissue_type()`
**Purpose**: Determine tissue type at given position
**Returns**: Tissue classification (WM, GM, CSF, or outside brain)
**Logic**:
- Round position to nearest voxel
- Check bounds
- Query tissue masks in priority order
- Return tissue type enumeration

### Step 4.4: Implement ACT Termination Criteria
**Location**: `track_fiber_fact()` main tracking loop (lines 473-544)
**Additions in Loop**:
- After computing new position, check tissue type
- If CSF detected: terminate immediately (invalid endpoint)
- If GM detected: add point and terminate (valid endpoint)
- If outside brain: terminate immediately
- If WM: continue tracking normally

**Modification**: Update track storage to mark termination type
**Purpose**: Enable filtering and analysis of termination reasons

### Step 4.5: Handle Missing ACT Masks Gracefully
**Location**: Throughout ACT-related code
**Fallback Logic**:
- If ACT enabled but masks unavailable, issue warning
- Disable ACT features automatically
- Fall back to FA-based termination only
- Document mask generation requirements in warnings

---

## Phase 5: Update Tracking Algorithm Integration

### Step 5.1: Modify Main Tracking Function
**Location**: `track_fiber_fact()` function (lines 424-548)
**Rename**: Change to `track_fiber_hinec()`
**Core Modifications**:
- Replace discrete voxel-based tracking with continuous sub-voxel tracking
- Use interpolation for all direction queries
- Use RK4 for position advancement
- Apply ACT constraints at each step
- Remove voxel boundary intersection logic entirely

### Step 5.2: Update Tracking Loop Structure
**Current Loop**: While loop with boundary intersection
**New Loop**: While loop with fixed step size advancement
**Changes**:
- Remove `find_voxel_boundary_exit()` calls
- Add `advance_position()` using selected integration method
- Add `check_tissue_type()` for ACT constraints
- Keep FA and angle threshold checks
- Keep max steps termination

### Step 5.3: Direction Consistency Management
**Location**: Within tracking loop after getting new direction
**Enhancement**:
- Check dot product between consecutive directions
- Flip direction if dot product negative
- Enforce smoother transitions with interpolation
- Apply angle threshold using cosine comparison

---

## Phase 6: Remove FACT-Specific Functions

### Step 6.1: Delete Voxel Boundary Intersection
**Action**: Remove `find_voxel_boundary_exit()` function (lines 550-608)
**Reason**: Not needed for continuous high-order tracking

### Step 6.2: Update Initial Direction Function
**Action**: Modify `get_initial_direction_fact()` (lines 610-616)
**Changes**:
- Rename to `get_initial_direction_hinec()`
- Use interpolation instead of voxel lookup
- Apply ACT constraint (must start in WM if enabled)

---

## Phase 7: Update Diagnostics and Output

### Step 7.1: Modify Timing Diagnostics
**Location**: Lines 96-149, 208-213, 280-294
**Changes**:
- Remove "FACT verification" timing section
- Add "Interpolation setup" timing
- Update tracking statistics labels (remove "voxel steps", add "integration steps")
- Add RK4-specific performance metrics

### Step 7.2: Enhance Failure Analysis
**Location**: Lines 197-307 (failure tracking section)
**Additions**:
- Track ACT-based terminations separately (GM termination, CSF termination)
- Count interpolation failures
- Report termination type distribution
- Add tissue-specific success rates

### Step 7.3: Update Track Statistics
**Location**: Lines 323-342 (track statistics calculation)
**Additions**:
- Calculate percentage GM terminations (valid endpoints)
- Calculate percentage CSF terminations (invalid)
- Report average track curvature
- Add ACT effectiveness metrics

### Step 7.4: Update Console Output Messages
**Locations**: Throughout file (fprintf statements)
**Changes**:
- Replace "FACT" references with "HINEC High-Order"
- Update parameter descriptions (interpolation method, integration order)
- Add ACT status reporting
- Update algorithm description in progress messages

---

## Phase 8: Handle Edge Cases and Validation

### Step 8.1: Boundary Handling for Interpolation
**Location**: Within `interpolate_direction_trilinear()`
**Additions**:
- Check position bounds before interpolation (must have buffer for trilinear)
- Return empty direction if too close to volume edges
- Add tolerance parameter for boundary proximity
- Handle numerical precision issues

### Step 8.2: Integration Stability Checks
**Location**: Within `rk4_integration_step()`
**Additions**:
- Validate intermediate positions are within bounds
- Check for NaN or Inf values in k1, k2, k3, k4
- Implement fallback to lower-order method if instability detected
- Add step size reduction on failure

### Step 8.3: ACT Mask Validation
**Location**: After mask loading
**Checks**:
- Verify masks have same dimensions as diffusion data
- Check for non-empty masks (at least some voxels in each tissue type)
- Warn if tissue masks have suspicious characteristics (e.g., no WM voxels)
- Validate mask overlap (WM and GM should not overlap significantly)

---

## Phase 9: Optimization

### Step 9.1: Pre-Compute Interpolation Grids
**Location**: After data loading, before tracking loop
**Action**: Extract and store interpolation-ready arrays
**Components**:
- Primary eigenvector components (v1_x, v1_y, v1_z)
- FA map
- Tissue masks (if ACT enabled)
**Benefit**: Reduce memory access overhead during repeated interpolation calls

### Step 9.2: Optimize Interpolation Calls
**Location**: Within interpolation functions
**Optimizations**:
- Cache frequently accessed volume data
- Minimize array indexing operations
- Use vectorized operations where possible
- Avoid repeated bounds checking for same position

### Step 9.3: Add Early Termination Optimizations
**Location**: Main tracking loop
**Optimizations**:
- Check FA threshold before expensive RK4 calculation
- Skip RK4 if direction hasn't changed significantly
- Use simpler integration near track endpoints
- Add maximum step count per track to prevent runaway tracks

---

## Phase 10: Testing and Validation

### Step 10.1: Add Interpolation Validation
**Purpose**: Verify interpolation produces reasonable results
**Tests**:
- Compare interpolated directions at voxel centers with actual voxel directions
- Check smoothness of direction field across voxel boundaries
- Validate FA interpolation matches expected values
- Test edge cases near volume boundaries

### Step 10.2: Add RK4 Accuracy Validation
**Purpose**: Ensure integration is mathematically correct
**Tests**:
- Compare RK4 with Euler on same seed (RK4 should be smoother)
- Verify k1, k2, k3, k4 are computed at correct positions
- Check final position computation uses correct weighting
- Validate step size effects on accuracy

### Step 10.3: Add ACT Logic Validation
**Purpose**: Confirm anatomical constraints work correctly
**Tests**:
- Verify seeds are only placed in WM when ACT enabled
- Confirm tracks terminate in GM (not passing through)
- Validate CSF avoidance (no tracks entering CSF)
- Check tracks respect brain boundaries

### Step 10.4: Compare with Standard FACT
**Purpose**: Validate high-order improvements
**Comparisons**:
- Track count (HINEC should generate more valid tracks)
- Track smoothness (HINEC should have smoother curvature)
- Anatomical plausibility (HINEC should respect tissue boundaries better)
- Major bundle reconstruction (both should find major pathways)

---

## Implementation Checklist

### Core Algorithm Components
- [ ] Trilinear interpolation of eigenvector fields
- [ ] Trilinear interpolation of FA values
- [ ] RK4 integration with four-stage evaluation
- [ ] Integration method selection (Euler, RK2, RK4)
- [ ] Direction consistency enforcement

### ACT Components
- [ ] Tissue mask loading and validation
- [ ] WM-constrained seed generation
- [ ] Tissue type checking function
- [ ] GM termination criteria
- [ ] CSF avoidance logic
- [ ] Termination type tracking

### Supporting Features
- [ ] Pre-extraction of eigenvector components
- [ ] Boundary handling for interpolation
- [ ] Fallback logic for missing masks
- [ ] Enhanced diagnostics and timing
- [ ] Track statistics with ACT metrics
- [ ] Parameter validation

### Quality Assurance
- [ ] Edge case handling (volume boundaries)
- [ ] Numerical stability checks
- [ ] Mask dimension validation
- [ ] Direction field smoothness validation
- [ ] Integration accuracy verification
- [ ] ACT logic correctness testing

---

## Key Differences from FACT

**FACT (nim_tractography_standard.m)**:
- Discrete voxel-based tracking
- Voxel boundary intersection for position updates
- No interpolation (uses voxel center directions only)
- Euler integration (straight line between voxels)
- FA-based termination only

**HINEC (nim_tractography_hinec.m)**:
- Continuous sub-voxel tracking
- Fixed step size with RK4 integration
- Trilinear interpolation of direction fields
- Higher-order numerical integration (RK4)
- Multi-criteria termination (FA, angle, tissue type)
- Anatomically constrained to white matter pathways

---

## Expected Outcomes

**Track Quality Improvements**:
- Smoother fiber trajectories (less angular artifacts)
- More accurate path estimation in curved regions
- Better sub-voxel precision

**Anatomical Validity**:
- Tracks confined to white matter
- Valid termination in gray matter
- CSF regions avoided
- Biologically plausible connectivity

**Trade-offs**:
- Slower execution (more computation per step)
- Requires tissue segmentation masks for full ACT
- More complex parameter tuning
- Higher memory usage for pre-computed interpolation grids

---

## Troubleshooting Guide

**Issue: Interpolation returns empty directions**
- Check position bounds (must have buffer from edges)
- Verify eigenvector data is valid (no NaN/Inf)
- Ensure correct array indexing (MATLAB uses 1-based)

**Issue: RK4 produces unstable tracks**
- Reduce step size for better numerical stability
- Check intermediate position bounds
- Verify direction consistency at k2, k3, k4
- Consider fallback to RK2 in problematic regions

**Issue: No tracks generated with ACT**
- Verify tissue masks are loaded correctly
- Check WM mask has sufficient coverage
- Ensure masks have correct dimensions
- Try disabling ACT to isolate issue

**Issue: Tracks terminate prematurely**
- Check FA termination threshold (may be too high)
- Verify angle threshold is reasonable
- Examine tissue masks for gaps
- Review ACT termination criteria

**Issue: Poor performance**
- Pre-compute eigenvector component arrays
- Optimize interpolation calls
- Reduce seed density
- Consider using RK2 instead of RK4
