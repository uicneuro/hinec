# High-Order Tractography Enhancements for MATLAB FACT Pipeline

This report describes three major improvements for extending a basic FACT (Fiber Assignment by Continuous Tracking) tractography algorithm into a **high-order tractography** method. The upgrades include **interpolation**, **Runge–Kutta 4th order integration**, and **anatomically constrained tractography (ACT)**. Each feature enhances realism, accuracy, and anatomical validity.

---

## 1. Interpolation of Diffusion Data

### **Concept and Purpose**
FACT normally uses voxel-centered principal diffusion directions. However, when step sizes are smaller than voxel dimensions, the streamline may cross between voxels with abrupt changes in orientation.  
**Interpolation** allows estimating the diffusion direction at sub-voxel positions, ensuring smooth transitions and continuous trajectories.

### **Common Methods**
- **Trilinear interpolation**: Uses eight neighboring voxels to estimate direction (most common).
- **Spherical interpolation**: Applied when dealing with ODFs (Orientation Distribution Functions).
- **Tensor interpolation (Log-Euclidean)**: Smoothly interpolates tensors instead of eigenvectors directly.

### **Implementation Outline (MATLAB)**
1. Load FA map and eigenvector fields.  
2. For each tracking step:
   - Compute new position:  
     ```matlab
     x_next = x + step_size * dir;
     ```
   - Interpolate the direction field at `x_next` using:
     ```matlab
     dir_interp = interp3(evec_x, evec_y, evec_z, x_next, y_next, z_next);
     dir = dir_interp / norm(dir_interp);
     ```
3. Use interpolated directions for subsequent integration steps.

---

## 2. Runge–Kutta 4th Order (RK4) Integration

### **Concept and Purpose**
The basic FACT uses **Euler integration**, which assumes constant direction between steps and can deviate from true paths in curved regions.  
The **RK4 method** improves numerical accuracy by evaluating intermediate points within each step, effectively producing smoother, higher-order trajectories.

### **Mathematical Formulation**
\[
\begin{aligned}
k_1 &= v(r_n), \\
k_2 &= v(r_n + \tfrac{1}{2}h k_1), \\
k_3 &= v(r_n + \tfrac{1}{2}h k_2), \\
k_4 &= v(r_n + h k_3), \\
r_{n+1} &= r_n + \tfrac{h}{6}(k_1 + 2k_2 + 2k_3 + k_4).
\end{aligned}
\]

### **Implementation Outline (MATLAB)**
1. At each streamline point `r_n`:
   ```matlab
   k1 = interpolate_direction(r_n);
   k2 = interpolate_direction(r_n + 0.5*h*k1);
   k3 = interpolate_direction(r_n + 0.5*h*k2);
   k4 = interpolate_direction(r_n + h*k3);
   r_next = r_n + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
    ```

2. Normalize the final direction vector.
3. Apply stopping criteria (FA threshold, curvature angle, or mask limits).

---

## 3. Anatomically Constrained Tractography (ACT)

### **Concept and Purpose**

ACT introduces **tissue-based constraints** to ensure streamlines follow anatomically valid white matter (WM) pathways and terminate in gray matter (GM), avoiding CSF or non-brain regions. This approach improves biological realism.

### **Required Maps**

* WM mask
* GM mask
* CSF mask (optional: boundary masks WM–GM and WM–CSF)
  These can be generated from T1-weighted MRI segmentation (e.g., FSL FAST or Freesurfer).

### **Implementation Outline (MATLAB)**

1. Initialize seed points within WM.
2. During propagation:

   ```matlab
   if in_CSF(r_next, csf_mask)
       break; % invalid termination
   elseif in_GM(r_next, gm_mask)
       streamline = [streamline; r_next];
       break; % valid termination
   end
   ```
3. Exclude streamlines leaving the brain mask or entering CSF.

---

## Integration Summary

| Step | Feature         | Function                                       | Dependency                        |
| ---- | --------------- | ---------------------------------------------- | --------------------------------- |
| 1    | Interpolation   | Smooth direction estimation at sub-voxel level | Replaces discrete voxel lookup    |
| 2    | RK4 Integration | Higher-order streamline propagation            | Uses interpolated direction field |
| 3    | ACT             | Anatomical plausibility enforcement            | Uses segmentation masks           |

---

## Example MATLAB Skeleton

```matlab
for seed = seed_points
    r = seed;
    streamline = r;
    while is_in_brain(r, brain_mask)
        % Runge–Kutta 4th order integration
        k1 = interpolate_direction(r, evecs);
        k2 = interpolate_direction(r + 0.5*h*k1, evecs);
        k3 = interpolate_direction(r + 0.5*h*k2, evecs);
        k4 = interpolate_direction(r + h*k3, evecs);
        r_next = r + (h/6)*(k1 + 2*k2 + 2*k3 + k4);

        % Apply ACT constraints
        if in_CSF(r_next, csf_mask)
            break;
        elseif in_GM(r_next, gm_mask)
            streamline = [streamline; r_next];
            break;
        end

        streamline = [streamline; r_next];
        r = r_next;
    end
    save_streamline(streamline);
end
```

---

## Summary

By combining:

* **Interpolation** for continuous direction fields,
* **RK4 integration** for precise streamline propagation, and
* **Anatomically Constrained Tractography (ACT)** for biological plausibility,

the MATLAB FACT tractography pipeline becomes a **high-order deterministic tractography system** suitable for advanced diffusion MRI analysis.

```