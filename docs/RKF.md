## Embedded Runge–Kutta (RKF / Dormand–Prince–type) Method

### 1. Problem Statement

We seek a numerical solution to the initial value problem (IVP)
for a (possibly vector-valued) ODE:
\[
\frac{dy}{dx} = f(x, y), \quad y(x_0) = y_0,
\]
where \(y \in \mathbb{R}^m\) and \(f: \mathbb{R} \times \mathbb{R}^m \to \mathbb{R}^m\).

### 2. Explicit s-Stage Runge–Kutta Method

An explicit \(s\)-stage Runge–Kutta (RK) method advances the solution
from \((x_n, y_n)\) to \((x_{n+1}, y_{n+1})\) with step size \(h\) by:

1. Compute stages \(k_i\):
   \[
   \begin{aligned}
   k_1 &= f(x_n, y_n), \\
   k_2 &= f\!\bigl(x_n + c_2 h,\; y_n + h a_{21} k_1 \bigr), \\
   k_3 &= f\!\bigl(x_n + c_3 h,\; y_n + h (a_{31}k_1 + a_{32}k_2) \bigr), \\
       &\quad \vdots \\
   k_s &= f\!\bigl(x_n + c_s h,\; y_n + h (a_{s1}k_1 + \cdots + a_{s,s-1}k_{s-1}) \bigr).
   \end{aligned}
   \]

2. Update solution:
   \[
   y_{n+1} = y_n + h \sum_{i=1}^s b_i k_i.
   \]

The coefficients \(a_{ij}, b_i, c_i\) define the scheme and are usually
presented in a **Butcher tableau**:

\[
\begin{array}{c|cccc}
c_1 & 0      &        &        &        \\
c_2 & a_{21} & 0      &        &        \\
c_3 & a_{31} & a_{32} & 0      &        \\
\vdots & \vdots & \vdots & \ddots & \vdots \\
c_s & a_{s1} & a_{s2} & \cdots & a_{s,s-1} \\
\hline
    & b_1    & b_2    & \cdots & b_s
\end{array}
\]

For explicit methods we have \(a_{ij} = 0\) for \(j \ge i\).

---

### 3. Embedded RK Pair (RKF / RK5(4) Type)

An **embedded Runge–Kutta pair** uses a single set of stages \(k_i\)
but two different sets of weights:

- One gives an approximation \(\tilde{y}_{n+1}\) of order \(p\) (lower order),
- The other gives an approximation \(\hat{y}_{n+1}\) of order \(q = p+1\) (higher order).

Given stages \(k_i\) from above, we form:
\[
\begin{aligned}
\hat{y}_{n+1} &= y_n + h \sum_{i=1}^s \hat{b}_i k_i \quad \text{(higher order, e.g. 5th)}, \\
\tilde{y}_{n+1} &= y_n + h \sum_{i=1}^s b_i k_i \quad \text{(lower order, e.g. 4th)}.
\end{aligned}
\]

Here \(\hat{b}_i\) and \(b_i\) are the two sets of weights (often called
“embedded weights”). Both approximations reuse the same \(k_i\).

In Dormand–Prince RK5(4)7M (the main method in the paper), we have:

- \(s = 7\) stages,
- a 5th-order formula (weights \(\hat{b}_i\)),
- an embedded 4th-order formula (weights \(b_i\)),
- the coefficients are chosen to:
  - satisfy the order conditions up to order 5,
  - minimize the principal truncation term for the 5th-order solution,
  - and provide a reasonably large region of absolute stability.

---

### 4. Local Error Estimate

The **difference** between the two embedded solutions is used
as an estimate of the local truncation error:

\[
e_{n+1} = \hat{y}_{n+1} - \tilde{y}_{n+1}.
\]

For an RK5(4) pair, this difference is \(O(h^5)\) and is therefore
a practical estimator of the local error incurred by the 4th-order solution.

In practice, for vector-valued \(y\), a norm is used, e.g.:

\[
\text{err} = \| e_{n+1} \|.
\]

The norm can be componentwise max-norm, Euclidean norm, or a weighted
norm (e.g. combining absolute and relative tolerances per component).

---

### 5. Adaptive Step Size Control

Given a user-specified tolerance `tol`, we compare:

- If `err <= tol`:
  - Accept the step:
    \[
    y_{n+1} \gets \hat{y}_{n+1} \quad (\text{use higher-order solution}),
    \]
  - Optionally increase the step size for the next step.
- If `err > tol`:
  - Reject the step:
    \[
    y_{n+1} \text{ not updated},
    \]
  - Decrease the step size and recompute.

A typical step size update rule for a RK5(4) pair is:

\[
h_{\text{new}} = \text{safety} \cdot h \left( \frac{\text{tol}}{\text{err}} \right)^{1/5},
\]

where:
- `safety` is a factor slightly less than 1 (e.g. 0.8–0.9),
- the exponent \(1/5\) corresponds to the fact that the local error
  scales like \(h^5\) for the 4th-order approximation.

Commonly, \(h_{\text{new}}\) is also constrained, e.g.:

\[
h_{\text{new}} = \max(h_{\min}, \min(h_{\max}, h_{\text{new}})),
\]
and sometimes limited in growth (e.g. at most doubled per step).

---

### 6. Summary of the RKF / Dormand–Prince Step

Given \((x_n, y_n)\) and current step size \(h\):

1. Compute stages \(k_1, \dots, k_s\) using coefficients \(a_{ij}, c_i\).
2. Compute:
   - higher-order solution \(\hat{y}_{n+1}\) using weights \(\hat{b}_i\),
   - lower-order solution \(\tilde{y}_{n+1}\) using weights \(b_i\).
3. Compute error estimate:
   \[
   e_{n+1} = \hat{y}_{n+1} - \tilde{y}_{n+1},
   \quad \text{err} = \| e_{n+1} \|.
   \]
4. If `err <= tol`:
   - accept step: \(x_{n+1} = x_n + h\), \(y_{n+1} = \hat{y}_{n+1}\).
5. Choose new step size:
   \[
   h_{\text{new}} = \text{safety} \cdot h \left( \frac{\text{tol}}{\text{err}} \right)^{1/5},
   \]
   with appropriate bounds and growth limits.
6. If the step was rejected, retry with \(h_{\text{new}}\).

This provides **automatic error control** and **adaptive step size**
using a single set of stage evaluations per step.
