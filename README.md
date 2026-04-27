# NUBSTrajectory

A lightweight C++ library for constructing and evaluating **non-uniform B-spline (NUBS) trajectories** with minimum-control-effort structure.

This repository focuses on the trajectory representation itself: given boundary states, intermediate waypoints, and segment durations, it constructs an open-clamped non-uniform B-spline whose control points are recovered from a structured interpolation system. The resulting trajectory supports state evaluation, energy computation, and gradient propagation interfaces for optimization.

> **NUBS** here means **Non-Uniform B-Spline**.  
> It is not NURBS: no rational weights are used.

---

## ✨ Overview

Core capabilities:

- Non-uniform B-spline trajectory construction from:
  - start boundary state
  - terminal boundary state
  - intermediate waypoints
  - segment durations
- Open-clamped knot vector generated from physical segment times
- System order based formulation:
  - `s = 2`: minimum-acceleration style trajectory
  - `s = 3`: minimum-jerk style trajectory
  - `s = 4`: minimum-snap style trajectory
- B-spline degree:
  $$
  p = 2s - 1
  $$
- Control point count:
  $$
  N_c = M + 2s - 1
  $$
  where `M` is the number of trajectory segments.
- Banded linear system for recovering B-spline control points
- Evaluation of position and derivatives
- Energy computation:
  $$
  E = \int_0^T \|p^{(s)}(t)\|^2 dt
  $$
- Gradient propagation interfaces for waypoint and time optimization
- Convex-hull-friendly representation for safety and dynamic constraints

---

## Why NUBS?

Classical minimum-control-effort trajectories are often represented by piecewise polynomials. This repository instead uses a **non-uniform B-spline basis** while preserving a similar low-dimensional parameterization:

$$
(P_{\text{inner}}, T) \longrightarrow C \longrightarrow p(t)
$$

where:

- $P_{\text{inner}}$ are intermediate waypoints,
- $T$ are segment durations,
- $C$ are B-spline control points,
- $p(t)$ is the continuous-time trajectory.

The key idea is to keep the external optimization variables physically meaningful, while using B-spline control points internally.

---

## Mathematical Formulation

For system order $s$, the B-spline degree is chosen as

$$
p = 2s - 1.
$$

Given $M$ time segments, the number of control points is

$$
N_c = M + 2s - 1.
$$

The trajectory is represented as

$$
p(t) = \sum_{i=0}^{N_c-1} N_{i,p}(t; u(T)) C_i,
$$

where:

- $N_{i,p}$ is the B-spline basis function,
- $u(T)$ is the open-clamped non-uniform knot vector generated from segment durations,
- $C_i$ are control points.

The control points are recovered from the linear system

$$
A(T) C = b(P),
$$

where $A(T)$ encodes boundary derivative constraints and waypoint interpolation constraints, and $b(P)$ contains the prescribed boundary states and intermediate waypoints.

The constraint layout is:

1. Start boundary derivatives:
   $$
   p^{(d)}(0), \quad d = 0,\dots,s-1
   $$

2. Intermediate waypoint interpolation:
   $$
   p(t_i)=P_i, \quad i=1,\dots,M-1
   $$

3. Terminal boundary derivatives:
   $$
   p^{(d)}(T_{\Sigma}), \quad d=s-1,\dots,0
   $$

Since B-spline basis functions have local support, $A(T)$ is a banded matrix.

---

## Knot Vector

The knot vector is generated directly from physical segment durations:

$$
T = [T_1, T_2, \dots, T_M].
$$

The cumulative time nodes are

$$
0,\quad T_1,\quad T_1+T_2,\quad \dots,\quad \sum_i T_i.
$$

The final knot vector is open-clamped:

$$
u_0=\cdots=u_p=0,
$$

$$
u_{p+i}=\sum_{j=1}^{i}T_j,
$$

$$
u_{N_c}=\cdots=u_{N_c+p}=\sum_{j=1}^{M}T_j.
$$

This makes the trajectory naturally parameterized by real physical time.

---

## Minimum-Control-Effort Energy

The default energy term is

$$
E = \int_0^{T_\Sigma} \|p^{(s)}(t)\|^2 dt.
$$

Typical cases:

| System order `s` | Degree $p = 2s - 1$ | Energy |
| --- | --- | --- |
| 2 | 3 | minimum acceleration |
| 3 | 5 | minimum jerk |
| 4 | 7 | minimum snap |

For `s = 3`, the trajectory is a quintic non-uniform B-spline and minimizes jerk energy under the imposed interpolation and boundary constraints.

---

## Convex Hull Property

A B-spline curve has local convex-hull containment.

For each knot span, the trajectory segment lies inside the convex hull of its active control points:

$$
p(t) \in \mathrm{conv}\{C_{i-p}, \dots, C_i\}.
$$

This property is useful for trajectory safety and dynamic feasibility:

- If all active control points of a span lie in the same convex safe set, the entire span lies in that safe set.
- Derivative trajectories are also B-spline curves, so velocity, acceleration, and jerk constraints can be conservatively enforced through derivative control points.

This repository focuses on the trajectory-level primitives needed for such constraints. A full motion planner can build corridor, ESDF, or occupancy-grid cost functions on top of these primitives.

---

## API Sketch

### Include

```cpp
#include "NUBSTrajectory.hpp"
```

### Fixed-order API

The original runtime-order API is still available:

```cpp
nubs::NUBSTrajectory<3> traj(3); // s = 3, quintic NUBS
```

For `s = 2`, `s = 3`, and `s = 4`, fixed-order aliases are provided so energy and gradient paths can use compile-time degree constants:

```cpp
nubs::CubicNUBS<3> cubic;     // NUBSTrajectoryT<3, 2>
nubs::QuinticNUBS<3> quintic; // NUBSTrajectoryT<3, 3>
nubs::SepticNUBS<3> septic;   // NUBSTrajectoryT<3, 4>
```

The fixed-order implementation uses compile-time Gauss rules, fixed-degree basis kernels, and fixed-degree matrix assembly for construction and gradient propagation. Centered-difference time gradients avoid the dense knot-Jacobian path; the analytic sensitivity path builds that Jacobian lazily only when requested.

## Tests

This repository includes separate CTest executables for focused checks:

- `test_basic_s2_d2`, `test_basic_s3_d3`, `test_basic_s4_d3`: focused construction checks
- `test_minco_s2_3d`, `test_minco_s3_3d`, `test_minco_s4_3d`: 3D comparisons against upstream MINCO
- `test_minco_gradient_s2_3d`, `test_minco_gradient_s3_3d`, `test_minco_gradient_s4_3d`: energy coefficient gradients and propagated point/time gradients against upstream MINCO
- `test_centered_gradient_s2_d2`, `test_centered_gradient_s3_d3`, `test_centered_gradient_s4_d3`: centered-difference gradient propagation checks
- `test_boundary_grid`: randomized boundary and waypoint checks for `s = 2, 3, 4` and `Dim = 1, 2, 3`
- `test_generic_specialized_equivalence`: generic/runtime-order vs fixed-order equivalence checks
- `test_external_gradient`: finite-difference gradient checks against regenerated trajectories
- `test_local_vs_full_fd`: optimized local finite-difference path vs full finite-difference reference
- `test_small_duration`: robustness check for small positive durations
- `bench_finite_diff_s3_d3`: `Dim = 3`, `s = 3` timing comparison for local and full finite-difference paths

Running the executables in `bin/` prints per-case numerical errors and timing statistics.

The MINCO comparison headers are vendored under `include/gcopter/` and keep their upstream MIT license notice. Test helper headers live under `include/tools/`.

```bash
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
```

Build files are generated under `build/`, and test executables are written to `bin/`; both directories are ignored by git.

## Acknowledgments

- B-spline theory provides the basis representation, knot-vector construction, derivative evaluation, and local support properties used by this implementation.
- MINCO is used as an important construction and gradient reference. See [Geometrically Constrained Trajectory Optimization for Multicopters](https://ieeexplore.ieee.org/document/9765821).
- The testing structure and comparison style were inspired by [Bziyue/SplineTrajectory](https://github.com/Bziyue/SplineTrajectory).
