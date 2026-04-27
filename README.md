# NUBSTrajectory

A lightweight C++ header-only implementation of **non-uniform B-spline (NUBS) trajectories** for minimum-control-effort trajectory generation.

> **NUBS** means **Non-Uniform B-Spline**. It is not NURBS: no rational weights are used.

## Features

- Construct non-uniform B-spline trajectories from boundary states, intermediate waypoints, and segment durations.
- Support common minimum-control-effort system orders:
  - `s = 2`: minimum acceleration
  - `s = 3`: minimum jerk
  - `s = 4`: minimum snap
- Use degree:

$$
p = 2s - 1
$$

- Use control point count:

$$
N_c = M + 2s - 1
$$

where `M` is the number of trajectory segments.

- Evaluate position and derivatives.
- Compute minimum-control-effort energy:

$$
E = \int_0^T \|p^{(s)}(t)\|^2 dt
$$

- Propagate energy gradients to intermediate points and segment times.
- Provide both runtime-order and fixed-order APIs.
- Include focused tests against MINCO construction, energy, and gradient propagation.

## Basic Form

The trajectory is represented as

$$
p(t) = \sum_{i=0}^{N_c-1} N_{i,p}(t; u(T)) C_i,
$$

where $N_{i,p}$ is the B-spline basis function, $u(T)$ is the open-clamped non-uniform knot vector generated from segment durations, and $C_i$ are control points.

The control points are recovered from a structured linear system:

$$
A(T) C = b(P).
$$

Here `A(T)` is assembled from the boundary and waypoint constraints, while `b(P)` contains the prescribed boundary states and intermediate waypoints. Since the B-spline basis has local support, this system is banded.

## Smaller Forward Construction System

For a trajectory with `M` segments and system order `s`, NUBSTrajectory solves for

$$
N_c = M + 2s - 1
$$

B-spline control points per dimension.

In a MINCO-style piecewise polynomial construction, each segment has `2s` polynomial coefficients per dimension, so the forward linear system size is

$$
2Ms.
$$

Therefore, for the same `M` and `s`, the NUBS forward construction system is smaller:

$$
M + 2s - 1 \quad \text{vs.} \quad 2Ms.
$$

This is one practical reason to use the B-spline representation: the external trajectory still satisfies the same boundary and waypoint constraints, while the forward construction solves a lower-dimensional banded system.

## Why It Is Equivalent to MINCO

MINCO usually represents the trajectory with piecewise polynomial coefficients. This project represents the same minimum-control-effort trajectory with a non-uniform B-spline basis.

For the same system order `s`, segment times, boundary states, and intermediate waypoint constraints, both formulations solve the same constrained minimum-control-effort problem. The difference is the internal basis:

- MINCO uses piecewise polynomial coefficients.
- NUBSTrajectory uses non-uniform B-spline control points.

Because both bases span the same polynomial trajectory space under the same constraints, the resulting trajectory, energy, and propagated gradients are numerically equivalent. The tests compare NUBS against MINCO for construction, energy, coefficient-space energy gradients, waypoint gradients, and time gradients.

## API

Include:

```cpp
#include "NUBSTrajectory.hpp"
```

Runtime-order API:

```cpp
nubs::NUBSTrajectory<3> traj(3); // Dim = 3, s = 3
```

Fixed-order API:

```cpp
nubs::CubicNUBS<3> cubic;     // NUBSTrajectoryT<3, 2>
nubs::QuinticNUBS<3> quintic; // NUBSTrajectoryT<3, 3>
nubs::SepticNUBS<3> septic;   // NUBSTrajectoryT<3, 4>
```

The fixed-order implementation uses compile-time Gauss rules, fixed-degree basis kernels, and fixed-degree matrix assembly for construction and gradient propagation. The default optimization path uses centered finite-difference time gradients with local affected-span and affected-row reduction. The analytic time-gradient path is kept mainly for validation.

## Tests

Tests are split into separate executables under `src/`.

Main groups:

- Basic construction checks for `s = 2, 3, 4`
- MINCO trajectory and energy comparisons
- MINCO energy-gradient propagation comparisons
- Centered finite-difference gradient checks
- Generic vs fixed-order equivalence checks
- Local finite-difference vs full finite-difference checks
- Small-duration robustness checks
- A small benchmark for `Dim = 3`, `s = 3`

Build and run:

```bash
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
```
## Repository Layout

- `include/NUBSTrajectory.hpp`: main implementation
- `include/gcopter/`: vendored MINCO-related headers used for comparison tests
- `include/tools/`: test helpers and MINCO adapter
- `src/`: test and benchmark entry points

## Acknowledgments

- B-spline theory provides the basis representation, derivative evaluation, and local support properties used by this implementation.
- MINCO is used as an important reference for construction and gradient propagation. See [Geometrically Constrained Trajectory Optimization for Multicopters](https://ieeexplore.ieee.org/document/9765821).
- This project was inspired by the minimum-control-effort / minimum-norm trajectory solving idea in [Bziyue/SplineTrajectory](https://github.com/Bziyue/SplineTrajectory).
