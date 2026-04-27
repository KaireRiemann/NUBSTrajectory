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
  \[
  p = 2s - 1
  \]
- Control point count:
  \[
  N_c = M + 2s - 1
  \]
  where `M` is the number of trajectory segments.
- Banded linear system for recovering B-spline control points
- Evaluation of position and derivatives
- Energy computation:
  \[
  E = \int_0^T \|p^{(s)}(t)\|^2 dt
  \]
- Gradient propagation interfaces for waypoint and time optimization
- Convex-hull-friendly representation for safety and dynamic constraints

---

## Why NUBS?

Classical minimum-control-effort trajectories are often represented by piecewise polynomials. This repository instead uses a **non-uniform B-spline basis** while preserving a similar low-dimensional parameterization:

\[
(P_{\text{inner}}, T) \longrightarrow C \longrightarrow p(t)
\]

where:

- \(P_{\text{inner}}\) are intermediate waypoints,
- \(T\) are segment durations,
- \(C\) are B-spline control points,
- \(p(t)\) is the continuous-time trajectory.

The key idea is to keep the external optimization variables physically meaningful, while using B-spline control points internally.

---

## Mathematical Formulation

For system order \(s\), the B-spline degree is chosen as

\[
p = 2s - 1.
\]

Given \(M\) time segments, the number of control points is

\[
N_c = M + 2s - 1.
\]

The trajectory is represented as

\[
p(t) = \sum_{i=0}^{N_c-1} N_{i,p}(t; u(T)) C_i,
\]

where:

- \(N_{i,p}\) is the B-spline basis function,
- \(u(T)\) is the open-clamped non-uniform knot vector generated from segment durations,
- \(C_i\) are control points.

The control points are recovered from the linear system

\[
A(T) C = b(P),
\]

where \(A(T)\) encodes boundary derivative constraints and waypoint interpolation constraints, and \(b(P)\) contains the prescribed boundary states and intermediate waypoints.

The constraint layout is:

1. Start boundary derivatives:
   \[
   p^{(d)}(0), \quad d = 0,\dots,s-1
   \]

2. Intermediate waypoint interpolation:
   \[
   p(t_i)=P_i, \quad i=1,\dots,M-1
   \]

3. Terminal boundary derivatives:
   \[
   p^{(d)}(T_{\Sigma}), \quad d=s-1,\dots,0
   \]

Since B-spline basis functions have local support, \(A(T)\) is a banded matrix.

---

## Knot Vector

The knot vector is generated directly from physical segment durations:

\[
T = [T_1, T_2, \dots, T_M].
\]

The cumulative time nodes are

\[
0,\quad T_1,\quad T_1+T_2,\quad \dots,\quad \sum_i T_i.
\]

The final knot vector is open-clamped:

\[
u_0=\cdots=u_p=0,
\]

\[
u_{p+i}=\sum_{j=1}^{i}T_j,
\]

\[
u_{N_c}=\cdots=u_{N_c+p}=\sum_{j=1}^{M}T_j.
\]

This makes the trajectory naturally parameterized by real physical time.

---

## Minimum-Control-Effort Energy

The default energy term is

\[
E = \int_0^{T_\Sigma} \|p^{(s)}(t)\|^2 dt.
\]

Typical cases:

| System order `s` | Degree `p=2s-1` | Energy |
| --- | --- | --- |
| 2 | 3 | minimum acceleration |
| 3 | 5 | minimum jerk |
| 4 | 7 | minimum snap |

For `s = 3`, the trajectory is a quintic non-uniform B-spline and minimizes jerk energy under the imposed interpolation and boundary constraints.

---

## Convex Hull Property

A B-spline curve has local convex-hull containment.

For each knot span, the trajectory segment lies inside the convex hull of its active control points:

\[
p(t) \in \mathrm{conv}\{C_{i-p}, \dots, C_i\}.
\]

This property is useful for trajectory safety and dynamic feasibility:

- If all active control points of a span lie in the same convex safe set, the entire span lies in that safe set.
- Derivative trajectories are also B-spline curves, so velocity, acceleration, and jerk constraints can be conservatively enforced through derivative control points.

This repository focuses on the trajectory-level primitives needed for such constraints. A full motion planner can build corridor, ESDF, or occupancy-grid cost functions on top of these primitives.

---

## API Sketch

### Include

```cpp
#include "NUBSTrajectory.hpp"