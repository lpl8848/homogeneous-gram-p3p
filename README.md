# Homogeneous Gram P3P

Research code for **Fast P3P through Homogeneous Gram Geometry and Direct Pose Recovery**.

This is an initial public research snapshot, accompanying a manuscript in preparation.
It contains the current core solver and standalone synthetic simulations. No publication
venue or acceptance is claimed.

The solver estimates camera-from-world poses from three calibrated 2D–3D correspondences:

```text
R * X[i] + t = s[i] * bearing[i],  s[i] > 0,  i = 0, 1, 2.
```

It retains unnormalized bearing Gram products, eliminates depths using a scalar
camera-triangle normal component, solves a quartic with a bounded resolvent-cubic
computation, and recovers rotation directly from corresponding triangle vectors.
The generic formulas impose no small-angle or positive-camera-z restriction. A
complementary coordinate path is retained for exceptional numerical configurations.

## Build

Requirements: C++17, CMake 3.16+, and Eigen 3.3+ (3.4 is recommended). The small
PoseLib polynomial-helper dependency is included with its original BSD license.
No OpenCV or downloaded RGB-D dataset is required.

Ubuntu/Debian:

```bash
sudo apt-get install cmake g++ libeigen3-dev
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
./build/p3p_example
```

Windows, with Visual Studio and Eigen available:

```powershell
cmake -S . -B build -A x64 -DEigen3_DIR="<Eigen installation>/share/eigen3"
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
.\build\Release\p3p_example.exe
```

MSVC builds use `/O2 /fp:precise`. Do not enable fast-math when comparing numerical
behavior. Link your application to `gram_p3p::gram_p3p`.

## Core API

```cpp
#include <gram_p3p/p3p.h>

std::vector<Eigen::Vector3d> bearings; // exactly three nonzero calibrated rays
std::vector<Eigen::Vector3d> world;    // corresponding world points, same order
std::vector<gram_p3p::Pose> poses;
int count = gram_p3p::solve(bearings, world, &poses);
// pose.R() and pose.t map world coordinates into camera coordinates.
```

Inputs must be finite; world points must be distinct and noncollinear. Bearing
vectors may have arbitrary positive magnitudes. The output vector is overwritten,
and the solver returns up to four positive-ray-scale poses. Negative camera z is
allowed. Additional correspondences are needed to select a pose in an application.
The API assumes its stated preconditions rather than adding input validation to
the timed solver path.

`include/gram_p3p/detail/gram_normal_solver.h` fixes the paper kernel's template
policy. The numerical code is taken from the frozen 2026-09-29 implementation;
packaging removes unused benchmark and pose-conversion dependencies. The remaining
coordinate and polynomial fallback code is part of the solver. Source hashes are
recorded in `SOURCE_SNAPSHOT.json`.

## Synthetic simulations

The accuracy model follows the statistical setup used in our manuscript:

- Four independent standard Gaussian values form a normalized rotation quaternion.
- Translation has three independent standard Gaussian components.
- Image coordinates are uniform in `[-1, 1]^2`; rays are `(u,v,1)` normalized.
- Euclidean ray ranges are uniform in `[0.1,10]`.
- World points are obtained by applying the inverse generating pose.
- Generating-pose error is the entrywise L1 rotation error plus L1 translation error,
  minimized over finite returned poses; failure threshold is `1e-6`.
- Returned poses are separately checked for ray positivity, rotation validity and
  reprojection consistency; duplicate valid poses use an L1 threshold of `1e-5`.

```bash
# Small, reproducible run; 100,000 inputs is the default.
./build/p3p_simulation --mode accuracy --samples 100000 --seed 1 --output results/accuracy

# Optional manuscript-scale sample count; no repeated timing loops.
./build/p3p_simulation --mode accuracy --samples 10000000 --seed 1 --output results/accuracy_10m

# 100,000 inputs per relative-height level: epsilon = 1e-2,...,1e-10.
# Seeds are 100 + exponent, matching the manuscript's seed convention.
./build/p3p_simulation --mode collinear --samples 100000 --seed 100 --output results/collinear

# Additional full-sphere check; the validation uses positive ray scale, not camera z.
./build/p3p_simulation --mode sphere --samples 100000 --seed 1 --output results/sphere
```

On Windows, use `build/Release/p3p_simulation.exe`. Each output directory contains
`summary.csv`, `errors.csv`, `failure_inputs.csv` and `config.json`. Existing result
directories are not overwritten. Error records are limited to the first 100,000
inputs by default; use `--save-samples N` to change this. The first 10,000 failure
inputs are retained for replay, while summary counts use every generated input.

Near-collinear points are constructed by moving the third point a relative height
`epsilon` away from the baseline, with a baseline fraction uniform in `[0.2,0.8]`.
Positive camera z is retained for this perspective-image experiment. Full-sphere
sampling is a separate validation mode and is not a claim about the paper's timing
distribution.

The portable driver uses explicit random draw order. It reproduces the statistical
model, not the exact input sequence of the original Windows paper harness. Standard
library random distributions can also differ across toolchains. Consequently these
runs must not be presented as bit-for-bit reproduction of the manuscript tables.
This initial release evaluates our solver; the full third-party baseline harness
and real-data pipeline are outside this snapshot.

## Numerical scope

This floating-point implementation does not guarantee recovery of every physical
branch. The manuscript experiment observed one generating-pose failure in ten
million ordinary inputs. Very small triangle heights and extreme numerical scales
remain ill conditioned. A small residual for returned poses does not prove that
all solutions have been returned. The supplied failure outputs support inspection
and further testing of these limits.

## Attribution and licensing

The formulation is related to the orientation-based P3P family, including
[Ke and Roumeliotis (CVPR 2017)](https://doi.org/10.1109/CVPR.2017.491).
Direct pose computation alone is not claimed as a new idea. The research concerns
the retained homogeneous quantities, depth lifting, dual-frame recovery and their
numerical evaluation.

Polynomial fallback helpers in `third_party/poselib` come from
[PoseLib](https://github.com/PoseLib/PoseLib) and retain the original BSD-3-Clause
license and copyright notices. Eigen is an external dependency under its upstream
licenses. See `THIRD_PARTY_NOTICES.md` for details.

A license for the original project code has not been assigned in this initial
snapshot. Third-party licenses continue to apply to their respective files.
Citation metadata will be added with the formal manuscript release.
