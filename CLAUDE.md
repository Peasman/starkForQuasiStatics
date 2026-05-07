# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repo is

A fork of the [Stark](https://github.com/InteractiveComputerGraphics/stark) physics simulation framework, modified to be more suitable for quasi-static simulation. Key changes from upstream:

- PSD projection now supports absolute eigenvalue clamping (see `Settings::NewtonsMethod::project_to_PD` and the referenced [paper](https://www.cs.columbia.edu/cg/abs-psd/paper.pdf)).
- Newton solver convergence is based on the max-norm of the displacement step (`du_norm_threshhold`) instead of only residual changes.
- Added MKL-backed direct solvers: `LinearSystemSolver::MKL_LU` and `LinearSystemSolver::MKL_LDLT` (Pardiso).

## Build

Requires CMake 3.15+, a C++17 compiler, OpenMP, and Intel MKL (linked as `MKL::MKL`).

```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . -j$(nproc)
```

Build targets produced:
- `stark` — the core static library
- `examples` — C++ example executable (`examples/main.cpp`)
- `stark_tests` — Catch2 test executable (`tests/rb_constraints.cpp`)
- `pystark` — Python extension module (nanobind, Python 3.8+), built when `STARK_ENABLE_PYTHON_BINDINGS=ON`

To disable AVX (e.g. for ARM): `-DSTARK_ENABLE_AVX=OFF`

### Running tests

```bash
./build/stark_tests
```

### Running examples

```bash
./build/examples
```

## Architecture

```
stark/src/
  core/           # Low-level engine
  models/         # High-level simulation objects
  utils/          # Mesh helpers
  extern/         # Bundled deps (symx, Eigen, TMCD, fmt, …)
```

### Core (`stark/src/core/`)

| File | Role |
|------|------|
| `Stark.{h,cpp}` | Owns the time loop, calls Newton's method each step |
| `NewtonsMethod.{h,cpp}` | Newton iteration + line search + linear solve dispatch |
| `Settings.{h,cpp}` | All tunable parameters (see below) |
| `Callbacks.h` | Hooks for before/after energy evaluation and validity checks |
| `EventDrivenScript.{h,cpp}` | Time-event scripting system |

**Time integration**: DOFs are *velocities* (`v1`), not displacements. Positions are computed as `x1 = x0 + dt * v1`.

### Models (`stark/src/models/`)

`Simulation` is the user-facing entry point. It owns:
- `Deformables` — point cloud + FEM/shell energies
- `RigidBodies` — rigid body dynamics and constraints
- `Interactions` — frictional contact (IPC) and attachments
- `Presets` — convenience factories for common objects

Each energy term is expressed symbolically via **symx**, which auto-generates and JIT-compiles gradient/Hessian code at first use. A C++17 compiler must be accessible at runtime for this.

### Python bindings (`pystark/`)

`pystark/cpp/` contains nanobind wrappers that mirror the C++ API 1-to-1. `pystark/examples/` has runnable Python scenes.

## Key settings

`stark::core::Settings` (used via `stark::Settings` at the user level):

```cpp
settings.newton.linear_system_solver    // CG | DirectLU | MKL_LU | MKL_LDLT
settings.newton.use_du_norm_threshhold  // converge when max|du|*dt < du_norm_threshhold
settings.newton.du_norm_threshhold      // default 1e-3
settings.newton.use_residual_threshhold // converge when max residual < residual.tolerance
settings.newton.residual                // {ResidualType::Force|Acceleration, tolerance}
settings.newton.project_to_PD          // enable PSD projection of element Hessians
settings.debug.symx_suppress_compiler_output  // set false to debug codegen
settings.output.codegen_directory       // cache for generated/compiled energy code
```

When both `use_du_norm_threshhold` and `use_residual_threshhold` are true, convergence requires *both* conditions to hold simultaneously.

## Adding new energy models

1. Define the energy symbolically in `stark/src/models/` using the symx API (see existing energies like `EnergyTriangleStrain.cpp` for the pattern).
2. Register it with `stark.global_energy` and hook callbacks via `stark.callbacks`.
3. Expose via the relevant aggregate class (`Deformables`, `RigidBodies`, etc.).
4. Add nanobind wrappers in `pystark/cpp/` to expose to Python.
