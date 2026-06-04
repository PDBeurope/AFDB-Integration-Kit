## ipSAE C++ implementation

### Dependencies

- A C++17 compiler with OpenMP (`g++` recommended)
- `curl` and `tar` (used by `make` to fetch Eigen on first build)
- [Eigen](https://eigen.tuxfamily.org/) **3.4.0** — fetched on demand into `deps/eigen-3.4.0/` (not vendored in this repo)
- [nlohmann/json](https://github.com/nlohmann/json) — vendored as the single-header `deps/json.hpp`

### Build

The provided `Makefile` prefers an already-available Eigen install (for example
`/usr/include/eigen3` or `/usr/local/include/eigen3`) and only fetches Eigen
into `deps/eigen-3.4.0/` when no usable local copy is available:

```bash
make            # fetch Eigen if needed, then build static portable binary
make dynamic    # build dynamic-linked binary
make clean      # remove built binary (keeps fetched Eigen)
```

To fetch Eigen explicitly (e.g., to pre-cache it before an offline build):

```bash
make deps
```

To use a system Eigen install instead of the auto-fetched copy, point `EIGEN_DIR` at it:

```bash
make EIGEN_DIR=/usr/include/eigen3
```

Or pin a different Eigen release:

```bash
make EIGEN_VERSION=3.4.0
```

### Dependency provenance and licensing

- `deps/json.hpp` is the upstream [nlohmann/json](https://github.com/nlohmann/json)
  single-header release, version **3.11.3**.
- The vendored header retains the upstream SPDX metadata:
  `SPDX-FileCopyrightText: 2013-2023 Niels Lohmann`
  and `SPDX-License-Identifier: MIT`.
- Step 8 integration review keeps `json.hpp` vendored because it is a
  header-only dependency that keeps the iPSAE build self-contained and avoids
  adding a separate package-manager requirement for the production pipeline.
- Eigen is **not** vendored. The build uses Eigen **3.4.0** headers either from
  a local/system install or from the upstream tarball fetched by `make deps`.
  Eigen is distributed under the MPL-2.0 license.

### Manual build (without `make`)

```bash
# One-time: fetch Eigen
curl -fsSL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz \
    | tar -xz -C deps

# Compile (static linking, ~3MB portable binary)
g++ -O3 -mtune=generic -fopenmp -std=c++17 -static \
    -I deps/eigen-3.4.0 -I deps \
    ipsae_cpp.cpp -o ipsae_cpp
```
