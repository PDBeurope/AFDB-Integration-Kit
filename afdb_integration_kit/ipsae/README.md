## ipSAE C++ implementation

### Dependencies

- A C++17 compiler with OpenMP (`g++` recommended)
- `curl` and `tar` (used by `make` to fetch Eigen on first build)
- [Eigen](https://eigen.tuxfamily.org/) **3.4.0** — fetched on demand into `deps/eigen-3.4.0/` (not vendored in this repo)
- [nlohmann/json](https://github.com/nlohmann/json) — vendored as the single-header `deps/json.hpp`

### Build

The provided `Makefile` will auto-fetch Eigen on first build:

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
