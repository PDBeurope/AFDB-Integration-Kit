### Build ipSAE C++ Implementation

```bash
# Compile with static linking (creates portable binary)
g++ -O3 -march=native -fopenmp -std=c++17 -static \
    -I deps/eigen-3.4.0 -I deps \
    ipsae_cpp.cpp -o ipsae_cpp

# The resulting binary (~3MB) has no dependencies and can be copied/shared
chmod +x ipsae_cpp
```


