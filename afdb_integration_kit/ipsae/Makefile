# Makefile for ipSAE C++ Implementation
# Used by production_pipeline.py for automatic compilation

CXX ?= g++
CXXFLAGS = -O3 -mtune=generic -fopenmp -std=c++17
LDFLAGS = -fopenmp

# Eigen is a header-only third-party dependency and is NOT vendored in this
# repository. It is fetched on first build and verified against a pinned
# checksum. Packagers may point at a system install instead, e.g.
#   make EIGEN_DIR=/usr/include/eigen3
EIGEN_VERSION ?= 3.4.0
EIGEN_URL ?= https://gitlab.com/libeigen/eigen/-/archive/$(EIGEN_VERSION)/eigen-$(EIGEN_VERSION).tar.gz
EIGEN_SHA256 ?= 8586084f71f9bde545ee7fa6d00288b264a2b7ac3607b974e54d13e7162c1c72
EIGEN_DIR ?= deps/eigen-$(EIGEN_VERSION)

INCLUDES = -I $(EIGEN_DIR) -I deps

SRC = ipsae_cpp.cpp
TARGET = ipsae_cpp

# Static linking (creates portable binary ~3MB)
STATIC_FLAGS = -static

.PHONY: all deps clean check help

all: $(TARGET)

# Fetch and verify Eigen unless it is already present (vendored, fetched
# earlier, or supplied via EIGEN_DIR). Safe to run offline once satisfied.
deps:
	@if [ -d "$(EIGEN_DIR)" ]; then \
	    echo "Eigen present at $(EIGEN_DIR)"; \
	else \
	    echo "Fetching Eigen $(EIGEN_VERSION)"; \
	    mkdir -p deps; \
	    curl -fsSL -o deps/eigen-$(EIGEN_VERSION).tar.gz "$(EIGEN_URL)"; \
	    echo "$(EIGEN_SHA256)  deps/eigen-$(EIGEN_VERSION).tar.gz" | sha256sum -c - ; \
	    tar -xzf deps/eigen-$(EIGEN_VERSION).tar.gz -C deps; \
	    rm -f deps/eigen-$(EIGEN_VERSION).tar.gz; \
	    echo "Eigen unpacked to $(EIGEN_DIR)"; \
	fi

$(TARGET): $(SRC) | deps
	$(CXX) $(CXXFLAGS) $(STATIC_FLAGS) $(INCLUDES) $< -o $@ $(LDFLAGS)
	@echo "Built $(TARGET) (static linking - portable)"

dynamic: $(SRC) | deps
	$(CXX) $(CXXFLAGS) $(INCLUDES) $< -o $(TARGET) $(LDFLAGS)
	@echo "Built $(TARGET) (dynamic linking)"

# Verify dependencies are present. Does not create anything.
check:
	@echo "Checking dependencies..."
	@test -d "$(EIGEN_DIR)" || (echo "ERROR: $(EIGEN_DIR) not found - run 'make deps' or set EIGEN_DIR" && exit 1)
	@test -f deps/json.hpp || (echo "ERROR: deps/json.hpp not found" && exit 1)
	@echo "All dependencies found."

clean:
	rm -f $(TARGET)
	@echo "Cleaned build artifacts"

help:
	@echo "ipSAE C++ Makefile"
	@echo ""
	@echo "Targets:"
	@echo "  all     - Fetch deps if needed and build the static binary (default)"
	@echo "  deps    - Fetch and verify Eigen only (useful for offline/CI pre-caching)"
	@echo "  dynamic - Build with dynamic linking"
	@echo "  check   - Verify dependencies exist (creates nothing)"
	@echo "  clean   - Remove built binary"
	@echo ""
	@echo "Overrides:"
	@echo "  make EIGEN_DIR=/usr/include/eigen3   # use a system Eigen"
	@echo "  make EIGEN_VERSION=3.4.0             # pin a different release"
