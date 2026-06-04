# Makefile for ipSAE C++ Implementation
# Used by production_pipeline.py for automatic compilation

CXX ?= g++
CXXFLAGS = -O3 -mtune=generic -fopenmp -std=c++17
LDFLAGS = -fopenmp

# Eigen is fetched on demand into deps/ (not vendored in the repo).
# Prefer an already-available local/system install before attempting a fetch.
# Override EIGEN_DIR, EIGEN_VERSION, or EIGEN_URL as needed.
EIGEN_VERSION ?= 3.4.0
EIGEN_FETCH_DIR = deps/eigen-$(EIGEN_VERSION)
SYSTEM_EIGEN_DIR := $(firstword $(foreach dir,/usr/include/eigen3 /usr/local/include/eigen3,$(if $(wildcard $(dir)/Eigen),$(dir),)))
EIGEN_DIR ?= $(if $(wildcard $(EIGEN_FETCH_DIR)/Eigen),$(EIGEN_FETCH_DIR),$(if $(SYSTEM_EIGEN_DIR),$(SYSTEM_EIGEN_DIR),$(EIGEN_FETCH_DIR)))
EIGEN_URL ?= https://gitlab.com/libeigen/eigen/-/archive/$(EIGEN_VERSION)/eigen-$(EIGEN_VERSION).tar.gz

# Include paths for dependencies
INCLUDES = -I $(EIGEN_DIR) -I deps

# Source and target
SRC = ipsae_cpp.cpp
TARGET = ipsae_cpp

# Static linking (creates portable binary ~3MB)
STATIC_FLAGS = -static

.PHONY: all clean static check help deps ensure-eigen

# Default target - static linking (portable binary ~3MB, no runtime dependencies)
all: $(TARGET)

$(TARGET): $(SRC) | ensure-eigen
	$(CXX) $(CXXFLAGS) $(STATIC_FLAGS) $(INCLUDES) $< -o $@ $(LDFLAGS)
	@echo "Built $(TARGET) (static linking - portable)"

# Dynamic target - faster compile, requires libs at runtime
dynamic: $(SRC) | ensure-eigen
	$(CXX) $(CXXFLAGS) $(INCLUDES) $< -o $(TARGET) $(LDFLAGS)
	@echo "Built $(TARGET) (dynamic linking)"

# Resolve Eigen from the selected directory, fetching only when the selected
# directory is the repo-local cache path.
ensure-eigen:
	@if [ -d "$(EIGEN_DIR)/Eigen" ]; then \
		echo "Using Eigen headers from $(EIGEN_DIR)"; \
	elif [ "$(EIGEN_DIR)" = "$(EIGEN_FETCH_DIR)" ]; then \
		$(MAKE) deps; \
	else \
		echo "ERROR: EIGEN_DIR=$(EIGEN_DIR) does not contain Eigen headers" >&2; \
		exit 1; \
	fi

# Fetch Eigen header tree into the repo-local cache path.
deps: $(EIGEN_FETCH_DIR)

$(EIGEN_FETCH_DIR):
	@echo "Fetching Eigen $(EIGEN_VERSION) into $(EIGEN_FETCH_DIR)..."
	@mkdir -p deps
	@curl -fsSL "$(EIGEN_URL)" -o deps/eigen-$(EIGEN_VERSION).tar.gz
	@tar -xzf deps/eigen-$(EIGEN_VERSION).tar.gz -C deps
	@rm -f deps/eigen-$(EIGEN_VERSION).tar.gz
	@test -d "$(EIGEN_FETCH_DIR)" || (echo "ERROR: extraction did not produce $(EIGEN_FETCH_DIR)" && exit 1)
	@echo "Eigen $(EIGEN_VERSION) ready."

# Check that dependencies are available (will not fetch).
check:
	@echo "Checking dependencies..."
	@test -d "$(EIGEN_DIR)/Eigen" || (echo "ERROR: Eigen headers not found under $(EIGEN_DIR). Set EIGEN_DIR or run 'make deps'." && exit 1)
	@test -f deps/json.hpp || (echo "ERROR: deps/json.hpp not found" && exit 1)
	@echo "All dependencies found."

# Clean build artifacts (does not delete fetched Eigen tree).
clean:
	rm -f $(TARGET)
	@echo "Cleaned build artifacts"

# Help target
help:
	@echo "ipSAE C++ Makefile"
	@echo ""
	@echo "Targets:"
	@echo "  all     - Build static portable binary (default; auto-fetches Eigen if missing)"
	@echo "  dynamic - Build with dynamic linking"
	@echo "  deps    - Fetch Eigen $(EIGEN_VERSION) into $(EIGEN_FETCH_DIR)"
	@echo "  check   - Verify the selected Eigen directory and json.hpp exist"
	@echo "  clean   - Remove built binary"
	@echo "  help    - Show this help message"
	@echo ""
	@echo "Variables:"
	@echo "  EIGEN_DIR     (default: $(EIGEN_DIR))"
	@echo "  EIGEN_VERSION (default: $(EIGEN_VERSION))"
	@echo "  EIGEN_URL     (default: $(EIGEN_URL))"
	@echo ""
	@echo "Usage:"
	@echo "  make          # Use system/repo-local Eigen, fetching only if needed"
	@echo "  make dynamic  # Build dynamic binary"
	@echo "  make clean    # Remove binary"
