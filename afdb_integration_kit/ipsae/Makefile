# Makefile for ipSAE C++ Implementation
# Used by production_pipeline.py for automatic compilation

CXX ?= g++
CXXFLAGS = -O3 -mtune=generic -fopenmp -std=c++17
LDFLAGS = -fopenmp

# Eigen is fetched on demand into deps/ (not vendored in the repo).
# Override EIGEN_VERSION or EIGEN_URL to use a different release.
EIGEN_VERSION ?= 3.4.0
EIGEN_DIR = deps/eigen-$(EIGEN_VERSION)
EIGEN_URL ?= https://gitlab.com/libeigen/eigen/-/archive/$(EIGEN_VERSION)/eigen-$(EIGEN_VERSION).tar.gz

# Include paths for dependencies
INCLUDES = -I $(EIGEN_DIR) -I deps

# Source and target
SRC = ipsae_cpp.cpp
TARGET = ipsae_cpp

# Static linking (creates portable binary ~3MB)
STATIC_FLAGS = -static

.PHONY: all clean static check help deps

# Default target - static linking (portable binary ~3MB, no runtime dependencies)
all: $(TARGET)

$(TARGET): $(SRC) | deps
	$(CXX) $(CXXFLAGS) $(STATIC_FLAGS) $(INCLUDES) $< -o $@ $(LDFLAGS)
	@echo "Built $(TARGET) (static linking - portable)"

# Dynamic target - faster compile, requires libs at runtime
dynamic: $(SRC) | deps
	$(CXX) $(CXXFLAGS) $(INCLUDES) $< -o $(TARGET) $(LDFLAGS)
	@echo "Built $(TARGET) (dynamic linking)"

# Fetch Eigen header tree if it isn't already present.
deps: $(EIGEN_DIR)

$(EIGEN_DIR):
	@echo "Fetching Eigen $(EIGEN_VERSION) into $(EIGEN_DIR)..."
	@mkdir -p deps
	@curl -fsSL "$(EIGEN_URL)" -o deps/eigen-$(EIGEN_VERSION).tar.gz
	@tar -xzf deps/eigen-$(EIGEN_VERSION).tar.gz -C deps
	@rm -f deps/eigen-$(EIGEN_VERSION).tar.gz
	@test -d "$(EIGEN_DIR)" || (echo "ERROR: extraction did not produce $(EIGEN_DIR)" && exit 1)
	@echo "Eigen $(EIGEN_VERSION) ready."

# Check that dependencies are available (will not fetch).
check:
	@echo "Checking dependencies..."
	@test -d $(EIGEN_DIR) || (echo "ERROR: $(EIGEN_DIR) not found. Run 'make deps' to fetch it." && exit 1)
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
	@echo "  deps    - Fetch Eigen $(EIGEN_VERSION) into $(EIGEN_DIR)"
	@echo "  check   - Verify dependencies exist (does not fetch)"
	@echo "  clean   - Remove built binary"
	@echo "  help    - Show this help message"
	@echo ""
	@echo "Variables:"
	@echo "  EIGEN_VERSION (default: $(EIGEN_VERSION))"
	@echo "  EIGEN_URL     (default: $(EIGEN_URL))"
	@echo ""
	@echo "Usage:"
	@echo "  make          # Fetch Eigen if needed, then build static binary"
	@echo "  make dynamic  # Build dynamic binary"
	@echo "  make clean    # Remove binary"
