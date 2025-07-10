# Stage 1: Build C++ tools
FROM ubuntu:24.04 AS builder

RUN apt-get update -y && apt-get install -y \
    build-essential \
    cmake \
    git \
    libeigen3-dev \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /local

# Build libcifpp
RUN git clone https://github.com/PDB-REDO/libcifpp.git
WORKDIR /local/libcifpp
RUN cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
RUN cmake --build build
RUN cmake --install build

# Build libmcfp
WORKDIR /local
RUN git clone https://github.com/mhekkel/libmcfp.git
WORKDIR /local/libmcfp
RUN cmake -S . -B build
RUN cmake --build build --config Release
RUN cmake --install build

# Build dssp
WORKDIR /local
RUN git clone https://github.com/PDB-REDO/dssp.git
WORKDIR /local/dssp
RUN cmake -S . -B build
RUN cmake --build build --config Release
RUN cmake --install build

# Stage 2: Main image
FROM ubuntu:24.04 AS runtime

# Install system dependencies
RUN apt-get update && \
    apt-get install -y \
        curl \
        openjdk-17-jre-headless \
        git \
        python3.12 \
        python3.12-venv \
        python3-pip \
        npm \
        nodejs \
        ca-certificates \
        && rm -rf /var/lib/apt/lists/*

# Set python3.12 as default python
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3.12 1

# Install Mol*
RUN npm install -g molstar

# Setup Nextflow
RUN mkdir -p /usr/local/bin && curl -s https://get.nextflow.io | bash -s -- && \
    mv nextflow /usr/local/bin/nextflow && chmod +x /usr/local/bin/nextflow

# Copy C++ tools from builder
COPY --from=builder /usr/local /usr/local

WORKDIR /app

# Copy Python dependency files
COPY pyproject.toml requirements.txt uv.lock ./
RUN pip install --break-system-packages uv
RUN uv pip install --system --no-cache-dir --break-system-packages -r requirements.txt

# Copy the rest of the application code
COPY . .

# Default command
CMD ["python", "main.py", "--help"]
