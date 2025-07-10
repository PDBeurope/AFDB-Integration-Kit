# 1. Use an official Python runtime as a parent image
FROM python:3.12-slim-bookworm

# 2. Install Node.js and npm
# We are using the official Node.js setup script for Node.js 20.x (LTS)
RUN apt-get update && \
    apt-get install -y curl && \
    curl -fsSL https://deb.nodesource.com/setup_20.x | bash - && \
    apt-get install -y nodejs && \
    rm -rf /var/lib/apt/lists/*

# 3. Set the working directory in the container
WORKDIR /app

# 4. Copy dependency definition files
# Copy Python dependency files first
COPY pyproject.toml requirements.txt uv.lock ./
# Install Python dependencies using uv
RUN pip install uv
RUN uv pip install --system --no-cache-dir -r requirements.txt

# 5. Copy the rest of the application code, including the molstar submodule
COPY . .

# 6. Install Node.js dependencies and build the molstar project
# The molstar project contains its own package.json
RUN cd molstar && npm install && npm run build

# 7. Define the entrypoint for the container
ENTRYPOINT ["python", "main.py"]

# 8. Set a default command (e.g., to show help)
CMD ["--help"]
