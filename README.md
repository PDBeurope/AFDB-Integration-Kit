AFDB Integration Kit

```bash
docker build -t af-toolkit .
```

This command will create a Docker image named `af-toolkit` on your local machine.

## How to Run Commands

The Docker container is designed to be used as an executable. You can pass arguments to it directly on the `docker run` command line, which will be forwarded to the Python script inside the container.

### Displaying Help

To see the available commands and options, you can run the container without any additional arguments. This will execute the default command, which is to display the help message.

```bash
docker run --rm af-toolkit
```

This is equivalent to running `docker run --rm af-toolkit --help`.

### Running the `cif2bcif` Command

To run the `cif2bcif` command, you need to provide the input and output file paths.

A critical aspect of using this Docker container is managing file access. The container runs in isolation from your local file system. To allow the script inside the container to read your input files and write the output, you must mount a local directory into the container using a Docker volume.

**Example:**

Let's say you have your input files (e.g., `file.cif`) in a local directory named `data`.

1.  **Create a `data` directory** in your project's root and place your input files there.

2.  **Run the container with a volume mount**:

    ```bash
    docker run --rm -v "$(pwd)/data:/app/data" af-toolkit cif2bcif -i /app/data/file.cif -o /app/data/file.bcif
    ```

    Or to output a gzipped BCIF file:

    ```bash
    docker run --rm -v "$(pwd)/data:/app/data" af-toolkit cif2bcif -i /app/data/file.cif -o /app/data/file.bcif.gz
    ```

    Let's break down this command:
    *   `docker run --rm`: Runs the container and automatically removes it when it exits.
    *   `-v "$(pwd)/data:/app/data"`: This is the volume mount.
        *   `$(pwd)/data`: This is the absolute path to your local `data` directory.
        *   `/app/data`: This is the path inside the container where your local directory will be accessible. `/app` is the working directory defined in the `Dockerfile`.
    *   `af-toolkit`: The name of the image to run.
    *   `cif2bcif -i /app/data/file.cif -o /app/data/file.bcif`: These are the arguments passed to the Python script. Note that the file paths are relative to the container's file system (`/app/data`), not your local machine.

The output file (`file.bcif` or `file.bcif.gz`) will be created in the `/app/data` directory inside the container, and because of the volume mount, it will appear in your local `data` directory.

### Running the `batch_cif2bcif` Command

To convert all `.cif` files in a directory to `.bcif` or `.bcif.gz` in batch mode, use the `batch_cif2bcif` command. You can also specify the number of parallel workers and whether to output gzipped files.

**Example:**

```bash
docker run --rm -v "$(pwd)/data:/app/data" af-toolkit batch_cif2bcif -id /app/data -od /app/data/converted
```

To output gzipped BCIF files:

```bash
docker run --rm -v "$(pwd)/data:/app/data" af-toolkit batch_cif2bcif -id /app/data -od /app/data/converted --gzip
```

- `-id /app/data`: Input directory containing `.cif` files.
- `-od /app/data/converted`: Output directory for `.bcif` or `.bcif.gz` files.
- `--gzip`: (Optional) Output files as `.bcif.gz` instead of `.bcif`.
- `--workers N`: (Optional) Number of parallel workers (default: 4).
