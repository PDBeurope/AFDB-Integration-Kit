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

### Running the `preprocess` Command

To run the `preprocess` command, you need to provide the input and output file paths.

A critical aspect of using this Docker container is managing file access. The container runs in isolation from your local file system. To allow the script inside the container to read your input files and write the output, you must mount a local directory into the container using a Docker volume.

**Example:**

Let's say you have your input files (e.g., `file.cif`) in a local directory named `data`.

1.  **Create a `data` directory** in your project's root and place your input files there.

2.  **Run the container with a volume mount**:

    ```bash
    docker run --rm -v "$(pwd)/data:/app/data" af-toolkit preprocess -i /app/data/file.cif -o /app/data/file.bcif
    ```

    Let's break down this command:
    *   `docker run --rm`: Runs the container and automatically removes it when it exits.
    *   `-v "$(pwd)/data:/app/data"`: This is the volume mount.
        *   `$(pwd)/data`: This is the absolute path to your local `data` directory.
        *   `/app/data`: This is the path inside the container where your local directory will be accessible. `/app` is the working directory defined in the `Dockerfile`.
    *   `af-toolkit`: The name of the image to run.
    *   `preprocess -i /app/data/file.cif -o /app/data/file.bcif`: These are the arguments passed to the Python script. Note that the file paths are relative to the container's file system (`/app/data`), not your local machine.

The output file (`file.bcif`) will be created in the `/app/data` directory inside the container, and because of the volume mount, it will appear in your local `data` directory.
