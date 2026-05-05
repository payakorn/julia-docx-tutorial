FROM julia:1.10

# Install required system dependencies: Python for CondaPkg/PythonCall
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app

# Copy Project.toml and Manifest.toml first to cache Julia package installation
COPY Project.toml Manifest.toml ./
COPY CondaPkg.toml ./

# Instantiate Julia packages without precompiling during the build
ENV JULIA_PKG_PRECOMPILE_AUTO=0
RUN julia --project=. -e 'using Pkg; Pkg.instantiate()'

# Copy the rest of the application
COPY . .

# Set default command to run the web server
CMD ["julia", "--project=.", "serv.jl", "9000"]
