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

# Precompile Julia packages
RUN julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

# Copy the rest of the application
COPY . .

# Set default command to generate the document
CMD ["julia", "--project=.", "generate_lecture_doc.jl"]
