# Use a minimal micromamba base image (modern 2026 standard for bioinformatics)
FROM mambaorg/micromamba:latest

# 1. Set the working directory
WORKDIR /app
# 2. Copy EVERYTHING first (including pyproject.toml, src/, etc.)
# Using chown ensures micromamba has permissions to write temporary pip files
COPY --chown=$MAMBA_USER:$MAMBA_USER . /app

# 3. Install dependencies directly into the base environment
# This avoids the need for 'conda activate' inside the container
RUN micromamba install -y -n base -f environment.yml && \
    micromamba clean --all --yes

# Permanently set the library path to resolve the CXXABI_1.3.15 issue
# This replaces your manual 'export' fix for the container environment
ENV LD_LIBRARY_PATH="/opt/conda/lib:${LD_LIBRARY_PATH}"
ENV PATH="/opt/conda/bin:${PATH}"

# Default command to run when the container starts
ENTRYPOINT ["stmut-hires"]
CMD ["--help"]