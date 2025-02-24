FROM ubuntu:18.04

LABEL org.opencontainers.image.authors="Sander Thierens"

ADD requirements.txt requirements.txt
ADD environment.yml environment.yml

RUN apt-get update
RUN apt-get install -y curl file
RUN curl -L micro.mamba.pm/install.sh | /bin/bash

# Set up environment variables for MicroMamba
ENV MAMBA_ROOT_PREFIX=/root/.micromamba
ENV PATH="/root/.local/bin:$MAMBA_ROOT_PREFIX/bin:$PATH"

# Update micromamba
RUN micromamba self-update

# Create mamba environemnt based on file
RUN micromamba env create -f environment.yml

# Ensure the environment is activated every time the container is used
SHELL ["micromamba", "run", "-n", "mini_ex", "/bin/bash", "-c"]

# Default command (optional, but useful for interactive usage)
ENTRYPOINT ["micromamba", "run", "-n", "mini_ex"]
CMD ["/bin/bash"]
