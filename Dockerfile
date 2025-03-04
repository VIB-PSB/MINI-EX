FROM mambaorg/micromamba

LABEL org.opencontainers.image.authors="Camilla FERRARI"
ARG MAMBA_DOCKERFILE_ACTIVATE=1

ADD requirements.txt requirements.txt

USER root
RUN apt-get update && apt-get install -y file procps

USER mambauser

# Installing python 3.6 and graph-tll
RUN micromamba install -y python=3.6
RUN micromamba install -c conda-forge graph-tool=2.43

# Installing remaining packages
RUN pip install -r requirements.txt

ENV PATH="/opt/conda/bin:$PATH"
