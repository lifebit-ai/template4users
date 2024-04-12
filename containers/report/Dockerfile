FROM continuumio/miniconda3@sha256:a2e6aa4cd0b6dd696ae9e3e5732943250a977ab3a42b2fe5fb7ef0c19d2d9f16

LABEL name="quay.io/lifebitaiorg/report" \
      description="Docker image with all packages to make reports" \
      version="1.2.0" \
      maintainer=""

# Use a non-root user for security reasons
USER 1001

ARG ENV_NAME="base"

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        procps=2:3.3.16-1ubuntu2 && \
    rm -rf /var/lib/apt/lists/*

# Install mamba for faster installation in the subsequent step
RUN conda install -c conda-forge mamba -y

# Install the conda environment
COPY environment.yml /
RUN mamba env update --quiet --name ${ENV_NAME} --file /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/${ENV_NAME}/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN mamba env export --name ${ENV_NAME} > ${ENV_NAME}_exported.yml
