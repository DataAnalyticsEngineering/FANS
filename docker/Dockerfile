# During build time, don't ask for user input (has to be included in every stage
# to take effect)
ARG DEBIAN_FRONTEND=noninteractive
ARG UBUNTU_VERSION=noble
ARG USER=fans

################################################################################

FROM ubuntu:${UBUNTU_VERSION} AS fans_base
ARG DEBIAN_FRONTEND
ARG USER

# Context: https://askubuntu.com/questions/1513927/ubuntu-24-04-docker-images-now-includes-user-ubuntu-with-uid-gid-1000
RUN bash -c 'if id "ubuntu" &>/dev/null; then \
      touch /var/mail/ubuntu && \
      chown ubuntu /var/mail/ubuntu && \
      userdel -r ubuntu && \
      echo "Deleted user ubuntu."; \
    fi'

# Create a non-root user
RUN useradd -m -s /bin/bash ${USER}

################################################################################

FROM fans_base AS fans_ci
ARG DEBIAN_FRONTEND

RUN apt-get update -qq && apt-get install -y --no-install-recommends \
# Build basics
    software-properties-common \
    build-essential \
# CMake + git for FetchContent + file for CPack
    cmake \
    git \
    file \
# FANS dependencies \
    libhdf5-openmpi-dev \
    libopenmpi-dev \
    libeigen3-dev \
    libfftw3-dev \
    libfftw3-mpi-dev \
    nlohmann-json3-dev \
# Required for preCICE Micro Manager Python bindings
    python3-dev \
# Clean up
    && apt-get clean \
    && apt-get autoremove --purge -y \
    && rm -rf /var/lib/apt/lists/*

################################################################################

FROM fans_ci AS fans_dev
ARG DEBIAN_FRONTEND
ARG USER
ARG FANS_venv=FANS_venv

RUN apt-get update -qq && apt-get install -y --no-install-recommends \
# Packages required for setting up the non-root user
    sudo \
    gosu \
# Some additional packages for convenience
    time \
    htop \
    vim \
    python3-pip \
    python3-venv \
    python-is-python3 \
# Clean up
    && apt-get clean \
    && apt-get autoremove --purge -y \
    && rm -rf /var/lib/apt/lists/*

# Create a python venv for test/h52xdmf.py script
USER ${USER}

RUN python -m venv /home/${USER}/venvs/${FANS_venv} && \
    echo "\nsource /home/${USER}/venvs/${FANS_venv}/bin/activate\n" >> /home/${USER}/.bashrc && \
    . /home/${USER}/venvs/${FANS_venv}/bin/activate && \
    python -m pip install --no-cache-dir h5py lxml

USER root

# Add fans user to sudoers
RUN echo ${USER} ALL=\(root\) NOPASSWD:ALL > /etc/sudoers.d/${USER} \
    && chmod 440 /etc/sudoers.d/${USER}

# Entrypoint script changes UID and GID to match given host UID and GID
COPY --chmod=755 docker/Dockerfile_user_env_entrypoint.sh /entrypoint.sh
ENTRYPOINT ["/entrypoint.sh"]

CMD ["bash"]
