# Dockerfile for building FANS (see README.md for instructions)

# during build time, don't ask for user input (has to be included in every stage
# to take effect
ARG DEBIAN_FRONTEND=noninteractive

################################################################################
FROM ubuntu:noble AS fans_base
ARG DEBIAN_FRONTEND

################################################################################
FROM fans_base AS fans_ci
ARG DEBIAN_FRONTEND

RUN apt-get update -qq && apt-get install -y --no-install-recommends \
# basic packages
    software-properties-common \
    git \
# build basics
    build-essential \
    cmake \
    file \
# packages required to build FANS \
    libhdf5-dev \
    libopenmpi-dev \
    libeigen3-dev \
    libfftw3-dev \
    libfftw3-mpi-dev \
    && rm -rf /var/lib/apt/lists/*

################################################################################

FROM fans_ci AS fans_dev
ARG DEBIAN_FRONTEND

RUN apt-get update -qq && apt-get install -y --no-install-recommends \
    sudo \
    gosu \
    time \
    && apt-get clean \
    && apt-get autoremove --purge -y \
    && rm -rf /var/lib/apt/lists/*

# Create a non-root user
ENV USER=develop
ENV UID=1000
ENV GID=100
ENV HOME=/home/${USER}
RUN adduser --disabled-password \
    --gecos "Non-root user" \
    --uid ${UID} \
    --gid ${GID} \
    --home ${HOME} \
    ${USER} \
    #
    && echo ${USER} ALL=\(root\) NOPASSWD:ALL > /etc/sudoers.d/${USER} \
    && chmod 0440 /etc/sudoers.d/${USER}

# entrypoint script changes uid and gid to match given host uid and gid
COPY --chmod=755 docker/Dockerfile_user_env_entrypoint.sh /entrypoint.sh
ENTRYPOINT ["/entrypoint.sh"]

CMD ["bash"]
