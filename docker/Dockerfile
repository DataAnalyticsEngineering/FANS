# Dockerfile for building FANS (see README.md for instructions)

# during build time, don't ask for user input (has to be included in every stage
# to take effect
ARG DEBIAN_FRONTEND=noninteractive
ARG UBUNTU_VERSION=noble

################################################################################

FROM ubuntu:${UBUNTU_VERSION} AS fans_base
ARG DEBIAN_FRONTEND

# context: https://askubuntu.com/questions/1513927/ubuntu-24-04-docker-images-now-includes-user-ubuntu-with-uid-gid-1000
RUN if id "ubuntu" &>/dev/null; then \
      touch /var/mail/ubuntu && \
      chown ubuntu /var/mail/ubuntu && \
      userdel -r ubuntu; \
    fi

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

################################################################################
# These stages are for the standalone FANS container
################################################################################

FROM fans_base AS fans
ARG DEBIAN_FRONTEND
ARG WORKSPACE=/workspace
ARG BUILD_DIR=build

COPY ${BUILD_DIR}/packages/fans_*.deb ${WORKSPACE}/

RUN [ $(ls ${WORKSPACE}/*.deb | wc -l) -eq 1 ] || (echo "Error: Expected exactly one .deb file" && exit 1)

RUN apt-get update -qq && apt-get install -y --no-install-recommends \
    ${WORKSPACE}/*.deb \
    sudo \
    gosu \
    time \
    && apt-get clean \
    && apt-get autoremove --purge -y \
    && rm -rf /var/lib/apt/lists/*

RUN rm -rf ${WORKSPACE}

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
