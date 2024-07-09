# Dockerfile for building FANS
#
# How to use:
# Build image (for example with tag FANS_build_image):
#   docker build -t FANS_build_image .
# Run image and map FANS directory into the container:
#   docker run -it -v $FE_DIR:/FANS --name FANS_build_container FANS_build_image
# Inside the container you can then invoke cmake and start building the source, 
# as described in the README.md.

# during build time, don't ask for user input (has to be included in every stage
# to take effect
ARG DEBIAN_FRONTEND=noninteractive

################################################################################
FROM ubuntu:22.04 AS fans_base
ARG DEBIAN_FRONTEND

################################################################################
FROM fans_base AS fans_build
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
FROM fans_build AS fans_built
ARG DEBIAN_FRONTEND

# copy whole repository into the image
COPY . /FANS
RUN rm -rf /FANS/build

WORKDIR /FANS/build

RUN cmake .. && cmake --build . -j && cpack -G "DEB"

# only keep the runtime package in order to
RUN find ./packages/ -type f ! -name '*runtime*.deb' -delete

################################################################################
FROM fans_base AS fans
ARG DEBIAN_FRONTEND

COPY --from=fans_built /FANS/build/packages/ /FANS/build/packages/
RUN [ $(ls /FANS/build/packages/*.deb | wc -l) -eq 1 ] || (echo "Error: Expected exactly one .deb file" && exit 1)

RUN apt-get update -qq && apt-get install -y --no-install-recommends \
    /FANS/build/packages/*runtime*.deb \
    sudo \
    gosu \
    && apt-get clean \
    && apt-get autoremove --purge -y \
    && rm -rf /var/lib/apt/lists/*

RUN rm -rf /FANS

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

CMD ["FANS"]

################################################################################
FROM fans_build AS fans_dev
ARG DEBIAN_FRONTEND

RUN apt-get update -qq && apt-get install -y --no-install-recommends \
    sudo \
    gosu \
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
