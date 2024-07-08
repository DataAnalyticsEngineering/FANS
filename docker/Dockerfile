# Dockerfile for building FANS
#
# How to use:
# Build image (for example with tag FANS_build_image):
#   docker build -t FANS_build_image .
# Run image and map FANS directory into the container:
#   docker run -it -v $FE_DIR:/FANS --name FANS_build_container FANS_build_image
# Inside the container you can then invoke cmake and start building the source, 
# as described in the README.md.

################################################################################
FROM ubuntu:22.04 AS fans_runtime

# during build time, don't ask for user input
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update -qq && apt-get install -y --no-install-recommends \
# packages specifically required by FANS
    libhdf5-dev \
    libopenmpi-dev \
    libeigen3-dev \
    libfftw3-dev \
    libfftw3-mpi-dev \
    && rm -rf /var/lib/apt/lists/*

################################################################################
FROM fans_runtime AS fans_build_env

RUN apt-get update -qq && apt-get install -y --no-install-recommends \
# basic packages
    software-properties-common \
    git \
# build basics
    build-essential \
    cmake \
    file \
    && rm -rf /var/lib/apt/lists/*

################################################################################
FROM fans_build_env AS fans_build

# copy whole repository into the image
COPY . /FANS
RUN rm -rf /FANS/build

WORKDIR /FANS/build

RUN cmake .. && cmake --build . -j && cpack -G "DEB"

################################################################################
FROM fans_runtime AS fans_publish
ARG pkg_name="fans-fans_runtime_0.0.0_amd64.deb"
COPY --from=fans_build /FANS/build/${pkg_name} /FANS/build/${pkg_name}
RUN apt install /FANS/build/${pkg_name}

RUN apt-get update -qq && DEBIAN_FRONTEND=noninteractive apt-get install -y \
    sudo \
    gosu \
    && apt-get clean \
    && apt-get autoremove --purge -y \
    && rm -rf /var/lib/apt/lists/*

# Create a non-root user
ENV USER develop
ENV UID 1000
ENV GID 100
ENV HOME /home/${USER}
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
FROM fans_build_env AS fans_dev_env

RUN apt-get update -qq && DEBIAN_FRONTEND=noninteractive apt-get install -y \
    sudo \
    gosu \
    && apt-get clean \
    && apt-get autoremove --purge -y \
    && rm -rf /var/lib/apt/lists/*

# Create a non-root user
ENV USER develop
ENV UID 1000
ENV GID 100
ENV HOME /home/${USER}
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
