# Docker
We provide a set of docker images for different use cases on our [Dockerhub profile](https://hub.docker.com/u/unistuttgartdae):
- **fans-ci**: Contains the minimum tools to build FANS (including dev packages of dependencies with the required headers), but does not include FANS itself. Meant for a CI workflow.
- **fans-dev**: Based upon fans-ci, but offers a non-root user (`develop`) and handling of UID and GID to not mess up permissions when volume mounting into the container. Meant as an quick to setup build environment for FANS.

Both images are built for linux/amd64 and linux/arm64 as well as for the three most recent Ubuntu LTS versions (focal, jammy, noble). The Ubuntu version can be selected through tags, e.g. `fans-dev:focal`; `noble` is equivalent to the `latest` tag. The architecture is selected automatically depending on your host platform.

## Set up a Container
To set up a development container with your current working directory (in there, use `git clone` to obtain the latest FANS version) mounted into it, type:
```bash
git clone https://github.tik.uni-stuttgart.de/DAE/FANS.git

docker create --name fans-dev -it \
  -e HOST_UID=$(id -u) \
  -e HOST_GID=$(id -g) \
  -v /etc/localtime:/etc/localtime:ro \
  -v /etc/timezone:/etc/timezone:ro \
  -v $PWD/:/workspace/ \
  unistuttgartdae/fans-dev
```
The `-e` options provide the entrypoint script of the container with your host user ID and GID, such that the user ID and GID inside the container can be adapted to match yours. This is done to not mess up file permissions in the mounted volumes. The two volume mounts of `/etc/localtime` and `/etc/timezone` are required to have the host date and time inside the container.

The following workflow is suggested: You would work on the code as usual on your host; and only to build and run FANS you would attach to the container:
```bash
docker start fans-dev
docker attach fans-dev
cd /workspace/FANS
mkdir build
cd build
cmake ..
cmake --build . -j
./FANS
../test/run.sh
```

## FANS Usage Within a Container
By building inside the container, FANS is linked against the container's libs and therefore must run inside the container. After attaching to the container you can then continue to use FANS as described in the main [README](../README.md#usage). Just remember that any input and output files need to visible to the container and thus must lie somewhere inside the mounted volumes.

Special care has to be taken if you need to use FANS within scripts, as Docker's interactive mode (`-i`) is not suitable in this case. Instead you need to use `docker exec`. One basically replaces the original `FANS` call by `docker exec -u develop -w /workspace/test fans-dev [original call]`. For example in conjunction with nohup:
```bash
docker start fans-dev
nohup /usr/bin/time -v docker exec -u develop -w /workspace/test fans-dev [original call] &
docker stop fans-dev
```