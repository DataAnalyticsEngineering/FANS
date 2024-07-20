# Docker
We provide a set of docker images for different use cases:
- **fans**: Contains the minimum environment for FANS to run and has the package 'fans' installed. Offers a non-root user (`develop`) and handling of UID and GID to not mess up permissions when volume mounting into the container. Meant for users of FANS that can't install the fans package directly.
- **fans-ci**: Contains the minimum tools to build FANS (including dev packages of dependencies with the required headers), but does not include FANS itself. Meant for a CI workflow.
- **fans-dev**: Based upon fans-ci, but offers a non-root user (`develop`) and handling of UID and GID to not mess up permissions when volume mounting into the container. Meant for developers that can't install the required tools on their machines.

The images are built for both linux/amd64 and linux/arm64 and are available on our [Dockerhub profile](https://hub.docker.com/u/unistuttgartdae).

## Set up a Container
To set up a development container (same procedure for the fans image) with your current working directory mounted into it, type:
```bash
docker create --name fans-dev -it \
  -e HOST_UID=$(id -u) \
  -e HOST_GID=$(id -g) \
  -v /etc/localtime:/etc/localtime:ro \
  -v /etc/timezone:/etc/timezone:ro \
  -v $PWD/:/workspace/ \
  unistuttgartdae/fans-dev
```
The `-e` options provide the entrypoint script of the container with your host user ID and GID, such that the user ID and GID inside the container can be adapted to match yours. This is done to not mess up file permissions in the mounted volumes. The two volume mounts of `/etc/localtime` and `/etc/timezone` are required to have the host date and time inside the container.

To start the container and attach a shell, run:
```bash
docker start -i fans
```
As the `fans-dev` image is meant for developers that can't or don't want to install the required dependencies directly on their machine, the following workflow is suggested: You would work on the code as usual on your host; and only to build and run FANS you would attach to the container (see next section).

## FANS Usage Within a Container
Since FANS needs to run inside the container, you first need to start and attach to the container (`docker start -i fans-dev`). Then you are in interactive mode and can continue using FANS as described in the main [README](../README.md#usage).

However, if you need to use FANS within scripts, Docker's interactive mode (`-i`) is not suitable. In such cases you need to use `docker exec`. You basically replace the original `FANS` call by `docker exec -u develop -w /workspace/test fans-dev [original call]`. For example in conjunction with nohup:
```bash
docker start fans-dev
nohup /usr/bin/time -v docker exec -u develop -w /workspace/test fans-dev [original call] &
docker stop fans-dev
```